#include "PHG4TpcDirectLaser.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4HitContainer.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>

#include <trackbase_historic/SvtxTrack_v2.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>

#include <cassert>

namespace
{
  // utility
  template <class T> inline constexpr T square(const T &x) { return x * x; }

  // unique detector id for all direct lasers
  static const int detId = PHG4HitDefs::get_volume_id( "PHG4TpcDirectLaser" );

  // hard-coded number of electrons created per cm of laser
  static constexpr double electrons_per_cm=300;

  /*
   * number of electrons deposited in gas per GeV for ionizing particle
   * it is needed to convert the number of electrons deposited by the laser into some equivalent energy loss,
   * which is what G4Hit expects
   * copied from PHG4TpcElectronDrift::SetDefaultParameters
   * TODO: should really use an independent set of parameters to get them from macro, rather than copy
   */
  static constexpr double Ne_dEdx = 1.56;   // keV/cm
  static constexpr double CF4_dEdx = 7.00;  // keV/cm
  static constexpr double Ne_NTotal = 43;    // Number/cm
  static constexpr double CF4_NTotal = 100;  // Number/cm
  static constexpr double Tpc_NTot = 0.5 * Ne_NTotal + 0.5 * CF4_NTotal;
  static constexpr double Tpc_dEdx = 0.5 * Ne_dEdx + 0.5 * CF4_dEdx;
  static constexpr double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
  static constexpr double Tpc_ElectronsPerGeV = 1e6*Tpc_ElectronsPerKeV;
  
  ///@name units
  //@{
  static constexpr double mm = 0.10;
  static constexpr double cm = 1.0;
  //@}

  /// length of generated G4Hits along laser track
  static constexpr double maxHitLength=1.0;//1cm.

  // inner and outer radii of field cages/TPC
  static constexpr double begin_CM = 20.*cm ;
  static constexpr double end_CM = 78.*cm;

  //half the thickness of the CM;
  static constexpr double halfwidth_CM = 0.5*cm;

  //_____________________________________________________________  
  TVector3 central_membrane_intersection(TVector3 start, TVector3 direction)
  {
    
    const double end = start.z() > 0 ? halfwidth_CM:-halfwidth_CM;
    const double dist=end-start.z();
    const double direction_scale=dist/direction.z();
    return start + direction * direction_scale;
  }
   
  //_____________________________________________________________
  std::pair<TVector3,bool>  cylinder_line_intersection(TVector3 s, TVector3 v, double radius)
  {
    
    const double R2=square(radius);
    
    //Generalized Parameters for collision with cylinder of radius R:
    //from quadratic formula solutions of when a vector intersects a circle:
    const double a = square(v.x())+ square(v.y());
    const double b = 2*(v.x()*s.x()+v.y()*s.y());
    const double c = square(s.x()) + square(s.y())-R2;
    
    const double rootterm=square(b)-4*a*c;
    
    /* 
     * if a==0 then we are parallel and will have no solutions.
     * if the rootterm is negative, we will have no real roots,
     * we are outside the cylinder and pointing skew to the cylinder such that we never cross.
     */
    if( rootterm < 0 || a == 0 ) return std::make_pair(TVector3(), false);
    
    //Find the (up to) two points where we collide with the cylinder:
    const double sqrtterm=std::sqrt(rootterm);
    const double t1 = (-b+sqrtterm)/(2*a);
    const double t2 = (-b-sqrtterm)/(2*a);
    
    /* 
    * if either of the t's are nonzero, we have a collision
    * the collision closest to the start (hence with the smallest t that is greater than zero) is the one that happens.
    */
    const double& min_t = (t2<t1 && t2>0) ? t2:t1;
    return std::make_pair(s+v*min_t, true );
  }
  
  //_____________________________________________________________
  std::pair<TVector3,bool> field_cage_intersection(TVector3 start, TVector3 direction)
  {
    const auto ofc_strike=cylinder_line_intersection(start,direction,end_CM);
    const auto ifc_strike=cylinder_line_intersection(start,direction,begin_CM);
    
    // if either of the two intersection is invalid, return the other
    if( !ifc_strike.second ) return ofc_strike;
    if( !ofc_strike.second ) return ifc_strike;
    
    // both intersection are valid, calculate signed distance to start z
    const auto ifc_dist=(ifc_strike.first.Z()-start.Z())/direction.Z();
    const auto ofc_dist=(ofc_strike.first.Z()-start.Z())/direction.Z();
   
    if(ifc_dist<0) return (ofc_dist > 0) ? ofc_strike:std::make_pair( TVector3(), false );
    else if( ofc_dist<0 ) return ifc_strike;
    else return (ifc_dist<ofc_dist) ? ifc_strike:ofc_strike;
    
  }
  
  /// TVector3 stream
  inline std::ostream& operator << (std::ostream& out, const TVector3& vector )
  {
    out << "( " << vector.x() << ", " << vector.y() << ", " << vector.z() << ")";
    return out;
  }

}

//_____________________________________________________________
PHG4TpcDirectLaser::PHG4TpcDirectLaser(const std::string &name)
  : SubsysReco(name)
{}

//_____________________________________________________________
int PHG4TpcDirectLaser::InitRun(PHCompositeNode *topNode)
{

  // load and check G4Hit node
  hitnodename = "G4HIT_" + detector;
  auto *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
  {
    std::cout << Name() << " Could not locate G4HIT node " << hitnodename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // find or create track map
  /* it is used to store laser parameters on a per event basis */
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, m_track_map_name );
  if( !m_track_map )
  {
    // find DST node and check
    PHNodeIterator iter(topNode);
    auto dstNode = static_cast<PHCompositeNode*>(iter.findFirst( "PHCompositeNode", "DST"));
    if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    // find or create SVTX node
    iter = PHNodeIterator(dstNode);
    auto node = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode",  "SVTX"));
    if( !node ) dstNode->addNode( node = new PHCompositeNode("SVTX") );

    // add track node
    m_track_map = new SvtxTrackMap_v1;
    node->addNode( new PHIODataNode<PHObject>(m_track_map, m_track_map_name, "PHObject") );
  }

  // setup lasers
  SetupLasers();

  // print configuration
  std::cout << "PHG4TpcDirectLaser::InitRun - m_autoAdvanceDirectLaser: " << m_autoAdvanceDirectLaser << std::endl;
  std::cout << "PHG4TpcDirectLaser::InitRun - phi steps: " << nPhiSteps << " min: " << minPhi << " max: " << maxPhi << std::endl;
  std::cout << "PHG4TpcDirectLaser::InitRun - theta steps: " << nThetaSteps << " min: " << minTheta << " max: " << maxTheta << std::endl;
  std::cout << "PHG4TpcDirectLaser::InitRun - nTotalSteps: " << nTotalSteps << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________
int PHG4TpcDirectLaser::process_event(PHCompositeNode *topNode)
{
  // load track map
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, m_track_map_name );
  assert( m_track_map );

  // load g4hit container
  m_g4hitcontainer = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if( !m_g4hitcontainer )
  {
    std::cout << PHWHERE << "Could not locate g4 hit node " << hitnodename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if( m_autoAdvanceDirectLaser )
  {

    AimToNextPatternStep();

  } else {

    // use arbitrary direction
    AimToThetaPhi( M_PI/180.*10., M_PI/180.*90 );

  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________
void PHG4TpcDirectLaser::SetPhiStepping(int n, double min,double max)
{
  if (n<0 || max<min)
  {
    std::cout << PHWHERE << " - invalid" << std::endl;
    return;
  }
  nPhiSteps=n;
  minPhi=min;
  maxPhi=max;
  nTotalSteps=nThetaSteps*nPhiSteps;
  return;
}
//_____________________________________________________________
void PHG4TpcDirectLaser::SetThetaStepping(int n, double min,double max)
{
  if (n<0 || max<min)
  {
    std::cout << PHWHERE << " - invalid" << std::endl;
    return;
  }
  nThetaSteps=n;
  minTheta=min;
  maxTheta=max;
  nTotalSteps=nThetaSteps*nPhiSteps;

  return;
}

//_____________________________________________________________
void PHG4TpcDirectLaser::SetupLasers()
{
  // clear previous lasers
  m_lasers.clear();

  // default position
  const TVector3 position_base( 60*cm, 0., 105.5*cm );

  // add lasers
  for( int i = 0; i<8; ++i )
  {
    Laser laser;

    // set laser direction
    /* 
     * first four lasers are on positive z readout plane, and shoot towards negative z
     * next four lasers are on negative z readout plane and shoot towards positive z 
     */
    laser.m_position = position_base;
    if( i < 4 )
    {

      laser.m_position.SetZ( position_base.z() );
      laser.m_direction = -1;

    } else {

      laser.m_position.SetZ( -position_base.z() );
      laser.m_direction = 1;

    }
    
    // rotate around z
    laser.m_phi = M_PI/2*i;
    laser.m_position.RotateZ( laser.m_phi );

    // append
    m_lasers.push_back( laser );
  }

}

//_____________________________________________________________
void PHG4TpcDirectLaser::AimToNextPatternStep()
{
  if( nTotalSteps>=1 )
  {
    AimToPatternStep(currentPatternStep);
    ++currentPatternStep;
  }
}

//_____________________________________________________________
void PHG4TpcDirectLaser::AimToThetaPhi(double theta, double phi)
{
  AppendLaserTrack(theta,phi,m_lasers[0]);
  
//   for( const auto& laser:m_lasers )
//   { if( laser.m_direction > 0 ) AppendLaserTrack(theta,phi,laser); }
}

//_____________________________________________________________
void PHG4TpcDirectLaser::AimToPatternStep(int n)
{
  //trim against overflows
  n=n%nTotalSteps;

  if( Verbosity() )
  { std::cout << "PHG4TpcDirectLaser::AimToPatternStep - step: " << n << "/" << nTotalSteps << std::endl; }

  // store as current pattern
  currentPatternStep=n;

  // calculate theta
  const int thetaStep = n/nPhiSteps;
  const double theta = minTheta + thetaStep*(maxTheta-minTheta)/nThetaSteps;

  // calculate phi
  const int phiStep = n%nPhiSteps;
  const double phi = minPhi + phiStep*(maxPhi-minPhi)/nPhiSteps;

  // generate laser tracks
  AimToThetaPhi(theta, phi );

  return;
}

//_____________________________________________________________
void PHG4TpcDirectLaser::AppendLaserTrack(double theta, double phi, const PHG4TpcDirectLaser::Laser& laser)
{

  if( !m_g4hitcontainer )
  {
    std::cout << PHWHERE << "invalid g4hit container. aborting" << std::endl;
    return;
  }

  // store laser position
  const auto& pos = laser.m_position;

  // define track direction
  const auto& direction = laser.m_direction;
  TVector3 dir( 0, 0, direction );

  //adjust direction
  dir.RotateY(theta*direction);
  dir.RotateZ(phi);

  // also rotate by laser azimuth
  dir.RotateZ(laser.m_phi );

  // print
  if( Verbosity() )
  { std::cout << "PHG4TpcDirectLaser::AppendLaserTrack - position: " << pos << " direction: " << dir << std::endl; }

  // store in SvtxTrack map
  if( m_track_map )
  {
    SvtxTrack_v2 track;
    track.set_x( pos.x() );
    track.set_y( pos.y() );
    track.set_z( pos.z() );

    // total momentum is irrelevant. What matters is the direction
    static constexpr double total_momentum = 1;
    track.set_px( total_momentum*dir.x() );
    track.set_py( total_momentum*dir.y() );
    track.set_pz( total_momentum*dir.z() );

    // insert in map
    m_track_map->insert( &track );
    
    if( Verbosity() )
    { std::cout << "PHG4TpcDirectLaser::AppendLaserTrack - position: " << pos << " direction: " << dir << std::endl; }

  }

  // find collision point
  // intersection to central membrane
  const auto cm_strike=central_membrane_intersection(pos,dir);
  
  // field cage intersection
  const auto fc_strike=field_cage_intersection(pos,dir);

  const auto& strike = ( fc_strike.second && (fc_strike.first.z()-pos.z())/dir.z()<(cm_strike.z()-pos.z())/dir.z()) ?
    fc_strike.first:
    cm_strike;

  //find length
  TVector3 delta=(strike-pos);
  double fullLength=delta.Mag();
  int nHitSteps=fullLength/maxHitLength+1;

  TVector3 start=pos;
  TVector3 end=start;
  TVector3 step=dir*(maxHitLength/(dir.Mag()));
  double stepLength=0;
  for (int i=0;i<nHitSteps;i++)
  {
    start=end;//new starting point is the previous ending point.
    if (i+1==nHitSteps){
      //last step is the remainder size
      end=strike;
      delta=start-end;
      stepLength=delta.Mag();
    } else {
      //all other steps are uniform length
      end=start+step;
      stepLength=step.Mag();
    }

    //from phg4tpcsteppingaction.cc
    auto hit = new PHG4Hitv1();
    hit->set_layer(99); // dummy number

    //here we set the entrance values in cm
    hit->set_x(0, start.X() / cm);
    hit->set_y(0, start.Y() / cm);
    hit->set_z(0, start.Z() / cm);
    //and the exist values
    hit->set_x(1, end.X() / cm);
    hit->set_y(1, end.Y() / cm);
    hit->set_z(1, end.Z() / cm);

    // momentum
    hit->set_px(0, dir.X()); // GeV
    hit->set_py(0, dir.Y());
    hit->set_pz(0, dir.Z());

    // time in ns
    hit->set_t(0, 0.0); // nanosecond
    hit->set_trkid(-1); // dummy number

    hit->set_px(1, dir.X());
    hit->set_py(1, dir.Y());
    hit->set_pz(1, dir.Z());

    hit->set_t(1, 0.0); // dummy number, nanosecond

    const double totalE=electrons_per_cm*stepLength/Tpc_ElectronsPerGeV;

    hit->set_eion(totalE);
    hit->set_edep(totalE);
    m_g4hitcontainer->AddHit( detId, hit );

  }

  return;
}
