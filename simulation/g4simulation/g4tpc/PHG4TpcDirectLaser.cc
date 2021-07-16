#include "PHG4TpcDirectLaser.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4HitContainer.h>

#include <phool/getClass.h>

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
   */
  static constexpr double Ne_dEdx = 1.56;   // keV/cm
  static constexpr double CF4_dEdx = 7.00;  // keV/cm
  static constexpr double Ne_NTotal = 43;    // Number/cm
  static constexpr double CF4_NTotal = 100;  // Number/cm
  static constexpr double Tpc_NTot = 0.5 * Ne_NTotal + 0.5 * CF4_NTotal;
  static constexpr double Tpc_dEdx = 0.5 * Ne_dEdx + 0.5 * CF4_dEdx;
  static constexpr double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
  static constexpr double Tpc_ElectronsPerGeV = 1e6*Tpc_ElectronsPerKeV;
  
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

  /// load and check G4Hit node
  hitnodename = "G4HIT_" + detector;
  auto *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
  {
    std::cout << Name() << " Could not locate G4HIT node " << hitnodename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  /// setup lasers
  setupLasers();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________
int PHG4TpcDirectLaser::process_event(PHCompositeNode *topNode)
{
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
void PHG4TpcDirectLaser::SetPhiStepping(int n, float min,float max)
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
void PHG4TpcDirectLaser::SetThetaStepping(int n, float min,float max)
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
void PHG4TpcDirectLaser::setupLasers()
{
  // clear previous lasers
  m_lasers.clear();
  
  /// default position
  const TVector3 position_base( 60*cm, 0., 105.5*cm );
  
  /// add lasers
  for( int i = 0; i<8; ++i )
  {
    Laser laser;
    
    // set laser direction 
    // TODO: sort out sign conventions
    laser.m_direction = (i < 4) ? 1:-1;
    laser.m_position = position_base;
    
    // adjust z 
    laser.m_position.SetZ( position_base.z()*laser.m_direction );

    // rotate around z
    laser.m_phi = M_PI/2*i; 
    laser.m_position.RotateZ( laser.m_phi );

    // append
    m_lasers.push_back( laser );
  }
  
}

//_____________________________________________________________
void PHG4TpcDirectLaser::AimToThetaPhi(float theta, float phi)
{
  for( const auto& laser:m_lasers )
  { AppendLaserTrack(theta,phi,laser); }
}

//_____________________________________________________________
void PHG4TpcDirectLaser::AimToPatternStep(int n)
{

  n=n%nTotalSteps;//trim against overflows

  std::cout << "PHG4TpcDirectLaser::AimToPatternStep - step: " << n << std::endl;
  currentPatternStep=n;
  int phiStep=n%nThetaSteps;
  int thetaStep=n/nPhiSteps;
  AimToThetaPhi((maxTheta-minTheta)/(nThetaSteps*1.)*thetaStep,(maxPhi-minPhi)/(nPhiSteps*1.)*phiStep);
  return;
}

//_____________________________________________________________
TVector3 PHG4TpcDirectLaser::GetCmStrike(TVector3 start, TVector3 direction) const
{
  TVector3 ret;
  float end=halfwidth_CM;
  if (start.Z()<0) end=-halfwidth_CM;
  float dist=end-start.Z();
  float direction_scale=dist/direction.Z();
  ret=start+direction*direction_scale;

  return ret;
}

//_____________________________________________________________
TVector3 PHG4TpcDirectLaser::GetFieldcageStrike(TVector3 start, TVector3 direction) const
{
  TVector3 ofc_strike=GetCylinderStrike(start,direction,end_CM);
  TVector3 ifc_strike=GetCylinderStrike(start,direction,begin_CM);
  TVector3 no_strike(999,999,999);

  //measure which one occurs 'first' along the track, by dividing by the trajectory z.

  float ifc_dist=(ifc_strike.Z()-start.Z())/direction.Z();
  float ofc_dist=(ofc_strike.Z()-start.Z())/direction.Z();

  if (ifc_dist<0){
    if (ofc_dist>0) return ofc_strike;
    return no_strike;
  }

  //now ifc strike is guaranteed to be positive or zero
  if (ofc_dist<0){
    return ifc_strike;
  }
  //now they're both positive
  if (ifc_dist<ofc_dist){
    return ifc_strike;
  } else return ofc_strike;

}

//_____________________________________________________________
TVector3  PHG4TpcDirectLaser::GetCylinderStrike(TVector3 s, TVector3 v, float radius) const
{

  const float R2=square(radius);

  //Generalized Parameters for collision with cylinder of radius R:
  //from quadratic formula solutions of when a vector intersects a circle:
  const float a = square(v.x())+ square(v.y());
  const float b = 2*(v.x()*s.x()+v.y()*s.y());
  const float c = square(s.x()) + square(s.y())-R2;

  float rootterm=square(b)-4*a*c;

  //if a==0 then we are parallel and will have no solutions.
  //if the rootterm is negative, we will have no real roots -- we are outside the cylinder and pointing skew to the cylinder such that we never cross.
  float t1=-1,t2=-1; //this is the distance, in units of v, we must travel to find a collision
  if (rootterm >= 0 && a > 0) {
    //Find the (up to) two points where we collide with the cylinder:
    float sqrtterm=sqrt(rootterm);
     t1 = (-b+sqrtterm)/(2*a);
     t2 = (-b-sqrtterm)/(2*a);
  }

  //if either of the t's are nonzero, we have a collision.  the collision closest to the start (hence with the smallest t that is greater than zero) is the one that happens.
  float min_t=t1;
  if (t2>t1 && t2>0) min_t=t2;
  TVector3 ret=s+v*min_t;
  return ret;
}

//_____________________________________________________________
void PHG4TpcDirectLaser::AppendLaserTrack(float theta, float phi, const PHG4TpcDirectLaser::Laser& laser)
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
  TVector3 dir( 0, 0, -1.*direction );

  //adjust direction:
  dir.RotateY(theta*direction);
  dir.RotateZ(phi*direction);
   
  // also rotate by laser azimuth
  dir.RotateZ(laser.m_phi );

  // print
  // if( Verbosity() )
  { std::cout << "PHG4TpcDirectLaser::AppendLaserTrack - position: " << pos << " direction: " << dir << std::endl; }
  
  //find collision point
  TVector3 cm_strike=GetCmStrike(pos,dir);
  TVector3 fc_strike=GetFieldcageStrike(pos,dir);
  TVector3 strike=cm_strike;
  if( fc_strike.Z()!=999){
    if ((fc_strike.Z()-pos.Z())/dir.Z()<(cm_strike.Z()-pos.Z())/dir.Z()){
      strike=fc_strike;
    }
  }

  //find length
  TVector3 delta=(strike-pos);
  float fullLength=delta.Mag();
  int nHitSteps=fullLength/maxHitLength+1;

  TVector3 start=pos;
  TVector3 end=start;
  TVector3 step=dir*(maxHitLength/(dir.Mag()));
  float stepLength=0;
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
