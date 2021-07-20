/**
 * \file TpcDirectLaserReconstruction.cc
 * \brief performs the reconstruction of TPC direct laser tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcDirectLaserReconstruction.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState_v1.h>

#include <TVector3.h>

#include <cassert>


namespace
{

  //! range adaptor to be able to use range-based for loop
  template<class T> class range_adaptor
  {
    public:
    range_adaptor( const T& range ):m_range(range){}
    inline const typename T::first_type& begin() {return m_range.first;}
    inline const typename T::second_type& end() {return m_range.second;}
    private:
    T m_range;
  };
  
  //! convenience square method
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }

  //! get radius from x and y
  template<class T>
    inline constexpr T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }
    
  //_____________________________________________________________
  std::pair<TVector3,bool> cylinder_line_intersection(TVector3 s, TVector3 v, double radius)
  {
    
    const double R2=square(radius);
    
    //Generalized Parameters for collision with cylinder of radius R:
    //from quadratic formula solutions of when a vector intersects a circle:
    const double a = square(v.x())+ square(v.y());
    const double b = 2*(v.x()*s.x()+v.y()*s.y());
    const double c = square(s.x()) + square(s.y())-R2;
    const double rootterm=square(b)-4*a*c;
    
    //if a==0 then we are parallel and will have no solutions.
    //if the rootterm is negative, we will have no real roots -- we are outside the cylinder and pointing skew to the cylinder such that we never cross.
    if( rootterm <0 || a == 0 ) return std::make_pair( TVector3(), false );
    
    //Find the (up to) two points where we collide with the cylinder:
    const double sqrtterm=std::sqrt(rootterm);
    const double t1 = (-b+sqrtterm)/(2*a);
    const double t2 = (-b-sqrtterm)/(2*a);
    
    // chose smaller value of t, in absolute value
    const double& min_t = (t2<t1 && t2>0) ? t2:t1;
    return std::make_pair( s+v*min_t, true );
  }
  
  /// TVector3 stream
  inline std::ostream& operator << (std::ostream& out, const TVector3& vector )
  {
    out << "( " << vector.x() << ", " << vector.y() << ", " << vector.z() << ")";
    return out;
  }

}

//_____________________________________________________________________
TpcDirectLaserReconstruction::TpcDirectLaserReconstruction( const std::string& name ):
  SubsysReco( name)
  , PHParameterInterface(name)
{ 
  InitializeParameters(); 
}

//_____________________________________________________________________
int TpcDirectLaserReconstruction::Init(PHCompositeNode* topNode )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int TpcDirectLaserReconstruction::InitRun(PHCompositeNode* )
{

  // load parameters
  UpdateParametersWithMacro(); 
  return Fun4AllReturnCodes::EVENT_OK;
}


//_____________________________________________________________________
int TpcDirectLaserReconstruction::process_event(PHCompositeNode* topNode)
{
  // load nodes
  const auto res = load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  process_tracks();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcDirectLaserReconstruction::End(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//___________________________________________________________________________
void TpcDirectLaserReconstruction::SetDefaultParameters()
{}

//_____________________________________________________________________
int TpcDirectLaserReconstruction::load_nodes( PHCompositeNode* topNode )
{
  // get necessary nodes
  // tracks
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  assert(m_track_map);

  // hitset container
  m_hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");

  // clusters  
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  assert(m_cluster_map);
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void TpcDirectLaserReconstruction::process_tracks()
{
  if( !( m_track_map && m_cluster_map ) ) return;
  for( auto iter = m_track_map->begin(); iter != m_track_map->end(); ++iter )
  { process_track( iter->second ); }
}


//_____________________________________________________________________
void TpcDirectLaserReconstruction::process_track( SvtxTrack* track )
{
    
  // get track parameters
  const TVector3 origin( track->get_x(), track->get_y(), track->get_z() );
  const TVector3 direction( track->get_px(), track->get_py(), track->get_pz() );

  if( Verbosity() )
  { std::cout << "TpcDirectLaserReconstruction::process_track - position: " << origin << " direction: " << direction << std::endl; }
  
  // loop over hitsets
  for( const auto& [hitsetkey,hitset]:range_adaptor(m_hitsetcontainer->getHitSets()))
  {

    // only check TPC hitsets
    if( TrkrDefs::getTrkrId( hitsetkey ) != TrkrDefs::tpcId ) continue;
    
    // get corresponding clusters
    for( const auto& [key,cluster]:range_adaptor(m_cluster_map->getClusters(hitsetkey)))
    {
      // get cluster radius
      const auto cluster_r = get_r( cluster->getX(), cluster->getY() );      
      const auto [intersection,valid] = cylinder_line_intersection( origin, direction, cluster_r );
      if( !valid ) continue;
      
      // path length
      const auto pathlength = (intersection - origin).Mag();
    
      // create relevant state vector
      SvtxTrackState_v1 state( pathlength );
      state.set_x( intersection.x() );
      state.set_y( intersection.y() );
      state.set_z( intersection.z() );
      
      state.set_px( direction.x());
      state.set_py( direction.y());
      state.set_pz( direction.z());
      track->insert_state( &state );
    
      // associate cluster to track
      track->insert_cluster_key( key );      
    }
  }
  
}
