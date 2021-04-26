/*!
 * \file PHG4MicromegasHitReco.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHG4MicromegasHitReco.h"

#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <phparameter/PHParameterInterface.h>       // for PHParameterInterface

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                     // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHObject.h>                         // for PHObject


#include <TVector3.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>                            // for gsl_rng_alloc

#include <cassert>
#include <cmath>                                   // for atan2, sqrt, M_PI
#include <cstdlib>                                 // for exit
#include <iostream>                                 // for operator<<, basic...
#include <map>                                      // for _Rb_tree_const_it...
#include <numeric>

namespace
{
  //! convenient square function
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }

  //! return normalized gaussian centered on zero and of width sigma
  template<class T>
    inline T gaus( const T& x, const T& sigma )
  { return std::exp( -square(x/sigma)/2 )/(sigma*std::sqrt(2*M_PI)); }

  //! bind angle to [-M_PI,+M_PI[. This is useful to avoid edge effects when making the difference between two angles
  template<class T>
    inline T bind_angle( const T& angle )
  {
    if( angle >= M_PI ) return angle - 2*M_PI;
    else if( angle < -M_PI ) return angle + 2*M_PI;
    else return angle;
  }

  // local version of std::clamp, which is only available for c++17
  template<class T>
    constexpr const T& clamp( const T& v, const T& lo, const T& hi )
  { return (v < lo) ? lo : (hi < v) ? hi : v; }

  // this corresponds to integrating a gaussian centered on zero and of width sigma from xloc - pitch/2 to xloc+pitch/2
  template<class T>
    inline T get_rectangular_fraction( const T& xloc, const T& sigma, const T& pitch )
  { return (std::erf( (xloc + pitch/2)/(M_SQRT2*sigma) ) - std::erf( (xloc - pitch/2)/(M_SQRT2*sigma) ))/2; }

  /*
  this corresponds to integrating a gaussian centered on zero and of width sigma
  convoluted with a zig-zag strip response function, which is triangular from xloc-pitch to xloc+pitch, with a maximum of 1 at xloc
  */
  template<class T>
    inline T get_zigzag_fraction( const T& xloc, const T& sigma, const T& pitch )
  {
    return
      // rising edge
      (pitch - xloc)*(std::erf(xloc/(M_SQRT2*sigma)) - std::erf((xloc-pitch)/(M_SQRT2*sigma)))/(pitch*2)
      + (gaus(xloc-pitch, sigma) - gaus(xloc, sigma))*square(sigma)/pitch

      // descending edge
      + (pitch + xloc)*(std::erf((xloc+pitch)/(M_SQRT2*sigma)) - std::erf(xloc/(M_SQRT2*sigma)))/(pitch*2)
      + (gaus(xloc+pitch, sigma) - gaus(xloc, sigma))*square(sigma)/pitch;
  }

}

//___________________________________________________________________________
PHG4MicromegasHitReco::PHG4MicromegasHitReco(const std::string &name, const std::string& detector)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , m_detector(detector)
{
  // initialize rng
  const uint seed = PHRandomSeed();
  m_rng.reset( gsl_rng_alloc(gsl_rng_mt19937) );
  gsl_rng_set( m_rng.get(), seed );

  InitializeParameters();
}

//___________________________________________________________________________
int PHG4MicromegasHitReco::InitRun(PHCompositeNode *topNode)
{

  UpdateParametersWithMacro();

  // load parameters
  m_tmin = get_double_param("micromegas_tmin" );
  m_tmax = get_double_param("micromegas_tmax" );
  m_electrons_per_gev = get_double_param("micromegas_electrons_per_gev" );
  m_gain = get_double_param("micromegas_gain");
  m_cloud_sigma = get_double_param("micromegas_cloud_sigma");
  m_diffusion_trans = get_double_param("micromegas_diffusion_trans");
  m_zigzag_strips = get_int_param("micromegas_zigzag_strips");

  // printout
  std::cout
    << "PHG4MicromegasHitReco::InitRun\n"
    << " m_tmin: " << m_tmin << "ns, m_tmax: " << m_tmax << "ns\n"
    << " m_electrons_per_gev: " << m_electrons_per_gev << "\n"
    << " m_gain: " << m_gain << "\n"
    << " m_cloud_sigma: " << m_cloud_sigma << "cm\n"
    << " m_diffusion_trans: " << m_diffusion_trans << "cm/sqrt(cm)\n"
    << " m_zigzag_strips: " << std::boolalpha << m_zigzag_strips << "\n"
    << std::endl;

  // setup tiles
  setup_tiles( topNode );

  // get dst node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  // create hitset container if needed
  auto hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitsetcontainer)
  {
    std::cout << PHWHERE << "creating TRKR_HITSET." << std::endl;

    // find or create TRKR node
    PHNodeIterator dstiter(dstNode);
    auto trkrnode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!trkrnode)
    {
      trkrnode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrnode);
    }

    // create container and add to the tree
    hitsetcontainer = new TrkrHitSetContainer();
    auto newNode = new PHIODataNode<PHObject>(hitsetcontainer, "TRKR_HITSET", "PHObject");
    trkrnode->addNode(newNode);
  }

  // create hit truth association if needed
  auto hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!hittruthassoc)
  {
    std::cout << PHWHERE << "creating TRKR_HITTRUTHASSOC." << std::endl;

    // find or create TRKR node
    PHNodeIterator dstiter(dstNode);
    auto trkrnode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!trkrnode)
    {
      trkrnode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrnode);
    }

    hittruthassoc = new TrkrHitTruthAssoc();
    auto newNode = new PHIODataNode<PHObject>(hittruthassoc, "TRKR_HITTRUTHASSOC", "PHObject");
    trkrnode->addNode(newNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
int PHG4MicromegasHitReco::process_event(PHCompositeNode *topNode)
{
  // load relevant nodes
  // G4Hits
  const std::string g4hitnodename = "G4HIT_" + m_detector;
  auto g4hitcontainer = findNode::getClass<PHG4HitContainer>(topNode, g4hitnodename);
  assert(g4hitcontainer);

  // geometry
  const auto geonodename = full_geonodename();
  auto geonode =  findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
  assert(geonode);

  // Get the TrkrHitSetContainer node
  auto trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(trkrhitsetcontainer);

  // Get the TrkrHitTruthAssoc node
  auto hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  assert(hittruthassoc);

  // loop over layers in the g4hit container
  auto layer_range = g4hitcontainer->getLayers();
  for( auto layer_it = layer_range.first; layer_it != layer_range.second; ++layer_it )
  {
    // get layer
    const auto layer = *layer_it;

    // get relevant geometry
    auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(geonode->GetLayerGeom(layer));
    assert( layergeom );

    /* 
     * get the radius of the detector mesh. It depends on the drift direction 
     * it is used to calculate the drift distance of the primary electrons, and the
     * corresponding transverse diffusion
     */
    const auto mesh_radius = layergeom->get_drift_direction() == MicromegasDefs::DriftDirection::OUTWARD ?
      (layergeom->get_radius() + layergeom->get_thickness()/2): 
      (layergeom->get_radius() - layergeom->get_thickness()/2);
        
    // get corresponding hits
    const PHG4HitContainer::ConstRange g4hit_range = g4hitcontainer->getHits(layer);

    // loop over hits
    for( auto g4hit_it = g4hit_range.first; g4hit_it != g4hit_range.second; ++g4hit_it )
    {

      // get hit
      PHG4Hit* g4hit = g4hit_it->second;

      // check time window
      if(g4hit->get_t(0) > m_tmax) continue;
      if(g4hit->get_t(1) < m_tmin) continue;

      // get world coordinates
      TVector3 world_in( g4hit->get_x(0), g4hit->get_y(0), g4hit->get_z(0) );
      TVector3 world_out( g4hit->get_x(1), g4hit->get_y(1), g4hit->get_z(1) );

      // make sure that the mid point is in one of the tiles
      /*
       * at this point we do not check the strip validity.
       * This will be done when actually distributing electrons along the G4Hit track segment
       */
      const auto world_mid = (world_in+world_out)*0.5;
      const int tileid = layergeom->find_tile( world_mid );
      if( tileid < 0 ) continue;

      // number of primary elections
      const auto nprimary = get_primary_electrons( g4hit );
      if( !nprimary ) continue;

      // create hitset
      const TrkrDefs::hitsetkey hitsetkey = MicromegasDefs::genHitSetKey( layer, layergeom->get_segmentation_type(), tileid );
      const auto hitset_it = trkrhitsetcontainer->findOrAddHitSet(hitsetkey);

      // keep track of all charges
      using charge_map_t = std::map<int,double>;
      charge_map_t total_charges;

      // loop over primaries
      for( uint ie = 0; ie < nprimary; ++ie )
      {
        // put the electron at a random position along the g4hit path
        const auto t = gsl_ran_flat(m_rng.get(), 0.0, 1.0);
        auto world =  world_in*t + world_out*(1.0-t);
        
        if( m_diffusion_trans > 0 )
        {
          // add transeverse diffusion
          // first convert to polar coordinates
          const double radius = std::sqrt(square(world.x())+square(world.y()));
          const double phi = std::atan2(world.y(),world.x());
          const double drift_distance = std::abs(radius - mesh_radius);
          const double diffusion = gsl_ran_gaussian(m_rng.get(), m_diffusion_trans*std::sqrt(drift_distance));
          const double diffusion_angle = gsl_ran_flat(m_rng.get(), -M_PI, M_PI);
                    
          /*
           * diffusion happens in the phi,z plane (plane perpendicular to the radius direction)
           * with a magnitude 'diffusion' and an angle 'diffusion angle'
           * rotate back to cartesian coordinates
           */
          const auto sphi = std::sin(phi);
          const auto cphi = std::cos(phi);
          const auto salpha = std::sin(diffusion_angle);
          const auto calpha = std::cos(diffusion_angle);
          world += TVector3(-sphi*calpha*diffusion, cphi*calpha*diffusion, salpha*diffusion); 
        }
        
        // distribute charge among adjacent strips
        const auto fractions = distribute_charge( layergeom, tileid, world, m_cloud_sigma );

        // make sure fractions adds up to unity
        if( Verbosity() > 0 )
        {
          const auto sum = std::accumulate( fractions.begin(), fractions.end(), double( 0 ),
            []( double value, const charge_pair_t& pair ) { return value + pair.second; } );
          std::cout << "PHG4MicromegasHitReco::process_event - sum: " << sum << std::endl;
        }

        // generate gain for this electron
        const auto gain = get_single_electron_amplification();

        // merge to total charges
        for( const auto& pair: fractions )
        {
          const int strip = pair.first;
          if( strip < 0 || strip >= (int) layergeom->get_strip_count( tileid ) ) continue;

          const auto it = total_charges.lower_bound( strip );
          if( it != total_charges.end() && it->first == strip ) it->second += pair.second*gain;
          else total_charges.insert( it, std::make_pair( strip, pair.second*gain ) );
        }

      }

      // generate the key for this hit
      // loop over strips in list
      for( const auto pair:total_charges )
      {
        // get strip and bound check
        const int strip = pair.first;

        // get hit from hitset
        TrkrDefs::hitkey hitkey = MicromegasDefs::genHitKey(strip);
        auto hit = hitset_it->second->getHit(hitkey);
        if( !hit )
        {
          // create hit and insert in hitset
          hit = new TrkrHitv2;
          hitset_it->second->addHitSpecificKey(hitkey, hit);
        }

        // add energy from g4hit
        hit->addEnergy( pair.second );

        // associate this hitset and hit to the geant4 hit key
        hittruthassoc->addAssoc(hitsetkey, hitkey, g4hit_it->first);
      }

    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
void PHG4MicromegasHitReco::SetDefaultParameters()
{
  // default timing window (ns)
  /*
   * see https://indico.bnl.gov/event/8548/contributions/37753/attachments/28212/43343/2020_05_Proposal_sPhenixMonitoring_update_19052020.pptx slide 10
   * small negative time for tmin is set to catch out of time, same-bunch pileup events
   * similar value is used in PHG4InttReco
  */
  set_default_double_param("micromegas_tmin", -20 );
  set_default_double_param("micromegas_tmax", 800 );

  // gas data from
  // http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
  // assume Ar/iC4H10 90/10, at 20C and 1atm
  // dedx (KeV/cm) for MIP
  static constexpr double Ar_dEdx = 2.44;
  static constexpr double iC4H10_dEdx = 5.93;
  static constexpr double mix_dEdx = 0.9*Ar_dEdx + 0.1*iC4H10_dEdx;

  // number of electrons per MIP (cm-1)
  static constexpr double Ar_ntot = 94;
  static constexpr double iC4H10_ntot = 195;
  static constexpr double mix_ntot = 0.9*Ar_ntot + 0.1*iC4H10_ntot;

  // number of electrons per gev
  static constexpr double mix_electrons_per_gev = 1e6*mix_ntot / mix_dEdx;
  set_default_double_param("micromegas_electrons_per_gev", mix_electrons_per_gev );

  // gain
  set_default_double_param("micromegas_gain", 2000 );

  // electron cloud sigma, after avalanche (cm)
  set_default_double_param("micromegas_cloud_sigma", 0.04 );

  // transverse diffusion (cm/sqrt(cm))
  set_default_double_param("micromegas_diffusion_trans", 0.03 );
  
  // zigzag strips
  set_default_int_param("micromegas_zigzag_strips", true );
}

//___________________________________________________________________________
void PHG4MicromegasHitReco::setup_tiles(PHCompositeNode* topNode)
{

  // get geometry
  const auto geonodename_full = full_geonodename();
  auto geonode_full = findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename_full);
  if (!geonode_full)
  {
    // if full geometry (cylinder + tiles) do not exist, try create it from bare geometry (cylinder only)
    const auto geonodename_bare = bare_geonodename();
    auto geonode_bare =  findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename_bare);
    if( !geonode_bare )
    {
      std::cout << PHWHERE << "Could not locate geometry node " << geonodename_bare << std::endl;
      exit(1);
    }

    // create new node
    if( Verbosity() )
    { std::cout << "PHG4MicromegasHitReco::setup_tiles - " << PHWHERE << " creating node " << geonodename_full << std::endl; }

    geonode_full = new PHG4CylinderGeomContainer();
    PHNodeIterator iter(topNode);
    auto runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
    auto newNode = new PHIODataNode<PHObject>(geonode_full, geonodename_full, "PHObject");
    runNode->addNode(newNode);

    // copy cylinders
    PHG4CylinderGeomContainer::ConstRange range = geonode_bare->get_begin_end();
    for( auto iter = range.first; iter != range.second; ++iter )
    {
      const auto layer = iter->first;
      const auto cylinder = static_cast<CylinderGeomMicromegas*>(iter->second);
      geonode_full->AddLayerGeom( layer, new CylinderGeomMicromegas( *cylinder ) );
    }
  }

  // get cylinders
  PHG4CylinderGeomContainer::ConstRange range = geonode_full->get_begin_end();
  for( auto iter = range.first; iter != range.second; ++iter )
  {
    std::cout << "PHG4MicromegasHitReco::setup_tiles - processing layer " << iter->first << std::endl;
    auto cylinder = static_cast<CylinderGeomMicromegas*>(iter->second);

    // assign tiles
    cylinder->set_tiles( m_tiles );

    /*
     * asign segmentation type and pitch
     * assume first layer in phi, other(s) are z
     */
    const bool is_first( iter == range.first );
    cylinder->set_segmentation_type( is_first ?
      MicromegasDefs::SegmentationType::SEGMENTATION_PHI :
      MicromegasDefs::SegmentationType::SEGMENTATION_Z );

    /*
     * assign drift direction
     * assume first layer is outward, with readout plane at the top, and others are inward, with readout plane at the bottom
     * this is used to properly implement transverse diffusion in ::process_event
     */
    cylinder->set_drift_direction( is_first ? 
      MicromegasDefs::DriftDirection::OUTWARD :
      MicromegasDefs::DriftDirection::INWARD );      
    
    // pitch
    /* they correspond to 256 channels along the phi direction, and 256 along the z direction, assuming 25x50 tiles */
    cylinder->set_pitch( is_first ? 25./256 : 50./256 );
    cylinder->identify( std::cout );
  }
}

//___________________________________________________________________________
uint PHG4MicromegasHitReco::get_primary_electrons( PHG4Hit* g4hit ) const
{ return gsl_ran_poisson(m_rng.get(), g4hit->get_eion() * m_electrons_per_gev); }

//___________________________________________________________________________
uint PHG4MicromegasHitReco::get_single_electron_amplification() const
{
  /*
   * to handle gain fluctuations, an exponential distribution is used, similar to what used for the GEMS
   * (simulations/g4simulations/g4tpc/PHG4TpcPadPlaneReadout::getSingleEGEMAmplification)
   * One must get a different random number for each primary electron for this to be valid
   */
  return gsl_ran_exponential(m_rng.get(), m_gain);
}

//___________________________________________________________________________
PHG4MicromegasHitReco::charge_list_t PHG4MicromegasHitReco::distribute_charge(
  CylinderGeomMicromegas* layergeom,
  uint tileid,
  const TVector3& location,
  double sigma ) const
{

  // find tile and strip matching center position
  auto stripnum = layergeom->find_strip( tileid, location );

  // check tile and strip
  if( stripnum < 0 ) return charge_list_t();

  // store pitch and radius
  const auto pitch = layergeom->get_pitch();
  const auto radius = layergeom->get_radius();

  // find relevant strip indices
  const auto strip_count = layergeom->get_strip_count( tileid );
  const auto stripnum_min = clamp<int>( stripnum - 5.*sigma/pitch - 1, 0, strip_count );
  const auto stripnum_max = clamp<int>( stripnum + 5*sigma/pitch + 1, 0, strip_count );

  // prepare charge list
  charge_list_t charge_list;

  // store azimuthal angle
  const auto phi = std::atan2( location.y(), location.x() );

  // loop over strips
  for( int strip = stripnum_min; strip <= stripnum_max; ++strip )
  {
    // get strip center
    const TVector3 strip_location = layergeom->get_world_coordinate( tileid, strip );
    const auto phi_strip = std::atan2( strip_location.y(), strip_location.x() );

    // find relevant strip coordinate with respect to location
    const auto xloc = layergeom->get_segmentation_type() == MicromegasDefs::SegmentationType::SEGMENTATION_PHI ?
      radius * bind_angle( phi_strip - phi ):
      strip_location.z() - location.z();

    // calculate charge fraction
    const auto fraction = m_zigzag_strips ?
      get_zigzag_fraction( xloc, sigma, pitch ):
      get_rectangular_fraction( xloc, sigma, pitch );

    // store
    charge_list.push_back( std::make_pair( strip, fraction ) );
  }

  return charge_list;
}
