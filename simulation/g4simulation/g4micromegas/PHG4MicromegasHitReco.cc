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
#include <g4main/PHG4Utils.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TVector3.h>

#include <gsl/gsl_randist.h>
#include <cassert>
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

  PHNodeIterator iter(topNode);

  // store parameters
  m_tmin = get_double_param("micromegas_tmin" );
  m_tmax = get_double_param("micromegas_tmax" );
  m_electrons_per_gev = get_double_param("micromegas_electrons_per_gev" );
  m_gain = get_double_param("micromegas_gain");
  m_cloud_sigma = get_double_param("micromegas_cloud_sigma" );

  std::cout
    << "PHG4MicromegasHitReco::InitRun\n"
    << " m_tmin: " << m_tmin << "ns, m_tmax: " << m_tmax << "ns\n"
    << " m_electrons_per_gev: " << m_electrons_per_gev << "\n"
    << " m_gain: " << m_gain << "\n"
    << " m_cloud_sigma: " << m_cloud_sigma << "cm\n"
    << std::endl;

  // setup tiles
  setup_tiles( topNode );

  // get dst node
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
  const std::string geonodename = "CYLINDERGEOM_" + m_detector;
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

    // get corresponding hits
    PHG4HitContainer::ConstRange g4hit_range = g4hitcontainer->getHits(layer);

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
      auto world_mid = (world_in + world_out)*0.5;

      // distribute charge among adjacent strips
      // todo: will need to split charge for different segments along drift space
      // calculate the fired strips for each segment, including diffusion
      // get the complete cluster by resumming
      int tileid;
      charge_list_t fractions;
      std::tie(tileid, fractions) = distribute_charge( layergeom, world_mid, m_cloud_sigma );
      if( tileid < 0 || fractions.empty() ) continue;

      // make sure fractions adds up to unity
      if( Verbosity() > 0 )
      {
        const auto sum = std::accumulate( fractions.begin(), fractions.end(), double( 0 ),
          []( double value, const charge_pair_t& pair ) { return value + pair.second; } );
        std::cout << "PHG4MicromegasHitReco::process_event - sum: " << sum << std::endl;
      }

      /*
      calculate the total number of electrons used for this hit
      this is what will be stored as 'energy' in the hits
      this accounts for number of 'primary', detector gain, and fluctuations thereof
      */
      const auto nelectrons = get_electrons( g4hit );      
      
      // create hitset
      TrkrDefs::hitsetkey hitsetkey = MicromegasDefs::genHitSetKey( layer, tileid );
      auto hitset_it = trkrhitsetcontainer->findOrAddHitSet(hitsetkey);

      // generate the key for this hit
      // loop over strips in list
      for( const auto pair:fractions )
      {
        // get strip and bound check
        const int strip = pair.first;
        if( strip < 0 || strip >= layergeom->get_strip_count( tileid ) ) continue;

        // get hit from hitset
        TrkrDefs::hitkey hitkey = MicromegasDefs::genHitKey(strip);
        TrkrHit* hit = hitset_it->second->getHit(hitkey);
        if( !hit )
        {
          // create hit and insert in hitset
          hit = new TrkrHit;
          hitset_it->second->addHitSpecificKey(hitkey, hit);
        }

        // add energy from g4hit
        const double fraction = pair.second;
        hit->addEnergy( fraction*nelectrons );

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
  set_default_double_param("micromegas_tmin", -5000 );
  set_default_double_param("micromegas_tmax", 5000 );

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

}

//___________________________________________________________________________
void PHG4MicromegasHitReco::setup_tiles(PHCompositeNode* topNode)
{

  // get geometry
  const std::string geonodename = "CYLINDERGEOM_" + m_detector;
  auto geonode = findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename);
  if (!geonode)
  {
    std::cout << PHWHERE << "Could not locate geometry node " << geonodename << std::endl;
    exit(1);
  }

  // get cylinders
  PHG4CylinderGeomContainer::ConstRange range = geonode->get_begin_end();
  for( auto iter = range.first; iter != range.second; ++iter )
  {
    std::cout << "PHG4MicromegasHitReco::setup_tiles - processing layer " << iter->first << std::endl;
    auto cylinder = static_cast<CylinderGeomMicromegas*>(iter->second);

    // assign tiles
    cylinder->set_tiles( m_tiles );

    // asign segmentation type and pitch
    // assume first layer in phi, other(s) are z
    const bool is_first( iter == range.first );
    cylinder->set_segmentation_type( is_first ?
      MicromegasDefs::SegmentationType::SEGMENTATION_PHI :
      MicromegasDefs::SegmentationType::SEGMENTATION_Z );

    // pitch
    /* they correspond to 256 channels along the phi direction, and 256 along the z direction, assuming 25x50 tiles */
    cylinder->set_pitch( is_first ? 25./256 : 50./256 );
    cylinder->identify( std::cout );
  }
}

//___________________________________________________________________________
uint PHG4MicromegasHitReco::get_electrons( PHG4Hit* g4hit ) const
{
  // number of primary electrons
  uint nprimary = gsl_ran_poisson(m_rng.get(), g4hit->get_eion() * m_electrons_per_gev);
  if( !nprimary ) return 0;

  /*
   * to handle gain fluctuations, an exponential distribution is used, similar to what used for the GEMS
   * (simulations/g4simulations/g4tpc/PHG4TpcPadPlaneReadout::getSingleEGEMAmplification)
   * One must get a different random number for each primary electron for this to be valid
   */
  uint ntot = 0;
  for( uint i = 0; i < nprimary; ++i ) { ntot += gsl_ran_exponential(m_rng.get(), m_gain); }
  
  if( Verbosity() > 0 )
  {
    std::cout 
      << "PHG4MicromegasHitReco::get_electrons -"
      << " nprimary: " << nprimary 
      << " average gain: " << static_cast<double>(ntot)/nprimary 
      << std::endl;
  }
  
  return ntot;
}

//___________________________________________________________________________
PHG4MicromegasHitReco::charge_info_t PHG4MicromegasHitReco::distribute_charge(
  CylinderGeomMicromegas* layergeom, const TVector3& location, double sigma ) const
{

  // find tile and strip matching center position
  int tileid, stripnum;
  std::tie(tileid, stripnum) = layergeom->find_strip( location );

  // check tile and strip
  if( tileid < 0 || stripnum < 0 ) return std::make_pair( -1, charge_list_t() );

  // store pitch and radius
  const auto pitch = layergeom->get_pitch();
  const auto radius = layergeom->get_radius();

  // find relevant strip indices
  const auto strip_count = layergeom->get_strip_count( tileid );
  const auto stripnum_min = std::clamp<int>( stripnum - 5.*sigma/pitch - 1, 0, strip_count );
  const auto stripnum_max = std::clamp<int>( stripnum + 5*sigma/pitch + 1, 0, strip_count );

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
    const auto fraction = get_zigzag_fraction( xloc, sigma, pitch );

    // store
    charge_list.push_back( std::make_pair( strip, fraction ) );
  }

  return std::make_pair( tileid, charge_list );
}
