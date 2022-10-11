/*!
 * \file PHG4MicromegasHitReco.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHG4MicromegasHitReco.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                     // for SubsysReco

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4HitContainer.h>

#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>


#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHObject.h>                         // for PHObject

#include <phparameter/PHParameterInterface.h>       // for PHParameterInterface

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitv2.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitTruthAssocv1.h>

#include <TVector2.h>
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

  // TVector3 streamer
  [[maybe_unused]] inline std::ostream& operator << (std::ostream& out, const TVector3& position)
  {
    out << "(" << position.x() << ", " << position.y() << ", " << position.z() << ")";
    return out;
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
  m_added_smear_sigma_z = get_double_param("micromegas_added_smear_sigma_z");
  m_added_smear_sigma_rphi = get_double_param("micromegas_added_smear_sigma_rphi");

  // printout
  std::cout
    << "PHG4MicromegasHitReco::InitRun\n"
    << " m_tmin: " << m_tmin << "ns, m_tmax: " << m_tmax << "ns\n"
    << " m_electrons_per_gev: " << m_electrons_per_gev << "\n"
    << " m_gain: " << m_gain << "\n"
    << " m_cloud_sigma: " << m_cloud_sigma << "cm\n"
    << " m_diffusion_trans: " << m_diffusion_trans << "cm/sqrt(cm)\n"
    << " m_added_smear_sigma_z: " << m_added_smear_sigma_z << "cm\n"
    << " m_added_smear_sigma_rphi: " << m_added_smear_sigma_rphi << "cm\n"
    << std::endl;

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
    hitsetcontainer = new TrkrHitSetContainerv1;
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

    hittruthassoc = new TrkrHitTruthAssocv1;
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

  // acts geometry
  m_acts_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  assert( m_acts_geometry );

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
     * get the position of the detector mesh in local coordinates
     * in local coordinate the mesh is a plane perpendicular to the z axis
     * Its position along z depends on the drift direction
     * it is used to calculate the drift distance of the primary electrons, and the
     * corresponding transverse diffusion
     */
    const auto mesh_local_z = layergeom->get_drift_direction() == MicromegasDefs::DriftDirection::OUTWARD ?
      layergeom->get_thickness()/2:
      -layergeom->get_thickness()/2;

//     /*
//      * get the radius of the detector mesh. It depends on the drift direction
//      * it is used to calculate the drift distance of the primary electrons, and the
//      * corresponding transverse diffusion
//      */
//       const auto mesh_radius = layergeom->get_drift_direction() == MicromegasDefs::DriftDirection::OUTWARD ?
//       (layergeom->get_radius() + layergeom->get_thickness()/2):
//       (layergeom->get_radius() - layergeom->get_thickness()/2);

    // get hits
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

      // get tile id from g4hit
      const int tileid = g4hit->get_property_int( PHG4Hit::prop_index_i );
      
      // convert to local coordinate
      const auto local_in = layergeom->get_local_from_world_coords( tileid, m_acts_geometry, world_in );
      const auto local_out = layergeom->get_local_from_world_coords( tileid, m_acts_geometry, world_out );

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
        // in local reference frame, drift occurs along the y axis, from local y to mesh_local_z

        const auto t = gsl_ran_flat(m_rng.get(), 0.0, 1.0);
        auto local = local_in*t + local_out*(1.0-t);
        
        if( m_diffusion_trans > 0 )
        {
          // add transeverse diffusion
          // first convert to polar coordinates
          const double z = local.z();
          const double drift_distance = std::abs(z - mesh_local_z);
          const double diffusion = gsl_ran_gaussian(m_rng.get(), m_diffusion_trans*std::sqrt(drift_distance));
          const double diffusion_angle = gsl_ran_flat(m_rng.get(), -M_PI, M_PI);

          // diffusion occurs in x,z plane with a magnitude 'diffusion' and an angle 'diffusion angle'
          local += TVector3( diffusion*std::cos(diffusion_angle), diffusion*std::sin(diffusion_angle), 0 );
        }

        const auto& added_smear_sigma =  layergeom->get_segmentation_type() == MicromegasDefs::SegmentationType::SEGMENTATION_PHI ?
          m_added_smear_sigma_rphi: m_added_smear_sigma_z;

        if( added_smear_sigma > 0 )
        {
          // additional ad hoc smearing
          const double added_smear_trans = gsl_ran_gaussian(m_rng.get(), added_smear_sigma);
          const double added_smear_angle = gsl_ran_flat(m_rng.get(), -M_PI, M_PI);
          local += TVector3( added_smear_trans*std::cos(added_smear_angle), added_smear_trans*std::sin(added_smear_angle), 0 );
        }
        
        // distribute charge among adjacent strips
        const auto fractions = distribute_charge( layergeom, tileid, { local.x(), local.y() }, m_cloud_sigma );

        // make sure fractions adds up to unity
        if( Verbosity() > 10 )
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
          if( strip < 0 || strip >= (int) layergeom->get_strip_count( tileid, m_acts_geometry ) ) continue;

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

  // additional smearing (cm)
  set_default_double_param("micromegas_added_smear_sigma_z", 0);
  set_default_double_param("micromegas_added_smear_sigma_rphi", 0);
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
  const TVector2& local_coords,
  double sigma ) const
{

  // find tile and strip matching center position
  auto stripnum = layergeom->find_strip_from_local_coords( tileid, m_acts_geometry, local_coords );

  // check tile and strip
  if( stripnum < 0 ) return charge_list_t();

  // store pitch and radius
  const auto pitch = layergeom->get_pitch();

  // find relevant strip indices
  const auto strip_count = layergeom->get_strip_count( tileid, m_acts_geometry );
  const int nstrips = 5.*sigma/pitch + 1;
  const auto stripnum_min = std::clamp<int>( stripnum - nstrips, 0, strip_count );
  const auto stripnum_max = std::clamp<int>( stripnum + nstrips, 0, strip_count );

  // prepare charge list
  charge_list_t charge_list;

  // loop over strips
  for( int strip = stripnum_min; strip <= stripnum_max; ++strip )
  {
    // get strip center
    const auto strip_location = layergeom->get_local_coordinates( tileid, m_acts_geometry, strip );

    /*
     * find relevant strip coordinate with respect to location
     * in local coordinate, phi segmented view has strips along z and measures along x
     * in local coordinate, z segmented view has strips along phi and measures along y
     */
    const auto xloc = layergeom->get_segmentation_type() == MicromegasDefs::SegmentationType::SEGMENTATION_PHI ?
      (strip_location.X() - local_coords.X()):
      (strip_location.Y() - local_coords.Y());

    // decide of whether zigzag or straight strips are used depending on segmentation type
    /*
     * for the real detector SEGMENTATION_Z view has zigzag strip due to large pitch (2mm)
     * whereas SEGMENTATION_PHI has straight strips
     */
    const bool zigzag_strips = (layergeom->get_segmentation_type() == MicromegasDefs::SegmentationType::SEGMENTATION_Z );
    
    // calculate charge fraction
    const auto fraction = zigzag_strips ?
      get_zigzag_fraction( xloc, sigma, pitch ):
      get_rectangular_fraction( xloc, sigma, pitch );

    // store
    charge_list.push_back( std::make_pair( strip, fraction ) );
  }

  return charge_list;
}
