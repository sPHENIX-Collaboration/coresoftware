/*!
 *  \file PHSimpleKFProp.cc
 *  \brief		kalman filter based propagator
 *  \author Michael Peters & Christof Roland
 */

#include "PHSimpleKFProp.h"
#include "ALICEKF.h"
#include "GPUTPCTrackLinearisation.h"
#include "GPUTPCTrackParam.h"
#include "PHGhostRejection.h"
#include "nanoflann.hpp"

#include <ffamodules/CDBInterface.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <phfield/PHField.h>
#include <phfield/PHFieldConfig.h>
#include <phfield/PHFieldConfigv1.h>
#include <phfield/PHFieldUtility.h>

// tpc distortion correction
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

#include <TSystem.h>

#include <Geant4/G4SystemOfUnits.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <omp.h>

#include <filesystem>
#include <iostream>
#include <syncstream>
#include <vector>

// anonymous namespace for local functions
namespace
{
  // square
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }
}  // namespace

using keylist = std::vector<TrkrDefs::cluskey>;

PHSimpleKFProp::PHSimpleKFProp(const std::string& name)
  : SubsysReco(name)
{}

//______________________________________________________
PHSimpleKFProp::~PHSimpleKFProp()
{
  if( m_own_fieldmap )
  { delete _field_map; }
}

//______________________________________________________
int PHSimpleKFProp::End(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSimpleKFProp::InitRun(PHCompositeNode* topNode)
{
  int ret = get_nodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  PHFieldConfigv1 fcfg;
  fcfg.set_field_config(PHFieldConfig::FieldConfigTypes::Field3DCartesian);

  if (std::filesystem::path(m_magField).extension() != ".root")
  { m_magField = CDBInterface::instance()->getUrl(m_magField); }

  if (!_use_const_field)
  {
    if (!std::filesystem::exists(m_magField))
    {
      if (m_magField.empty())
      { m_magField = "empty string"; }
      std::cout << PHWHERE << "Fieldmap " << m_magField << " does not exist" << std::endl;
      gSystem->Exit(1);
    }

    fcfg.set_filename(m_magField);
  }
  else
  {
    fcfg.set_field_config(PHFieldConfig::FieldConfigTypes::kFieldUniform);
    fcfg.set_magfield_rescale(_const_field);
  }

  // compare field config from that on node tree
  /*
   * if the magnetic field is already on the node tree PHFieldUtility::GetFieldConfigNode returns the existing configuration.
   * One must then check wheter the two configurations are identical, to decide whether one must use the field from node tree or create our own.
   * Otherwise the configuration passed as argument is stored on the node tree.
   */
  const auto node_fcfg = PHFieldUtility::GetFieldConfigNode(&fcfg, topNode);
  if( fcfg == *node_fcfg )
  {
    // both configurations are identical, use field map from node tree
    std::cout << "PHSimpleKFProp::InitRun - using node tree field map" << std::endl;
    _field_map = PHFieldUtility::GetFieldMapNode(&fcfg, topNode);
    m_own_fieldmap = false;
  } else {
    // both configurations differ. Use our own field map
    std::cout << "PHSimpleKFProp::InitRun - using own field map" << std::endl;
    _field_map = PHFieldUtility::BuildFieldMap(&fcfg);
    m_own_fieldmap = true;
  }

  // alice kalman filter
  fitter = std::make_unique<ALICEKF>(topNode, _cluster_map, _field_map, _fieldDir, _min_clusters_per_track, _max_sin_phi, Verbosity());
  fitter->setNeonFraction(Ne_frac);
  fitter->setArgonFraction(Ar_frac);
  fitter->setCF4Fraction(CF4_frac);
  fitter->setNitrogenFraction(N2_frac);
  fitter->setIsobutaneFraction(isobutane_frac);
  fitter->useConstBField(_use_const_field);
  fitter->setConstBField(_const_field);
  fitter->useFixedClusterError(_use_fixed_clus_err);
  fitter->setFixedClusterError(0, _fixed_clus_err.at(0));
  fitter->setFixedClusterError(1, _fixed_clus_err.at(1));
  fitter->setFixedClusterError(2, _fixed_clus_err.at(2));
  // _field_map = PHFieldUtility::GetFieldMapNode(nullptr,topNode);
  // m_Cache = magField->makeCache(m_tGeometry->magFieldContext);

  // assign number of threads
  std::cout << "PHSimpleKFProp::InitRun - m_num_threads: " << m_num_threads << std::endl;
  if( m_num_threads >= 1 ) { omp_set_num_threads( m_num_threads ); }

  return Fun4AllReturnCodes::EVENT_OK;
}

double PHSimpleKFProp::get_Bz(double x, double y, double z) const
{
  if (_use_const_field || std::abs(z) > 105.5)
  {
    return _const_field;
  }
  double p[4] = {x * cm, y * cm, z * cm, 0. * cm};
  double bfield[3];

  // check thread number. Use uncached field accessor for all but thread 0.
  if( omp_get_thread_num() == 0 )
  {
    _field_map->GetFieldValue(p, bfield);
  } else {
    _field_map->GetFieldValue_nocache(p, bfield);
  }

  /*  Acts::Vector3 loc(0,0,0);
  int mfex = (magField != nullptr);
  int tgex = (m_tGeometry != nullptr);
  std::cout << " getting acts field " << mfex << " " << tgex << std::endl;
  auto Cache =  m_tGeometry->magField->makeCache(m_tGeometry->magFieldContext);
  auto bf = m_tGeometry->magField->getField(loc,Cache);
  if(bf.ok())
    {
      Acts::Vector3 val = bf.value();
      std::cout << "bz big: " <<  bfield[2]/tesla << " bz acts: " << val(2) << std::endl;
    }
  */
  return bfield[2] / tesla;
}

int PHSimpleKFProp::get_nodes(PHCompositeNode* topNode)
{

  // acts geometry
  m_tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tgeometry)
  {
    std::cout << "No Acts tracking geometry, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc global position wrapper
  m_globalPositionWrapper.loadNodes(topNode);

  // clusters
  if (_use_truth_clusters)
  {
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_TRUTH");
  }
  else
  {
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  }
  if (!_cluster_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tracks
  _track_map = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!_track_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find TrackSeedContainer " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc grometry
  auto geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cerr << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  for (int i = 7; i <= 54; i++)
  {
    radii.push_back(geom_container->GetLayerCellGeom(i)->get_radius());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSimpleKFProp::process_event(PHCompositeNode* topNode)
{
  if (_n_iteration != 0)
  {
    _iteration_map = findNode::getClass<TrkrClusterIterationMapv1>(topNode, "CLUSTER_ITERATION_MAP");
    if (!_iteration_map)
    {
      std::cerr << PHWHERE << "Cluster Iteration Map missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  PHTimer timer("KFPropTimer");
  timer.restart();

  // check number of seeds against maximum
  if(_max_seeds > 0 && _track_map->size() > _max_seeds)
  {
    std::cout << "number of TPC seeds: " << _track_map->size() << std::endl;
    std::cout << PHWHERE << "number of TPC seeds > " << _max_seeds << " aborting event." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  const auto globalPositions = PrepareKDTrees();
  if (Verbosity())
  { std::cout << "PHSimpleKFProp::process_event - PrepareKDTrees time: " << timer.elapsed() << " ms" << std::endl; }

  // list of cluster chains
  std::vector<std::vector<TrkrDefs::cluskey>> new_chains;
  std::vector<TrackSeed_v2> unused_tracks;

  timer.restart();
  #pragma omp parallel
  {
    if (Verbosity())
    {
      std::osyncstream(std::cout)
        << "PHSimpleKFProp -"
        << " num_threads: " << omp_get_num_threads()
        << " this thread: " << omp_get_thread_num()
        << std::endl;
    }

    PHTimer timer_mp("KFPropTimer_parallel");

    std::vector<std::vector<TrkrDefs::cluskey>> local_chains;
    std::vector<TrackSeed_v2> local_unused;

    #pragma omp for schedule(static)
    for (size_t track_it = 0; track_it != _track_map->size(); ++track_it)
    {
      if (Verbosity())
      {
        std::osyncstream(std::cout)
          << "PHSimpleKFProp -"
          << " num_threads: " << omp_get_num_threads()
          << " this thread: " << omp_get_thread_num()
          << " processing seed " << track_it << std::endl;
      }

      // if not a TPC track, ignore
      auto track = _track_map->get(track_it);
      const bool is_tpc = std::any_of(
        track->begin_cluster_keys(),
        track->end_cluster_keys(),
        [](const TrkrDefs::cluskey& key)
      { return TrkrDefs::getTrkrId(key) == TrkrDefs::tpcId; });

      if (is_tpc)
      {

        // copy list of seed cluster keys
        std::vector<std::vector<TrkrDefs::cluskey>> keylist_A(1);
        std::copy(track->begin_cluster_keys(), track->end_cluster_keys(), std::back_inserter(keylist_A[0]));

        // copy seed clusters position into local map
        std::map<TrkrDefs::cluskey, Acts::Vector3> trackClusPositions;
        std::transform(track->begin_cluster_keys(), track->end_cluster_keys(), std::inserter(trackClusPositions, trackClusPositions.end()),
          [globalPositions](const auto& key)
        { return std::make_pair(key, globalPositions.at(key)); });

        /// Can't circle fit a seed with less than 3 clusters, skip it
        if (keylist_A[0].size() < 3)
        {
          continue;
        }

        /// This will by definition return a single pair with each vector
        /// in the pair length 1 corresponding to the seed info
        std::vector<float> trackChi2;

        timer_mp.restart();
        auto seedpair = fitter->ALICEKalmanFilter(keylist_A, false, trackClusPositions, trackChi2);

        if (Verbosity() > 3)
        {
          std::cout << "PHSimpleKFProp::process_event - single track ALICEKF time " << timer_mp.elapsed() << " ms" << std::endl;
        }

        timer_mp.restart();

        /// circle fit back to update track parameters
        TrackSeedHelper::circleFitByTaubin(track, trackClusPositions, 7, 55);
        TrackSeedHelper::lineFit(track, trackClusPositions, 7, 55);
        track->set_phi(TrackSeedHelper::get_phi(track, trackClusPositions));
        if (Verbosity() > 3)
        {
          std::cout << "PHSimpleKFProp::process_event - single track circle fit time " << timer_mp.elapsed() << " ms" << std::endl;
        }

        if (seedpair.first.empty()|| seedpair.second.empty())
        {
          continue;
        }

        if (Verbosity())
        {
          std::cout << "is tpc track" << std::endl;
        }

        timer_mp.restart();

        if (Verbosity())
        {
          std::cout << "propagate first round" << std::endl;
        }

        auto preseed = PropagateTrack(track, PropagationDirection::Inward, seedpair.second.at(0), globalPositions);
        if (Verbosity())
        {
          std::cout << "preseed size " << preseed.size() << std::endl;
        }

        std::vector<std::vector<TrkrDefs::cluskey>> kl = {preseed};
        if (Verbosity())
        {
          std::cout << "kl size " << kl.size() << std::endl;
        }
        std::vector<float> pretrackChi2;

        auto prepair = fitter->ALICEKalmanFilter(kl, false, globalPositions, pretrackChi2);
        if (prepair.first.empty() || prepair.second.empty())
        {
          continue;
        }

        std::reverse(kl.at(0).begin(), kl.at(0).end());

        auto pretrack = prepair.first.at(0);

        // copy seed clusters position into local map
        std::map<TrkrDefs::cluskey, Acts::Vector3> pretrackClusPositions;
        std::transform(pretrack.begin_cluster_keys(), pretrack.end_cluster_keys(), std::inserter(pretrackClusPositions, pretrackClusPositions.end()),
          [globalPositions](const auto& key)
          { return std::make_pair(key, globalPositions.at(key)); });

        // fit seed
        TrackSeedHelper::circleFitByTaubin(&pretrack,pretrackClusPositions, 7, 55);
        TrackSeedHelper::lineFit(&pretrack, pretrackClusPositions, 7, 55);
        pretrack.set_phi(TrackSeedHelper::get_phi(&pretrack, pretrackClusPositions));

        prepair.second.at(0).SetDzDs(-prepair.second.at(0).GetDzDs());
        const auto finalchain = PropagateTrack(&pretrack, kl.at(0), PropagationDirection::Outward, prepair.second.at(0), globalPositions);

        if (finalchain.size() > kl.at(0).size())
        {
          local_chains.push_back(std::move(finalchain));
        }
        else
        {
          local_chains.push_back(std::move(kl.at(0)));
        }

        if (Verbosity() > 3)
        {
          std::cout << "PHSimpleKFProp::process_event - propagate track time " << timer_mp.elapsed() << " ms" << std::endl;
        }
      }
      else
      {
        if (Verbosity())
        {
          std::cout << "is NOT tpc track" << std::endl;
        }
        local_unused.emplace_back(*track);
      }
    }

    // sort list and remove duplicates
    std::sort(local_chains.begin(),local_chains.end());
    local_chains.erase(std::unique(local_chains.begin(),local_chains.end()),local_chains.end());

    // Critical sections to merge thread-local results
    #pragma omp critical
    {
      new_chains.reserve(new_chains.size()+local_chains.size());
      new_chains.insert(new_chains.end(), std::make_move_iterator(local_chains.begin()), std::make_move_iterator(local_chains.end()));

      unused_tracks.reserve(unused_tracks.size()+local_unused.size());
      unused_tracks.insert(unused_tracks.end(), std::make_move_iterator(local_unused.begin()), std::make_move_iterator(local_unused.end()));
    }
  }
  if (Verbosity())
  { std::cout << "PHSimpleKFProp::process_event - first seed loop time: " << timer.elapsed() << " ms" << std::endl; }

  // sort merged list and remove duplicates
  timer.restart();
  std::sort(new_chains.begin(),new_chains.end());
  new_chains.erase(std::unique(new_chains.begin(),new_chains.end()),new_chains.end());

  if (Verbosity())
  { std::cout << "PHSimpleKFProp::process_event - first cleanup time: " << timer.elapsed() << " ms" << std::endl; }

  if( Verbosity() )
  {
    // print all seeds
    std::cout << "PHSimpleKFProp::process_event - new_chains size: " << new_chains.size() << std::endl;
    for( const auto& chain:new_chains )
    {
      std::cout << "PHSimpleKFProp::process_event - { ";
      for( const auto& key:chain )
      {
        std::cout << key << " ";
      }
      std::cout << "}" << std::endl << std::endl;
    }
  }

  // erase all seeds for size 2 or less
  new_chains.erase( std::remove_if( new_chains.begin(), new_chains.end(),
    [](const auto& chain) { return chain.size()<3; } ),
    new_chains.end() );

  // re-run ALICE Kalman Filter on completed chains
  timer.restart();
  std::vector<float> trackChi2;
  auto seeds = fitter->ALICEKalmanFilter(new_chains, true, globalPositions, trackChi2);
  if (Verbosity())
  {  std::cout << "PHSimpleKFProp::process_event - ALICEKalmanFilter time: " << timer.elapsed() << " ms" << std::endl; }

  // reset track map
  _track_map->Reset();

  //  Move ghost rejection into publishSeeds, so that we don't publish rejected seeds
  timer.restart();
  if (m_ghostrejection)
  {
    rejectAndPublishSeeds(seeds.first, globalPositions, trackChi2);
  }
  else
  {
    publishSeeds(seeds.first);
  }

  // also publish unused seeds (not TPC)
  publishSeeds(unused_tracks);

  if (Verbosity())
  { std::cout << "PHSimpleKFProp::process_event - publishSeeds time: " << timer.elapsed() << " ms" << std::endl; }

  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 PHSimpleKFProp::getGlobalPosition(TrkrDefs::cluskey key, TrkrCluster* cluster) const
{
  // get global position from Acts transform
  return _pp_mode ?
    m_tgeometry->getGlobalPosition(key, cluster):
    m_globalPositionWrapper.getGlobalPositionDistortionCorrected( key, cluster, 0 );
}

PositionMap PHSimpleKFProp::PrepareKDTrees()
{
  PositionMap globalPositions;
  //***** convert clusters to kdhits, and divide by layer
  std::vector<std::vector<std::vector<double>>> kdhits;
  kdhits.resize(58);
  if (!_cluster_map)
  {
    std::cout << "WARNING: (tracking.PHTpcTrackerUtil.convert_clusters_to_hits) cluster map is not provided" << std::endl;
    return globalPositions;
  }

  for (const auto& hitsetkey : _cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
  {
    auto range = _cluster_map->getClusters(hitsetkey);
    for (TrkrClusterContainer::ConstIterator it = range.first; it != range.second; ++it)
    {
      const auto& [cluskey,cluster] = *it;
      if (!cluster)
      {
        continue;
      }

      // skip hits used in a previous iteration
      if (_n_iteration && _iteration_map && _iteration_map->getIteration(cluskey) > 0)
      {
        continue;
      }

      const auto globalpos = getGlobalPosition(cluskey, cluster);
      globalPositions.emplace(cluskey, globalpos);

      const int layer = TrkrDefs::getLayer(cluskey);
      std::vector<double> kdhit{ globalpos.x(), globalpos.y(),  globalpos.z(), 0 };
      const uint64_t key = cluskey;
      std::memcpy(&kdhit[3], &key, sizeof(key));

      //      HINT: way to get original uint64_t value from double:
      //
      //      LOG_DEBUG("tracking.PHTpcTrackerUtil.convert_clusters_to_hits")
      //        << "orig: " << cluster->getClusKey() << ", readback: " << (*((int64_t*)&kdhit[3]));

      kdhits[layer].push_back(std::move(kdhit));
    }
  }
  _ptclouds.resize(kdhits.size());
  _kdtrees.resize(kdhits.size());
  for (size_t l = 0; l < kdhits.size(); ++l)
  {
    if (Verbosity() > 1)
    {
      std::cout << "l: " << l << std::endl;
    }
    _ptclouds[l] = std::make_shared<KDPointCloud<double>>();
    _ptclouds[l]->pts = std::move(kdhits[l]);

    _kdtrees[l] = std::make_shared<nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, KDPointCloud<double>>, KDPointCloud<double>, 3>>(3, *(_ptclouds[l]), nanoflann::KDTreeSingleIndexAdaptorParams(10));
    _kdtrees[l]->buildIndex();
  }

  return globalPositions;
}

bool PHSimpleKFProp::TransportAndRotate(
  double old_radius,
  double new_radius,
  double& phi,
  GPUTPCTrackParam& kftrack,
  GPUTPCTrackParam::GPUTPCTrackFitParam& fp) const
{
  if (Verbosity() > 2)
  {
    std::cout << "old_radius " << old_radius << ", new_radius " << new_radius << std::endl;
  }
  const float transport_spacing = .05;  // max radial distance between transport points
  const int Ndivisions = floor(std::abs(new_radius - old_radius) / transport_spacing);
  if (Verbosity() > 2)
  {
    std::cout << "Ndivisions: " << Ndivisions << std::endl;
  }
  for (int i = 1; i <= Ndivisions + 1; i++)
  {
    if (std::isnan(kftrack.GetX()) ||
        std::isnan(kftrack.GetY()) ||
        std::isnan(kftrack.GetZ()))
    {
      if (Verbosity() > 1)
      {
        std::cout << "position is NaN, exiting" << std::endl;
      }
      return false;
    }

    if (new_radius > 78.)
    {
      if (Verbosity() > 1)
      {
        std::cout << "outside TPC, exiting" << std::endl;
      }
      return false;
    }

    double r_div;
    if (i == Ndivisions + 1)
    {
      r_div = new_radius;
    }
    else if (old_radius < new_radius)
    {
      r_div = old_radius + transport_spacing * i;
    }
    else
    {
      r_div = old_radius - transport_spacing * i;
    }

    if (Verbosity() > 2)
    {
      std::cout << "transporting to intermediate radius " << r_div << std::endl;
    }

    // track state position relative to radial direction
    const double tX = kftrack.GetX();
    const double tY = kftrack.GetY();
    const double tz = kftrack.GetZ();

    // track state global position (including tz above)
    const double tx = tX * cos(phi) - tY * sin(phi);
    const double ty = tX * sin(phi) + tY * cos(phi);

    const double Bz = _Bzconst * get_Bz(tx, ty, tz);

    kftrack.CalculateFitParameters(fp);

    // transport to radius
    if (!kftrack.TransportToXWithMaterial(r_div, fp, Bz, 1.))
    {
      if (Verbosity() > 1)
      {
        std::cout << "transport failed" << std::endl;
      }
      return false;
    }
    // rotate track state reference frame
    const double new_tX = kftrack.GetX();
    const double new_tY = kftrack.GetY();

    const double new_tx = new_tX * cos(phi) - new_tY * sin(phi);
    const double new_ty = new_tX * sin(phi) + new_tY * cos(phi);
    const double new_phi = atan2(new_ty, new_tx);
    const double alpha = new_phi - phi;

    if (!kftrack.Rotate(alpha, 1.))
    {
      if (Verbosity() > 1)
      {
        std::cout << "rotate failed" << std::endl;
      }
      return false;
    }
    phi = new_phi;
  }

  if (Verbosity() > 1)
  {
    const double final_tX = kftrack.GetX();
    const double final_tY = kftrack.GetY();
    const double final_tx = final_tX * cos(phi) - final_tY * sin(phi);
    const double final_ty = final_tX * sin(phi) + final_tY * cos(phi);
    const double final_tz = kftrack.GetZ();
    std::cout << "track position after transport: (" << final_tx << ", " << final_ty << ", " << final_tz << ")" << std::endl;
  }
  return true;
}

//___________________________________________________________________________________________
bool PHSimpleKFProp::PropagateStep(
  unsigned int& current_layer,
  double& current_phi,
  PropagationDirection& direction,
  std::vector<TrkrDefs::cluskey>& propagated_track,
  const std::vector<TrkrDefs::cluskey>& ckeys,
  GPUTPCTrackParam& kftrack,
  GPUTPCTrackParam::GPUTPCTrackFitParam& fp,
  const PositionMap& globalPositions) const
{
  // give up if position vector is NaN (propagation failed)
  if (std::isnan(kftrack.GetX()) ||
      std::isnan(kftrack.GetY()) ||
      std::isnan(kftrack.GetZ()))
  {
    if (Verbosity() > 1)
    {
      std::cout << "position is NaN, exiting loop" << std::endl;
    }
    return false;
  }

  if (Verbosity() > 1)
  {
    std::cout << "current layer: " << current_layer << std::endl;
  }
  if (Verbosity() > 5)
  {
    std::cout << "original seed size: " << ckeys.size() << std::endl;
  }

  // CosPhi is the projection of the pt onto the radial direction
  // based on the sign of CosPhi, decide whether to propagate outward or inward
  int next_layer;
  if (direction == PropagationDirection::Outward)
  {
    next_layer = current_layer + 1;
  }
  else
  {
    next_layer = current_layer - 1;
  }

  if (Verbosity() > 1)
  {
    std::cout << "next layer: " << next_layer << std::endl;
  }

  // give up if next layer is outside the TPC (propagation complete)
  if (next_layer < 7 || next_layer > 54)
  {
    if (Verbosity() > 1)
    {
      std::cout << "reached end of TPC, exiting loop" << std::endl;
    }
    return false;
  }

  // track state position relative to radial direction
  const double tX = kftrack.GetX();
  const double tY = kftrack.GetY();
  const double tz = kftrack.GetZ();

  // track state global position (including tz above)
  const double tx = tX * cos(current_phi) - tY * sin(current_phi);
  const double ty = tX * sin(current_phi) + tY * cos(current_phi);

  if (Verbosity() > 1)
  {
    std::cout << std::endl;
    std::cout << "track position: (" << tx << ", " << ty << ", " << tz << ")" << std::endl;
  }

  if (Verbosity() > 1)
  {
    std::cout << "transporting to radius: " << radii[next_layer - 7] << std::endl;
  }

  if (!TransportAndRotate(kftrack.GetX(), radii[next_layer - 7], current_phi, kftrack, fp))
  {
    if (std::isnan(kftrack.GetX()) ||
        std::isnan(kftrack.GetY()) ||
        std::isnan(kftrack.GetZ()))
    {
      return false;
    }

    if (Verbosity() > 1)
    {
      std::cout << "track turned around near X=" << kftrack.GetX() << std::endl;
    }

    // TransportAndRotate failure indicates that track has turned around before it reaches next layer
    // The track parameters don't update in that last step, so we're likely somewhere between two layers
    // So, first, we turn the track around ourselves, setting it to its next intersection at its current radius

    // basically circle project in xy, linear project in z
    double pt = 1. / std::abs(kftrack.GetQPt());
    double end_tx = kftrack.GetX() * cos(current_phi) - kftrack.GetY() * sin(current_phi);
    double end_ty = kftrack.GetX() * sin(current_phi) + kftrack.GetY() * cos(current_phi);
    double end_tz = kftrack.GetZ();
    if (Verbosity() > 1)
    {
      std::cout << "current parameters: pt=" << pt << ", (x, y, z) = (" << end_tx << ", " << end_ty << ", " << end_tz << ")" << std::endl;
    }
    // pt[GeV] = 0.3 B[T] R[m]
    double R = 100. * pt / (0.3 * get_Bz(end_tx, end_ty, end_tz));
    if (Verbosity() > 1)
    {
      std::cout << "R=" << R << std::endl;
    }
    double pX = pt * kftrack.GetCosPhi();
    double pY = pt * kftrack.GetSinPhi();
    if (Verbosity() > 1)
    {
      std::cout << "(pX, pY) = (" << pX << ", " << pY << ")" << std::endl;
    }
    double px = pX * cos(current_phi) - pY * sin(current_phi);
    double py = pX * sin(current_phi) + pY * cos(current_phi);
    double tangent_phi = std::atan2(py, px);
    double center_phi = 0;
    if (kftrack.GetQPt() > 0)
    {
      center_phi = tangent_phi + M_PI / 2.;
    }
    else
    {
      center_phi = tangent_phi - M_PI / 2.;
    }
    if (center_phi > M_PI)
    {
      center_phi -= 2. * M_PI;
    }
    if (center_phi < -M_PI)
    {
      center_phi += 2. * M_PI;
    }
    double xc = end_tx - R * cos(center_phi);
    double yc = end_ty - R * sin(center_phi);
    if (Verbosity() > 1)
    {
      std::cout << "(xc, yc) = (" << xc << ", " << yc << ")" << std::endl;
    }
    auto circle_output = TrackFitUtils::circle_circle_intersection(
      std::sqrt(square(kftrack.GetX()) + square(kftrack.GetY())), R, xc, yc);
    // pick the furthest point from current track position
    double new_tx;
    double new_ty;
    double xplus = std::get<0>(circle_output);
    double yplus = std::get<1>(circle_output);
    double xminus = std::get<2>(circle_output);
    double yminus = std::get<3>(circle_output);
    if (Verbosity() > 1)
    {
      std::cout << "circle-circle intersection: (" << xplus << ", " << yplus << "), (" << xminus << ", " << yminus << ")" << std::endl;
    }

    if (std::sqrt(square(end_tx - xplus) + square(end_ty - yplus)) > std::sqrt(square(end_tx - xminus) + square(end_ty - yminus)))
    {
      new_tx = xplus;
      new_ty = yplus;
    }
    else
    {
      new_tx = xminus;
      new_ty = yminus;
    }

    if (Verbosity() > 1)
    {
      std::cout << "track now at (" << new_tx << ", " << new_ty << ")" << std::endl;
    }
    const double rot_phi = atan2(new_ty, new_tx);
    // double rot_alpha = rot_phi - current_phi;

    // new track point is existing track point rotated by alpha
    const double new_tX = new_tx * cos(rot_phi) + new_ty * sin(rot_phi);
    const double new_tY = -new_tx * sin(rot_phi) + new_ty * cos(rot_phi);
    const double new_centerphi = atan2(new_ty - yc, new_tx - xc);
    double dcenterphi = new_centerphi - center_phi;
    if (dcenterphi > M_PI)
    {
      dcenterphi = 2. * M_PI - dcenterphi;
    }
    if (dcenterphi < -M_PI)
    {
      dcenterphi = 2. * M_PI + dcenterphi;
    }
    const double ds = R * std::abs(dcenterphi);
    const double dz = kftrack.GetDzDs() * ds;

    current_phi = rot_phi;
    kftrack.SetX(new_tX);
    kftrack.SetY(new_tY);
    kftrack.SetZ(end_tz + dz);

    kftrack.SetSignCosPhi(-kftrack.GetSignCosPhi());

    // now finish transport back down to current layer
    if (!TransportAndRotate(kftrack.GetX(), radii[current_layer - 7], current_phi, kftrack, fp))
    {
      return false;
    }

    if (direction == PropagationDirection::Outward)
    {
      direction = PropagationDirection::Inward;
    }
    else
    {
      direction = PropagationDirection::Outward;
    }

    // track landed in same layer for this step
    next_layer = current_layer;
  }

  // account for distortions
  if(!_pp_mode)
  {
    if (Verbosity() > 1)
    {
      std::cout << "doing distortion correction" << std::endl;
    }

    // The distortion corrected cluster positions in globalPos are not at the layer radius
    // We want to project to the radius appropriate for the globalPos values
    // Get the distortion correction for the projection point, and calculate the radial increment

    const double uncorr_tX = kftrack.GetX();
    const double uncorr_tY = kftrack.GetY();
    const double uncorr_tx = uncorr_tX * cos(current_phi) - uncorr_tY * sin(current_phi);
    const double uncorr_ty = uncorr_tX * sin(current_phi) + uncorr_tY * cos(current_phi);
    const double uncorr_tz = kftrack.GetZ();

    Acts::Vector3 proj_pt(uncorr_tx, uncorr_ty, uncorr_tz);

    const double proj_radius = sqrt(square(uncorr_tx) + square(uncorr_ty));
    if (proj_radius > 78.0)
    {
      // projection is bad, no cluster will be found
      return false;
    }

    if (Verbosity() > 2)
    {
      std::cout << "distortion correction for (" << tx << ", " << ty << ", " << tz << "), layer " << next_layer << ", radius " << proj_radius << std::endl;
    }

    // apply distortion corrections
    proj_pt = m_globalPositionWrapper.applyDistortionCorrections(proj_pt);

    // calculate radius
    const double dist_radius = sqrt(square(proj_pt[0]) + square(proj_pt[1]));

    if (Verbosity() > 2)
    {
      std::cout << "after correction: (" << proj_pt[0] << ", " << proj_pt[1] << ", " << proj_pt[2] << "), radius " << dist_radius << std::endl;
    }

    if (!TransportAndRotate(kftrack.GetX(), dist_radius, current_phi, kftrack, fp))
    {
      return false;
    }

    if (std::isnan(kftrack.GetX()) ||
        std::isnan(kftrack.GetY()) ||
        std::isnan(kftrack.GetZ()))
    {
      return false;
    }
  }

  const double new_tX = kftrack.GetX();
  const double new_tY = kftrack.GetY();
  const double new_tx = new_tX * cos(current_phi) - new_tY * sin(current_phi);
  const double new_ty = new_tX * sin(current_phi) + new_tY * cos(current_phi);
  const double new_tz = kftrack.GetZ();

  // search for closest available cluster within window
  double query_pt[3] = {new_tx, new_ty, new_tz};
  std::vector<long unsigned int> index_out(1);
  std::vector<double> distance_out(1);
  int n_results = _kdtrees[next_layer]->knnSearch(&query_pt[0], 1, &index_out[0], &distance_out[0]);

  // if no results, then no cluster to add, but propagation is not necessarily done
  if (!n_results)
  {
    if (Verbosity() > 1)
    {
      std::cout << "no clusters found in search window, moving on" << std::endl;
    }
    current_layer = next_layer;
    return true;
  }
  const std::vector<double>& point = _ptclouds[next_layer]->pts[index_out[0]];
  TrkrDefs::cluskey closest_ckey = (*((int64_t*) &point[3]));
  TrkrCluster* clusterCandidate = _cluster_map->findCluster(closest_ckey);
  const auto candidate_globalpos = globalPositions.at(closest_ckey);
  const double cand_x = candidate_globalpos(0);
  const double cand_y = candidate_globalpos(1);
  const double cand_z = candidate_globalpos(2);

  if (Verbosity() > 1)
  {
    std::cout << "found closest cluster candidate at (" << cand_x << ", " << cand_y << ", " << cand_z << ")" << std::endl;
  }

  // get cluster and track state position errors
  const double tYerr = std::sqrt(kftrack.GetCov(0));
  const double txerr = std::abs(tYerr * sin(current_phi));
  const double tyerr = std::abs(tYerr * cos(current_phi));
  const double tzerr = sqrt(kftrack.GetCov(5));

  if (Verbosity() > 1)
  {
    std::cout << "track state position errors: (" << txerr << ", " << tyerr << ", " << tzerr << ")" << std::endl;
  }

  const double cand_xerr = sqrt(fitter->getClusterError(clusterCandidate, closest_ckey, candidate_globalpos, 0, 0));
  const double cand_yerr = sqrt(fitter->getClusterError(clusterCandidate, closest_ckey, candidate_globalpos, 1, 1));
  const double cand_zerr = sqrt(fitter->getClusterError(clusterCandidate, closest_ckey, candidate_globalpos, 2, 2));

  if (Verbosity() > 1)
  {
    std::cout << "cluster position errors: (" << cand_xerr << ", " << cand_yerr << ", " << cand_zerr << ")" << std::endl;
  }

  // add cluster if its position is within error-based window
  if (std::abs(new_tx - cand_x) < _max_dist * std::sqrt(square(txerr) + square(cand_xerr)) &&
      std::abs(new_ty - cand_y) < _max_dist * std::sqrt(square(tyerr) + square(cand_yerr)) &&
      std::abs(new_tz - cand_z) < _max_dist * std::sqrt(square(tzerr) + square(cand_zerr)))
  {
    if (Verbosity() > 1)
    {
      std::cout << "added cluster" << std::endl;
    }
    propagated_track.push_back(closest_ckey);

    // don't re-filter clusters that are already in original seed
    if (std::find(ckeys.begin(), ckeys.end(), closest_ckey) == ckeys.end())
    {
      const double cand_Y = -cand_x * std::sin(current_phi) + cand_y * std::cos(current_phi);
      const double cand_xycov2 = fitter->getClusterError(clusterCandidate, closest_ckey, candidate_globalpos, 0, 1);
      const double cand_Yerr2 = cand_xerr * cand_xerr * sin(current_phi) * sin(current_phi) + cand_xycov2 * sin(current_phi) * cos(current_phi) + cand_yerr * cand_yerr * cos(current_phi) * cos(current_phi);
      const double cand_zerr2 = cand_zerr * cand_zerr;

      if (Verbosity() > 1)
      {
        std::cout << "Filtering cluster with Y=" << cand_Y << ", z=" << cand_z << ", Yerr2=" << cand_Yerr2 << ", zerr2=" << cand_zerr2 << std::endl;
      }

      kftrack.Filter(cand_Y, cand_z, cand_Yerr2, cand_zerr2, _max_sin_phi);
    }
  }

  // update layer
  current_layer = next_layer;

  if (Verbosity() > 1)
  {
    std::cout << "track current parameters:" << std::endl;
    std::cout << "(X, Y, Z) = (" << kftrack.GetX() << ", " << kftrack.GetY() << ", " << kftrack.GetZ() << ")" << std::endl;
    std::cout << "pt = " << 1. / fabs(kftrack.GetQPt()) << std::endl;
    std::cout << "sinPhi = " << kftrack.GetSinPhi() << std::endl;
    std::cout << "cosPhi = " << kftrack.GetCosPhi() << std::endl;
    std::cout << "dzds = " << kftrack.GetDzDs() << std::endl;
  }

  // propagation step done
  return true;
}

//_________________________________________________________________
std::vector<TrkrDefs::cluskey> PHSimpleKFProp::PropagateTrack(TrackSeed* track, PropagationDirection direction, GPUTPCTrackParam& aliceSeed, const PositionMap& globalPositions) const
{
  // extract cluster list
  std::vector<TrkrDefs::cluskey> ckeys;
  if (direction == PropagationDirection::Inward)
  {
    std::reverse_copy(track->begin_cluster_keys(), track->end_cluster_keys(), std::back_inserter(ckeys));
  } else {
    std::copy(track->begin_cluster_keys(), track->end_cluster_keys(), std::back_inserter(ckeys));
  }
  return PropagateTrack(track, ckeys, direction, aliceSeed, globalPositions);
}

//_________________________________________________________________
std::vector<TrkrDefs::cluskey> PHSimpleKFProp::PropagateTrack(TrackSeed* track, std::vector<TrkrDefs::cluskey>& ckeys, PropagationDirection direction, GPUTPCTrackParam& aliceSeed, const PositionMap& globalPositions) const
{
  if (direction == PropagationDirection::Inward)
  {
    aliceSeed.SetDzDs(-aliceSeed.GetDzDs());
  }

  const auto track_pos = TrackSeedHelper::get_xyz(track);
  double track_x = track_pos.x();
  double track_y = track_pos.y();
  double track_z = track_pos.z();

  if (Verbosity() > 1)
  {
    std::cout << "layers:" << std::endl;
    for (auto& c : ckeys)
    {
      std::cout << (int) TrkrDefs::getLayer(c) << ", ";
    }
    std::cout << std::endl;
  }

  if (Verbosity() > 1)
  {
    std::cout << "track (x,y,z) = (" << track_x << ", " << track_y << ", " << track_z << ")" << std::endl;
  }

  double track_px = NAN;
  double track_py = NAN;
  double track_pz = NAN;
  if (_use_const_field)
  {
    float pt = fabs(1. / track->get_qOverR()) * (0.3 / 100) * _const_field;
    float phi = track->get_phi();
    track_px = pt * std::cos(phi);
    track_py = pt * std::sin(phi);
    track_pz = pt * std::cosh(track->get_eta()) * std::cos(track->get_theta());
  }
  else
  {
    track_px = track->get_px();
    track_py = track->get_py();
    track_pz = track->get_pz();
  }

  if (Verbosity() > 1)
  {
    std::cout << "track (px,py,pz) = (" << track_px << ", " << track_py << ", " << track_pz << ")" << std::endl;
  }

  std::vector<Acts::Vector3> trkGlobPos;
  for (const auto& ckey : ckeys)
  {
    if (TrkrDefs::getTrkrId(ckey) == TrkrDefs::tpcId)
    {
      trkGlobPos.push_back(globalPositions.at(ckey));
    }
  }

  // want angle of tangent to circle at innermost (i.e. last) cluster
  size_t inner_index;
  size_t second_index;
  if (TrkrDefs::getLayer(ckeys[0]) > TrkrDefs::getLayer(ckeys.back()))
  {
    inner_index = ckeys.size() - 1;
    second_index = ckeys.size() - 2;
  }
  else
  {
    inner_index = 0;
    second_index = 1;
  }

  double xc = track->get_X0();
  double yc = track->get_Y0();
  double cluster_x = trkGlobPos.at(inner_index)(0);
  double cluster_y = trkGlobPos.at(inner_index)(1);
  double dy = cluster_y - yc;
  double dx = cluster_x - xc;
  double phi = atan2(dy, dx);
  double second_dx = trkGlobPos.at(second_index)(0) - xc;
  double second_dy = trkGlobPos.at(second_index)(1) - yc;
  double second_phi = atan2(second_dy, second_dx);
  double dphi = second_phi - phi;

  if (dphi > 0)
  {
    phi += M_PI / 2.0;
  }
  else
  {
    phi -= M_PI / 2.0;
  }

  double pt = sqrt(track_px * track_px + track_py * track_py);
  // rotate track momentum vector (pz stays the same)
  track_px = pt * cos(phi);
  track_py = pt * sin(phi);
  track_x = trkGlobPos.at(0)(0);
  track_y = trkGlobPos.at(0)(1);
  track_z = trkGlobPos.at(0)(2);

  if (Verbosity() > 1)
  {
    std::cout << "new track (x,y,z) = " << track_x << ", " << track_y << ", " << track_z << ")" << std::endl;
  }
  if (Verbosity() > 1)
  {
    std::cout << "new track (px,py,pz) = " << track_px << ", " << track_py << ", " << track_pz << ")" << std::endl;
  }

  //    }

  double track_pt = sqrt(track_px * track_px + track_py * track_py);
  if (Verbosity() > 1)
  {
    std::cout << "track pt: " << track_pt << std::endl;
  }

  // get track parameters
  GPUTPCTrackParam kftrack{};
  kftrack.InitParam();
  kftrack.setNeonFraction(Ne_frac);
  kftrack.setArgonFraction(Ar_frac);
  kftrack.setCF4Fraction(CF4_frac);
  kftrack.setNitrogenFraction(N2_frac);
  kftrack.setIsobutaneFraction(isobutane_frac);

  float track_phi = atan2(track_y, track_x);
  if (Verbosity() > 1)
  {
    std::cout << "track phi: " << track_phi << std::endl;
    std::cout << "track charge: " << track->get_charge() << std::endl;
  }
  kftrack.SetQPt(track->get_charge() / track_pt);

  float track_pX = track_px * cos(track_phi) + track_py * sin(track_phi);
  float track_pY = -track_px * sin(track_phi) + track_py * cos(track_phi);

  kftrack.SetSignCosPhi(track_pX / track_pt);
  kftrack.SetSinPhi(track_pY / track_pt);
  kftrack.SetDzDs(track_pz / sqrt(track_pt * track_pt + track_pz * track_pz));

  if (Verbosity() > 1)
  {
    std::cout << "track pX = " << track_pX << ", pY = " << track_pY << ", CosPhi = " << track_pX / track_pt << ", signCosPhi = " << kftrack.GetSignCosPhi() << std::endl;
  }
  /*
    // Y = y
    // Z = z
    // SinPhi = py/sqrt(px^2+py^2)
    // DzDs = pz/sqrt(px^2+py^2)
    // QPt = 1/sqrt(px^2+py^2)

    const double track_pt3 = std::pow(track_pt, 3.);

    Eigen::Matrix<double, 6, 5> Jrot;
    Jrot(0, 0) = 0;  // dY/dx
    Jrot(1, 0) = 1;  // dY/dy
    Jrot(2, 0) = 0;  // dY/dz
    Jrot(3, 0) = 0;  // dY/dpx
    Jrot(4, 0) = 0;  // dY/dpy
    Jrot(5, 0) = 0;  // dY/dpz

    Jrot(0, 1) = 0;  // dZ/dx
    Jrot(1, 1) = 0;  // dZ/dy
    Jrot(2, 1) = 1;  // dZ/dz
    Jrot(3, 1) = 0;  // dZ/dpx
    Jrot(4, 1) = 0;  // dZ/dpy
    Jrot(5, 1) = 0;  // dZ/dpz

    Jrot(0, 2) = 0;                                 // dSinPhi/dx
    Jrot(1, 2) = 0;                                 // dSinPhi/dy
    Jrot(2, 2) = 0;                                 // dSinPhi/dz
    Jrot(3, 2) = -track_py * track_px / track_pt3;  // dSinPhi/dpx
    Jrot(4, 2) = track_px * track_px / track_pt3;   // dSinPhi/dpy
    Jrot(5, 2) = 0;                                 // dSinPhi/dpz

    Jrot(0, 3) = 0;                                 // dDzDs/dx
    Jrot(1, 3) = 0;                                 // dDzDs/dy
    Jrot(2, 3) = 0;                                 // dDzDs/dz
    Jrot(3, 3) = -track_px * track_pz / track_pt3;  // dDzDs/dpx
    Jrot(4, 3) = -track_py * track_pz / track_pt3;  // dDzDs/dpy
    Jrot(5, 3) = 1. / track_pt;                     // dDzDs/dpz

    Jrot(0, 4) = 0;                      // dQPt/dx
    Jrot(1, 4) = 0;                      // dQPt/dy
    Jrot(2, 4) = 0;                      // dQPt/dz
    Jrot(3, 4) = -track_px / track_pt3;  // dQPt/dpx
    Jrot(4, 4) = -track_py / track_pt3;  // dQPt/dpy
    Jrot(5, 4) = 0;                      // dQPt/dpz

    Eigen::Matrix<double, 5, 5> kfCov = Jrot.transpose() * xyzCov * Jrot;

    int ctr = 0;
    for (int i = 0; i < 5; i++)
    {
      for (int j = 0; j < 5; j++)
      {
        if (i >= j)
        {
          kftrack.SetCov(ctr, kfCov(i, j));
          ctr++;
        }
      }
    }
  */
  std::vector<TrkrDefs::cluskey> propagated_track;

  kftrack.SetX(track_x * cos(track_phi) + track_y * sin(track_phi));
  kftrack.SetY(-track_x * sin(track_phi) + track_y * cos(track_phi));
  kftrack.SetZ(track_z);

  if (kftrack.GetSignCosPhi() < 0)
  {
    kftrack.SetSinPhi(-kftrack.GetSinPhi());
    kftrack.SetDzDs(-kftrack.GetDzDs());
    // kftrack.SetQPt(-kftrack.GetQPt());
    kftrack.SetCov(3, -kftrack.GetCov(3));
    kftrack.SetCov(4, -kftrack.GetCov(4));
    kftrack.SetCov(6, -kftrack.GetCov(6));
    kftrack.SetCov(7, -kftrack.GetCov(7));
    kftrack.SetCov(10, -kftrack.GetCov(10));
    kftrack.SetCov(11, -kftrack.GetCov(11));
  }

  if (Verbosity() > 1)
  {
    std::cout << "initial track params:" << std::endl;
    std::cout << "X: " << kftrack.GetX() << std::endl;
    std::cout << "Y: " << kftrack.GetY() << std::endl;
    std::cout << "Z: " << kftrack.GetZ() << std::endl;
    std::cout << "SinPhi: " << kftrack.GetSinPhi() << std::endl;
    std::cout << "DzDs: " << kftrack.GetDzDs() << std::endl;
    std::cout << "QPt: " << kftrack.GetQPt() << std::endl;
    std::cout << "cov: " << std::endl;
    for (int i = 0; i < 15; i++)
    {
      std::cout << kftrack.GetCov(i) << ", ";
    }
    std::cout << std::endl;
  }

  // get layer for each cluster
  std::vector<unsigned int> layers;
  std::transform(ckeys.begin(), ckeys.end(), std::back_inserter(layers), [](const TrkrDefs::cluskey& key)
                 { return TrkrDefs::getLayer(key); });

  double old_phi = track_phi;
  unsigned int old_layer = TrkrDefs::getLayer(ckeys[0]);
  if (Verbosity() > 1)
  {
    std::cout << "first layer: " << old_layer << std::endl;
    std::cout << "cluster (x,y,z) = (" << globalPositions.at(ckeys[0])(0) << ", " << globalPositions.at(ckeys[0])(1) << ", " << globalPositions.at(ckeys[0])(2) << ")" << std::endl;
  }

  propagated_track.push_back(ckeys[0]);

  GPUTPCTrackLinearisation kfline(aliceSeed);
  GPUTPCTrackParam::GPUTPCTrackFitParam fp{};
  aliceSeed.CalculateFitParameters(fp);

  for (int step = 0; step <= _max_propagation_steps; ++step)
  {
    if (Verbosity() > 1)
    {
      std::cout << std::endl
                << "------------------------" << std::endl
                << "step " << step << std::endl;
    }
    if (!PropagateStep(old_layer, old_phi, direction, propagated_track, ckeys, aliceSeed, fp, globalPositions))
    {
      break;
    }
  }

  if (Verbosity() > 1)
  {
    std::cout << "propagated track has " << propagated_track.size() << " clusters, at layers: ";
    for (auto c : propagated_track)
    {
      std::cout << (int) TrkrDefs::getLayer(c) << ", ";
    }
    std::cout << std::endl;
  }
  return propagated_track;
}

std::vector<keylist> PHSimpleKFProp::RemoveBadClusters(const std::vector<keylist>& chains, const PositionMap& /* globalPositions */) const
{
  if (Verbosity())
  {
    std::cout << "removing bad clusters" << std::endl;
  }
  std::vector<keylist> clean_chains;
  /*
    Hugo: this whole code just copy chains of size >= 3 to the output,
    because of the cut on outliers being commented out.
    Replaced with a simpler version
  */
  //   for(const keylist& chain : chains)
  //   {
  //     if(chain.size()<3) continue;
  //
  //     keylist clean_chain;
  //
  //
  //     std::vector<std::pair<double,double>> xy_pts;
  //     std::vector<std::pair<double,double>> rz_pts;
  //     for(const TrkrDefs::cluskey& ckey : chain)
  //     {
  //       auto global = globalPositions.at(ckey);
  //       xy_pts.push_back(std::make_pair(global(0),global(1)));
  //       float r = sqrt(square( global(0) ) + square( global(1) ));
  //       rz_pts.emplace_back( r,global(2) );
  //     }
  //     if(Verbosity()) std::cout << "chain size: " << chain.size() << std::endl;
  //
  //     //fit a circle through x,y coordinates and calculate residuals
  //     const auto [R, X0, Y0] = TrackFitUtils::circle_fit_by_taubin( xy_pts );
  //     const std::vector<double> xy_resid = TrackFitUtils::getCircleClusterResiduals(xy_pts,R,X0,Y0);
  //
  //     // fit a line through r,z coordinates and calculate residuals
  //     const auto [A, B] = TrackFitUtils::line_fit( rz_pts );
  //     const std::vector<double> rz_resid = TrackFitUtils::getLineClusterResiduals(rz_pts,A,B);
  //
  //     for(size_t i=0;i<chain.size();i++)
  //     {
  //       //if(xy_resid[i]>_xy_outlier_threshold) continue;
  //       clean_chain.push_back(chain[i]);
  //     }
  //
  //     clean_chains.push_back(clean_chain);
  //     if(Verbosity()) std::cout << "pushed clean chain with " << clean_chain.size() << " clusters" << std::endl;
  //   }

  // simpler version
  std::copy_if(chains.begin(), chains.end(), std::back_inserter(clean_chains),
               [](const keylist& chain)
               { return chain.size() >= 3; });

  return clean_chains;
}

void PHSimpleKFProp::rejectAndPublishSeeds(std::vector<TrackSeed_v2>& seeds, const PositionMap& positions, std::vector<float>& trackChi2)
{

  PHTimer timer("KFPropTimer");

  // now do the ghost rejection *before* publishing the seeds to the _track_map
  timer.restart();

  // testing with presets for rejection
  PHGhostRejection rejector(Verbosity(), seeds);
  rejector.set_phi_cut(_ghost_phi_cut);
  rejector.set_eta_cut(_ghost_eta_cut);
  rejector.set_x_cut(_ghost_x_cut);
  rejector.set_y_cut(_ghost_y_cut);
  rejector.set_z_cut(_ghost_z_cut);
  // If you want to reject tracks (before they are are made) can set them here:
  // rejector.set_min_pt_cut(0.2);
  // rejector.set_must_span_sectors(true);
  // rejector.set_min_clusters(8);

  // first path over seeds to marks those to be removed due to cluster size
  for (unsigned int itrack = 0; itrack < seeds.size(); ++itrack)
  { rejector.cut_from_clusters(itrack); }

  #pragma omp parallel
  {

    #pragma omp for schedule(static)
    for (unsigned int itrack = 0; itrack < seeds.size(); ++itrack)
    {
      // cut tracks with too-few clusters (or that don;t span a sector boundary, if desired)
      if (rejector.is_rejected(itrack))
      { continue; }

      auto& seed = seeds[itrack];
      /// The ALICEKF gives a better charge determination at high pT
      const int q = seed.get_charge();

      PositionMap local;
      std::transform(seed.begin_cluster_keys(), seed.end_cluster_keys(), std::inserter(local, local.end()),
        [positions](const auto& key)
        { return std::make_pair(key, positions.at(key)); });
      TrackSeedHelper::circleFitByTaubin(&seed,local, 7, 55);
      TrackSeedHelper::lineFit(&seed,local, 7, 55);
      seed.set_phi(TrackSeedHelper::get_phi(&seed,local));
      seed.set_qOverR(fabs(seed.get_qOverR()) * q);
    }

  }

  if (Verbosity())
  { std::cout << "PHSimpleKFProp::rejectAndPublishSeeds - circle fit: " << timer.elapsed() << " ms" << std::endl; }

  // now do the ghost rejection *before* publishing the seeds to the _track_map
  timer.restart();
  rejector.find_ghosts(trackChi2);
  if (Verbosity())
  { std::cout << "PHSimpleKFProp::rejectAndPublishSeeds - ghost rejection: " << timer.elapsed() << " ms" << std::endl; }

  timer.restart();
  for (unsigned int itrack = 0; itrack < seeds.size(); ++itrack)
  {
    if (rejector.is_rejected(itrack))
    {
      if (Verbosity() > 0)
      {
        std::cout << " Seed " << ((int) itrack) << " rejected. Not getting published." << std::endl;
      }
      continue;
    }
    auto& seed = seeds[itrack];
    _track_map->insert(&seed);

    int q = seed.get_charge();
    if (Verbosity() > 0)
    {
      std::cout << "Publishing seed " << ((int) itrack)
                << " q " << q
                << " qOverR " << fabs(seed.get_qOverR()) * q
                << " x " << TrackSeedHelper::get_x(&seed)
                << " y " << TrackSeedHelper::get_y(&seed)
                << " z " << TrackSeedHelper::get_z(&seed)
                << " pT " << seed.get_pt()
                << " eta " << seed.get_eta()
                << " phi " << seed.get_phi()
                << std::endl;
    }
  }
  if (Verbosity())
  { std::cout << "PHSimpleKFProp::rejectAndPublishSeeds - publication: " << timer.elapsed() << " ms" << std::endl; }

}

void PHSimpleKFProp::publishSeeds(const std::vector<TrackSeed_v2>& seeds)
{
  for (const auto& seed : seeds)
  {
    _track_map->insert(&seed);
  }
}
