/*!
 *  \File PrelimDistortionCorrection.cc
 *  \brief	Makes preliminary TPC distortion corrections in pp running and replaces TPC seed parameters with new fits
 *  \author Tony Frawley
 */

#include "PrelimDistortionCorrection.h"

#include "ALICEKF.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <phfield/PHField.h>
#include <phfield/PHFieldUtility.h>
#include <phfield/PHFieldConfig.h>
#include <phfield/PHFieldConfigv1.h>

#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>                       // for PHWHERE

// tpc distortion correction
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsGeometry.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <Geant4/G4SystemOfUnits.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>                            // for operator<<, basic_ostream
#include <vector>

//___________________________________________________________________________________________
PrelimDistortionCorrection::PrelimDistortionCorrection(const std::string& name)
  : SubsysReco(name)
{}

//___________________________________________________________________________________________
int PrelimDistortionCorrection::End(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________________
int PrelimDistortionCorrection::InitRun(PHCompositeNode* topNode)
{

  int ret = get_nodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) { return ret; }

  PHFieldConfigv1 fcfg;
  fcfg.set_field_config(PHFieldConfig::FieldConfigTypes::Field3DCartesian);

  char *calibrationsroot = getenv("CALIBRATIONROOT");
  assert(calibrationsroot);

  auto magField = std::string(calibrationsroot) +
    std::string("/Field/Map/sphenix3dtrackingmapxyz.root");

  fcfg.set_filename(magField);

  //  fcfg.set_rescale(1);
  _field_map = std::unique_ptr<PHField>(PHFieldUtility::BuildFieldMap(&fcfg));

  fitter = std::make_unique<ALICEKF>(topNode,_cluster_map,_field_map.get(), _fieldDir,
				     _min_clusters_per_track,_max_sin_phi,Verbosity());
  fitter->setNeonFraction(Ne_frac);
  fitter->setArgonFraction(Ar_frac);
  fitter->setCF4Fraction(CF4_frac);
  fitter->setNitrogenFraction(N2_frac);
  fitter->setIsobutaneFraction(isobutane_frac);
  fitter->useConstBField(_use_const_field);
  fitter->setConstBField(_const_field);
  fitter->useFixedClusterError(_use_fixed_clus_err);
  fitter->setFixedClusterError(0,_fixed_clus_err.at(0));
  fitter->setFixedClusterError(1,_fixed_clus_err.at(1));
  fitter->setFixedClusterError(2,_fixed_clus_err.at(2));
  //  _field_map = PHFieldUtility::GetFieldMapNode(nullptr,topNode);
  // m_Cache = magField->makeCache(m_tGeometry->magFieldContext);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________
int PrelimDistortionCorrection::get_nodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  _tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if(!_tgeometry)
    {
      std::cout << "No Acts tracking geometry, exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  // tpc distortion correction
  // input tpc distortion correction module edge
  m_dcc_module_edge = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerModuleEdge");
  if (m_dcc_module_edge)
  {
    std::cout << "PrelimDistortionCorrection::get_nodes - found TPC distortion correction container module edge" << std::endl;
  }

  // input tpc distortion correction static
  m_dcc_static = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerStatic");
  if (m_dcc_static)
  {
    std::cout << "PrelimDistortionCorrection::get_nodes - found TPC distortion correction container static" << std::endl;
  }

  // input tpc distortion correction average
  m_dcc_average = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerAverage");
  if (m_dcc_average)
  {
    std::cout << "PrelimDistortionCorrection::get_nodes - found TPC distortion correction container average" << std::endl;
  }


  if(_use_truth_clusters) {
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_TRUTH");
  } else {
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
}

  if (!_cluster_map)
  {
    std::cerr << PHWHERE << "PrelimDistortionCorrection::get_nodes - ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!_track_map)
  {
    std::cerr << PHWHERE << "PrelimDistortionCorrection::get_nodes - ERROR: Can't find TrackSeedContainer " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  PHG4TpcCylinderGeomContainer *geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cerr << PHWHERE << "PrelimDistortionCorrection::get_nodes - ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  for(int i=7;i<=54;i++)
  {
    radii.push_back(geom_container->GetLayerCellGeom(i)->get_radius());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//________________________________________________________________________________________________________________
int PrelimDistortionCorrection::process_event(PHCompositeNode* /*topNode*/)
{
  if(!_pp_mode)
  {
    std::cout << "PrelimDistortionCorrection called but not in pp_mode, do nothing!" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  PHTimer timer("PrelimDistortionCorrectionTimer");
  timer.stop();
  timer.restart();

  if(Verbosity()>0)
  {
    std::cout << "starting PrelimDistortionCorrection process_event" << std::endl;
  }

  // These accumulate trackseeds as we loop over tracks, keylist is a vector of vectors of seed cluskeys
  // The seeds are all given to the fitter at once
  std::vector<std::vector<TrkrDefs::cluskey>> keylist_A;
  PositionMap correctedOffsetTrackClusPositions;   //  this is an std::map<TrkrDefs::cluskey, Acts::Vector3>

  for(unsigned int track_it = 0; track_it != _track_map->size(); ++track_it )
  {
    if(Verbosity()>0) { std::cout << "TPC seed " << track_it << std::endl; }
    // if not a TPC track, ignore
    TrackSeed* track = _track_map->get(track_it);
    if(!track)
      {
	continue;
      }

    if(Verbosity() > 0)
      {
	std::cout << "Input seed pars for " << track_it
		  << " q " << track->get_charge()
		  << " qOverR " << fabs(track->get_qOverR()) * track->get_charge()
		  << " X0 " << TrackSeedHelper::get_x(track)
		  << " Y0 " << TrackSeedHelper::get_y(track)
		  << " Z0 " << TrackSeedHelper::get_z(track)
		  << " eta " << track->get_eta()
		  << " phi " << track->get_phi()
		  << std::endl;
      }

    const bool is_tpc = std::any_of(
      track->begin_cluster_keys(),
      track->end_cluster_keys(),
      []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getTrkrId(key) == TrkrDefs::tpcId; } );

    if(is_tpc)
    {
      // We want to move all clusters in this seed to point to Z0 = 0
      float Z0 = track->get_Z0();
      float offset_Z = 0.0 - Z0;

      if(Verbosity() > 0)
	{
	  std::cout << "  processing seed with offset_Z " << offset_Z
		    << " eta " << track->get_eta()
        << " x " << TrackSeedHelper::get_x(track)
        << " y " << TrackSeedHelper::get_y(track)
        << " z " << TrackSeedHelper::get_z(track)
		    << std::endl;
	}

      // We want to make  distortion corrections to all clusters in this seed after offsetting the z values
      std::vector<TrkrDefs::cluskey> dumvec;
      for(TrackSeed::ConstClusterKeyIter iter = track->begin_cluster_keys();
	  iter != track->end_cluster_keys();
	  ++iter)
	{
	  TrkrDefs::cluskey cluskey = *iter;
	  TrkrCluster *cluster = _cluster_map->findCluster(cluskey);

    // get global position
    const Acts::Vector3 pos = _tgeometry->getGlobalPosition(cluskey, cluster);

    // apply z offset
    Acts::Vector3 offsetpos(pos.x(), pos.y(), pos.z() + offset_Z);

    // apply distortion corrections
    if (m_dcc_module_edge)
    {
      offsetpos = m_distortionCorrection.get_corrected_position(offsetpos, m_dcc_module_edge);
    }

    if (m_dcc_static)
    {
      offsetpos = m_distortionCorrection.get_corrected_position(offsetpos, m_dcc_static);
    }

    if (m_dcc_average)
    {
      offsetpos = m_distortionCorrection.get_corrected_position(offsetpos, m_dcc_average);
    }

	  // now move the distortion corrected cluster back by offset_Z, to preserve the z measurement info
    const Acts::Vector3 corrpos(offsetpos(0), offsetpos(1), offsetpos(2) - offset_Z);
	  correctedOffsetTrackClusPositions.emplace(cluskey,corrpos);
	  dumvec.push_back(cluskey);

	  if(Verbosity() > 0)
	    {
	      std::cout << " cluskey " << cluskey << " input pos " << pos(0) << "  " << pos(1) << "  " << pos(2)
			<< "   corr. pos " << corrpos(0) << "  " << corrpos(1) << "  " << corrpos(2) << std::endl
			<< "distortion correction " <<  corrpos(0) - pos(0) << "  " << corrpos(1) - pos(1) << "  " << corrpos(2) - pos(2)
			<< std::endl;

	    }
	}

      /// Can't circle fit a seed with less than 3 clusters, skip it
      if(dumvec.size() < 3)
	{ continue; }
      keylist_A.push_back(dumvec);

      if(Verbosity() > 0)
	{
	  std::cout << "Added  input seed " << track_it << "  becomes output seed " << keylist_A.size() - 1 << std::endl;
	}

    } // end if TPC seed

  }  // end loop over tracks

  // reset the seed map for the TPC
  _track_map->Reset();

  // refit the corrected clusters
  std::vector<float> trackChi2;
  auto seeds = fitter->ALICEKalmanFilter(keylist_A, true, correctedOffsetTrackClusPositions, trackChi2);

  // update the seed parameters on the node tree
  // This calls circle fit and line fit (setting x, y, z, phi, eta) and sets qOverR explicitly using q from KF
  publishSeeds(seeds.first, correctedOffsetTrackClusPositions);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________________
void PrelimDistortionCorrection::publishSeeds(std::vector<TrackSeed_v2>& seeds, PositionMap& positions)
{
  int seed_index = 0;
  for(auto& seed: seeds )
  {
    /// The ALICEKF gives a better charge determination at high pT
    int q = seed.get_charge();
    TrackSeedHelper::circleFitByTaubin(&seed,positions, 7, 55);
    TrackSeedHelper::lineFit(&seed,positions, 7, 55);

    seed.set_qOverR(fabs(seed.get_qOverR()) * q);
    seed.set_phi(TrackSeedHelper::get_phi(&seed,positions));
    _track_map->insert(&seed);

    if(Verbosity() > 0)
      {
	std::cout << "Publishing seed " << seed_index
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
    if(Verbosity() > 5)
      {
	TrackSeed* readseed = _track_map->get(seed_index);
	readseed->identify();
      }

    seed_index++;
  }
}
