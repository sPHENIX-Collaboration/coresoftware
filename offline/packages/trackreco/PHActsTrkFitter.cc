/*!
 *  \file		PHActsTrkFitter.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */

#include "PHActsTrkFitter.h"
#include "MakeSourceLinks.h"

/// Tracking includes
#include <trackbase/Calibrator.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxAlignmentStateMap_v1.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <micromegas/MicromegasDefs.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <tpc/TpcDistortionCorrectionContainer.h>

#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>
#include <Acts/EventData/SourceLink.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>

#include <TDatabasePDG.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace
{
  // check vector validity
  inline bool is_valid(const Acts::Vector3 vec)
  {
    return !(std::isnan(vec.x()) || std::isnan(vec.y()) || std::isnan(vec.z()));
  }
  template <class T>
  inline T square(const T& x)
  {
    return x * x;
  }
}  // namespace

#include <trackbase/alignmentTransformationContainer.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>

PHActsTrkFitter::PHActsTrkFitter(const std::string& name)
  : SubsysReco(name)
  , m_trajectories(nullptr)
{
}

int PHActsTrkFitter::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "Setup PHActsTrkFitter" << std::endl;
  }

  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_alignStates.distortionContainers(_dcc_static, _dcc_average, _dcc_fluctuation);
  m_alignStates.actsGeometry(m_tGeometry);
  m_alignStates.clusters(m_clusterContainer);
  m_alignStates.stateMap(m_alignmentStateMap);
  m_alignStates.verbosity(Verbosity());
  m_alignStates.fieldMap(m_fieldMap);
  auto level = Acts::Logging::FATAL;
  if (Verbosity() > 5) level = Acts::Logging::VERBOSE;

  m_fitCfg.fit = ActsTrackFittingAlgorithm::makeKalmanFitterFunction(
      m_tGeometry->geometry().tGeometry,
      m_tGeometry->geometry().magField,
      true, true, 0.0, Acts::FreeToBoundCorrection(), *Acts::getDefaultLogger("Kalman", level));

  m_fitCfg.dFit = ActsTrackFittingAlgorithm::makeDirectedKalmanFitterFunction(
      m_tGeometry->geometry().tGeometry,
      m_tGeometry->geometry().magField);

  m_outlierFinder.verbosity = Verbosity();
  std::map<long unsigned int, float> chi2Cuts;
  chi2Cuts.insert(std::make_pair(10, 4));
  chi2Cuts.insert(std::make_pair(12, 4));
  chi2Cuts.insert(std::make_pair(14, 9));
  chi2Cuts.insert(std::make_pair(16, 4));
  m_outlierFinder.chi2Cuts = chi2Cuts;
  if (m_useOutlierFinder)
  {
    m_fitCfg.fit->outlierFinder(m_outlierFinder);
  }

  if (m_timeAnalysis)
  {
    m_timeFile = new TFile(std::string(Name() + ".root").c_str(),
                           "RECREATE");
    h_eventTime = new TH1F("h_eventTime", ";time [ms]",
                           100000, 0, 10000);
    h_fitTime = new TH2F("h_fitTime", ";p_{T} [GeV];time [ms]",
                         80, 0, 40, 100000, 0, 1000);
    h_updateTime = new TH1F("h_updateTime", ";time [ms]",
                            100000, 0, 1000);

    h_rotTime = new TH1F("h_rotTime", ";time [ms]",
                         100000, 0, 1000);
    h_stateTime = new TH1F("h_stateTime", ";time [ms]",
                           100000, 0, 1000);
  }

  if (m_actsEvaluator)
  {
    m_evaluator = std::make_unique<ActsEvaluator>(m_evalname);
    m_evaluator->Init(topNode);
    m_evaluator->verbosity(Verbosity());
  }

  _tpccellgeo =  findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");

  if (Verbosity() > 1)
  {
    std::cout << "Finish PHActsTrkFitter Setup" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::process_event(PHCompositeNode* topNode)
{
  PHTimer eventTimer("eventTimer");
  eventTimer.stop();
  eventTimer.restart();

  m_event++;

  auto logLevel = Acts::Logging::FATAL;

  if (m_actsEvaluator)
  {
    m_evaluator->next_event(topNode);
  }

  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
    std::cout << "Start PHActsTrkFitter::process_event" << std::endl;
    if (Verbosity() > 4)
      logLevel = Acts::Logging::VERBOSE;
  }

  /// Fill an additional track map if using the acts evaluator
  /// for proto track comparison to fitted track
  if (m_actsEvaluator)
  {
    /// wipe at the beginning of every new fit pass, so that the seeds
    /// are whatever is currently in SvtxTrackMap
    m_seedTracks->clear();
    for (const auto& [key, track] : *m_trackMap)
    {
      m_seedTracks->insert(track);
    }
  }

  loopTracks(logLevel);

  eventTimer.stop();
  auto eventTime = eventTimer.get_accumulated_time();

  if (Verbosity() > 1)
    std::cout << "PHActsTrkFitter total event time "
              << eventTime << std::endl;

  if (m_timeAnalysis)
    h_eventTime->Fill(eventTime);

  if (Verbosity() > 1)
    std::cout << "PHActsTrkFitter::process_event finished"
              << std::endl;

  // put this in the output file
  if (Verbosity() > 0)
  {
    std::cout << " SvtxTrackMap size is now " << m_trackMap->size()
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::ResetEvent(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "Reset PHActsTrkFitter" << std::endl;
  }

  m_trajectories->clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::End(PHCompositeNode* /*topNode*/)
{
  if (m_timeAnalysis)
  {
    m_timeFile->cd();
    h_fitTime->Write();
    h_eventTime->Write();
    h_rotTime->Write();
    h_stateTime->Write();
    h_updateTime->Write();
    m_timeFile->Write();
    m_timeFile->Close();
  }

  if (m_actsEvaluator)
  {
    m_evaluator->End();
  }

  if (Verbosity() > 0)
  {
    std::cout << "The Acts track fitter had " << m_nBadFits
              << " fits return an error" << std::endl;

    std::cout << "Finished PHActsTrkFitter" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrkFitter::loopTracks(Acts::Logging::Level logLevel)
{
  auto logger = Acts::getDefaultLogger("PHActsTrkFitter", logLevel);

  if (Verbosity() > 0)
  {
    std::cout << " seed map size " << m_seedMap->size() << std::endl;
  }

  for (auto trackiter = m_seedMap->begin(); trackiter != m_seedMap->end();
       ++trackiter)
  {
    TrackSeed* track = *trackiter;
    if (!track)
    {
      continue;
    }

    unsigned int tpcid = track->get_tpc_seed_index();
    unsigned int siid = track->get_silicon_seed_index();
    short int crossing_estimate = track->get_crossing_estimate();   // geometric crossing estimate

    
    /// A track seed is made for every tpc seed. Not every tpc seed
    /// has a silicon match, we skip those cases completely in pp running
    if (m_pp_mode && siid == std::numeric_limits<unsigned int>::max())
    {
      if (Verbosity() > 3) std::cout << " tpcid " << tpcid << " siid " << siid << " running in pp mode and SvtxSeedTrack has no silicon match, skip it" << std::endl;
      continue;
    }

    // get the INTT crossing number
    auto siseed = m_siliconSeeds->get(siid);
    short crossing = SHRT_MAX;
    if (siseed)
      crossing = siseed->get_crossing();
    else if (!m_pp_mode)
      crossing = 0;

    // if the crossing was not determined at all in pp running, skip this case completely
    if (m_pp_mode && crossing == SHRT_MAX && crossing_estimate == SHRT_MAX)
    {
      // Skip this in the pp case.
      if (Verbosity() > 3) std::cout << "tpcid " << tpcid << " siid " << siid << " crossing and crossing_estimate not determined, skipping track" << std::endl;
      continue;
    }

    if (Verbosity() > 1)
      {
	std::cout << "tpc and si id " << tpcid << ", " << siid << " crossing " << crossing << " crossing estimate " << crossing_estimate << std::endl;
      }

    // Can't do SC case without INTT crossing
    if( m_fitSiliconMMs && (crossing == SHRT_MAX) ) continue;

    auto tpcseed = m_tpcSeeds->get(tpcid);

    /// Need to also check that the tpc seed wasn't removed by the ghost finder
    if (!tpcseed)
    {
      std::cout << "no tpc seed" << std::endl;
      continue;
    }

    if (Verbosity() > 0)
    {
      if (siseed) std::cout << " silicon seed position is (x,y,z) = " << siseed->get_x() << "  " << siseed->get_y() << "  " << siseed->get_z() << std::endl;
      std::cout << " tpc seed position is (x,y,z) = " << tpcseed->get_x() << "  " << tpcseed->get_y() << "  " << tpcseed->get_z() << std::endl;
    }

    PHTimer trackTimer("TrackTimer");
    trackTimer.stop();
    trackTimer.restart();

    short int this_crossing = crossing;
    bool use_estimate = false;
    short int nvary = 0;
    std::vector<float> chisq_ndf;
    std::vector<SvtxTrack_v4> svtx_vec;

    if(Verbosity() > 1) { std::cout << " INTT crossing " << crossing << " crossing_estimate " << crossing_estimate << std::endl; }

    if(crossing == SHRT_MAX)
      {
	// If there is no INTT crossing, start with the crossing_estimate value, vary up and down, fit, and choose the best chisq/ndf
	use_estimate = true;
	nvary = max_bunch_search;
	if(Verbosity() > 1) { std::cout << " No INTT crossing: use crossing_estimate " << crossing_estimate << " with nvary " << nvary << std::endl; }
      }
    else
      {
	// use INTT crossing
	crossing_estimate = crossing;
      }

    // Fit this track assuming either:
    //    crossing = INTT value, if it exists (uses nvary = 0)
    //    crossing = crossing_estimate +/- max_bunch_search, if no INTT value exists

    for(short int ivary = -nvary; ivary <= nvary; ++ivary)
      {
	this_crossing = crossing_estimate + ivary;

	if(Verbosity() > 1) 
	  {
	    std::cout << "   nvary " << nvary << " trial fit with ivary " << ivary << " this_crossing = " << this_crossing << std::endl; 
	  }
 	
	ActsTrackFittingAlgorithm::MeasurementContainer measurements;

	SourceLinkVec sourceLinks;

	MakeSourceLinks makeSourceLinks;
	makeSourceLinks.initialize(_tpccellgeo);
	makeSourceLinks.setVerbosity(Verbosity());
	makeSourceLinks.set_pp_mode(m_pp_mode);

	// loop over modifiedTransformSet and replace transient elements modified for the previous track with the default transforms
	// does nothing if m_transient_id_set is empty
	makeSourceLinks.resetTransientTransformMap(
						  m_alignmentTransformationMapTransient,
						  m_transient_id_set,
						  m_tGeometry);
	if(m_use_clustermover)
	  {
	    if (siseed) 
      {
        sourceLinks = makeSourceLinks.getSourceLinksClusterMover(
								     siseed, 
								     measurements, 
								     m_clusterContainer, 
								     m_tGeometry, 
								     _dcc_static, _dcc_average, _dcc_fluctuation,
								     this_crossing);
      }
	    const auto tpcSourceLinks = makeSourceLinks.getSourceLinksClusterMover(
								       tpcseed, 
								       measurements, 
								       m_clusterContainer, 
								       m_tGeometry, 
								       _dcc_static, _dcc_average, _dcc_fluctuation,
								       this_crossing);

	    sourceLinks.insert(sourceLinks.end(), tpcSourceLinks.begin(), tpcSourceLinks.end());
	  }
	else
	  {
	    if (siseed) 
      {
        sourceLinks = makeSourceLinks.getSourceLinks(
								     siseed, 
								     measurements, 
								     m_clusterContainer, 
								     m_tGeometry, 
								     _dcc_static, _dcc_average, _dcc_fluctuation,
								     m_alignmentTransformationMapTransient, 
								     m_transient_id_set, 
								     this_crossing);
      }
	    const auto tpcSourceLinks = makeSourceLinks.getSourceLinks(
								       tpcseed, 
								       measurements, 
								       m_clusterContainer, 
								       m_tGeometry, 
								       _dcc_static, _dcc_average, _dcc_fluctuation,
								       m_alignmentTransformationMapTransient, 
								       m_transient_id_set, 
								       this_crossing);
	    sourceLinks.insert(sourceLinks.end(), tpcSourceLinks.begin(), tpcSourceLinks.end());
	  }

	// copy transient map for this track into transient geoContext
	m_transient_geocontext =  m_alignmentTransformationMapTransient;
	
	// position comes from the silicon seed, unless there is no silicon seed
	Acts::Vector3 position(0, 0, 0);
	if (siseed)
	  {
	    position(0) = siseed->get_x() * Acts::UnitConstants::cm;
	    position(1) = siseed->get_y() * Acts::UnitConstants::cm;
	    position(2) = siseed->get_z() * Acts::UnitConstants::cm;
	  }
	else
	  {
	    position(0) = tpcseed->get_x() * Acts::UnitConstants::cm;
	    position(1) = tpcseed->get_y() * Acts::UnitConstants::cm;
	    position(2) = tpcseed->get_z() * Acts::UnitConstants::cm;
	  }
	if (!is_valid(position)) continue;
	
	if (sourceLinks.empty())
	  {
	    continue;
	  }
	
	/// If using directed navigation, collect surface list to navigate
	SurfacePtrVec surfaces;
	if (m_fitSiliconMMs)
	  {
	    sourceLinks = getSurfaceVector(sourceLinks, surfaces);
	    
	    // skip if there is no surfaces
	    if (surfaces.empty()) continue;
	    
	    // make sure micromegas are in the tracks, if required
	    if (m_useMicromegas &&
		std::none_of(surfaces.begin(), surfaces.end(), [this](const auto& surface)
			     { return m_tGeometry->maps().isMicromegasSurface(surface); }))
	      {
		continue;
	      }
	  }
	
	float px = NAN;
	float py = NAN;
	float pz = NAN;
	if (m_fieldMap.find(".root") != std::string::npos)
	  {
	    px = tpcseed->get_px(m_clusterContainer, m_tGeometry);
	    py = tpcseed->get_py(m_clusterContainer, m_tGeometry);
	    pz = tpcseed->get_pz();
	  }
	else
	  {
	    float pt = fabs(1. / tpcseed->get_qOverR()) * (0.3 / 100) * std::stod(m_fieldMap);
	    float phi = tpcseed->get_phi(m_clusterContainer, m_tGeometry);
	    px = pt * std::cos(phi);
	    py = pt * std::sin(phi);
	    pz = pt * std::cosh(tpcseed->get_eta()) * std::cos(tpcseed->get_theta());
	  }

	Acts::Vector3 momentum(px, py, pz);
	if (!is_valid(momentum)) continue;

	auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
									position);

	auto actsFourPos = Acts::Vector4(position(0), position(1),
					 position(2),
					 10 * Acts::UnitConstants::ns);
	Acts::BoundSquareMatrix cov = setDefaultCovariance();

	int charge = tpcseed->get_charge();

	/// Reset the track seed with the dummy covariance
	auto seed = ActsTrackFittingAlgorithm::TrackParameters::create(
								       pSurface,
								       m_transient_geocontext,
								       actsFourPos,
								       momentum,
								       charge / momentum.norm(),
								       cov,
								       Acts::ParticleHypothesis::pion())
	                                                               .value();

	if (Verbosity() > 2)
	  {
	    printTrackSeed(seed);
	  }

	/// Set host of propagator options for Acts to do e.g. material integration
	Acts::PropagatorPlainOptions ppPlainOptions;

	auto calibptr = std::make_unique<Calibrator>();
	CalibratorAdapter calibrator{*calibptr, measurements};

	auto magcontext = m_tGeometry->geometry().magFieldContext;
	auto calibcontext = m_tGeometry->geometry().calibContext;

	ActsTrackFittingAlgorithm::GeneralFitterOptions
	  kfOptions{
	    m_transient_geocontext,
            magcontext,
            calibcontext,
            pSurface.get(),
            ppPlainOptions};

	PHTimer fitTimer("FitTimer");
	fitTimer.stop();
	fitTimer.restart();

	auto trackContainer =
	  std::make_shared<Acts::VectorTrackContainer>();
	auto trackStateContainer =
	  std::make_shared<Acts::VectorMultiTrajectory>();
	ActsTrackFittingAlgorithm::TrackContainer
	  tracks(trackContainer, trackStateContainer);

	auto result = fitTrack(sourceLinks, seed, kfOptions,
			       surfaces, calibrator, tracks);
	fitTimer.stop();
	auto fitTime = fitTimer.get_accumulated_time();

	if (Verbosity() > 1)
	  {
	    std::cout << "PHActsTrkFitter Acts fit time " << fitTime << std::endl;
	  }

	/// Check that the track fit result did not return an error
	if (result.ok()) 
	  {
	    if(use_estimate) // trial variation case
	      {
		// this is a trial variation of the crossing estimate for this track
		// Capture the chisq/ndf so we can choose the best one after all trials
		
		SvtxTrack_v4 newTrack;
		newTrack.set_tpc_seed(tpcseed);
		newTrack.set_crossing(this_crossing);
		newTrack.set_silicon_seed(siseed);
		
		if (getTrackFitResult(result, track, &newTrack, tracks, measurements))
		  {
		    float chi2ndf = newTrack.get_quality();
		    chisq_ndf.push_back(chi2ndf);
		    svtx_vec.push_back(newTrack);
		    if(Verbosity() > 1) { std::cout << "   tpcid " << tpcid << " siid " << siid << " ivary " << ivary << " this_crossing " << this_crossing << " chi2ndf " << chi2ndf << std::endl; }
		  }
		
		if(ivary != nvary)  { continue; } 
		
		// if we are here this is the last crossing iteration, evaluate the results
		if(Verbosity() > 1) { std::cout << "Finished with trial fits, chisq_ndf size is " << chisq_ndf.size() << " chisq_ndf values are:" << std::endl; }
		float best_chisq = 1000.0;
		short int best_ivary = 0;
		for(unsigned int i = 0; i<chisq_ndf.size(); ++i)
		  {
		    if(chisq_ndf[i] < best_chisq) 
		      {
			best_chisq = chisq_ndf[i];			
			best_ivary = i;
		      }
		    if(Verbosity() > 1) { std::cout << "  trial " << i  << " chisq_ndf " << chisq_ndf[i] << " best_chisq " << best_chisq << " best_ivary " << best_ivary << std::endl; }
		  }
		unsigned int trid = m_trackMap->size();
		svtx_vec[best_ivary].set_id(trid);

		m_trackMap->insertWithKey(&svtx_vec[best_ivary], trid);		    
	      } 
	    else   // case where INTT crossing is known 
	      {
		SvtxTrack_v4 newTrack;
		newTrack.set_tpc_seed(tpcseed);
		newTrack.set_crossing(this_crossing);
		newTrack.set_silicon_seed(siseed);
		
		if (m_fitSiliconMMs)
		  {
		    unsigned int trid = m_directedTrackMap->size();
		    newTrack.set_id(trid);
		    
		    if (getTrackFitResult(result, track, &newTrack, tracks, measurements))
		      {
			m_directedTrackMap->insertWithKey(&newTrack, trid);
		      }
		  }  // end insert track for SC calib fit
		else
		  {
		    unsigned int trid = m_trackMap->size();
		    newTrack.set_id(trid);
		    
		    if (getTrackFitResult(result, track, &newTrack, tracks, measurements))
		      {
			m_trackMap->insertWithKey(&newTrack, trid);
		      }
		  }  // end insert track for normal fit
	      }  // end case where INTT crossing is known
	  }
	else if (!m_fitSiliconMMs)
	  {
	    /// Track fit failed, get rid of the track from the map
	    m_nBadFits++;
	    if (Verbosity() > 1)
	      {
		std::cout << "Track fit failed for track " << m_seedMap->find(track)
			  << " with Acts error message "
			  << result.error() << ", " << result.error().message()
			  << std::endl;
	      }
	  }  // end fit failed case	
      }  // end ivary loop

    trackTimer.stop();
    auto trackTime = trackTimer.get_accumulated_time();
    
    if (Verbosity() > 1)
      {
	std::cout << "PHActsTrkFitter total single track time " << trackTime << std::endl;
      }
  }
  
  return;
}

bool PHActsTrkFitter::getTrackFitResult(FitResult& fitOutput,
                                        TrackSeed* seed, SvtxTrack* track,
                                        ActsTrackFittingAlgorithm::TrackContainer& tracks,
                                        const ActsTrackFittingAlgorithm::MeasurementContainer& measurements)
{
  /// Make a trajectory state for storage, which conforms to Acts track fit
  /// analysis tool
  std::vector<Acts::MultiTrajectoryTraits::IndexType> trackTips;
  trackTips.reserve(1);
  auto& outtrack = fitOutput.value();
  if (outtrack.hasReferenceSurface())
  {
    trackTips.emplace_back(outtrack.tipIndex());
    Trajectory::IndexedParameters indexedParams;
    indexedParams.emplace(std::pair{outtrack.tipIndex(),
                                    ActsExamples::TrackParameters{outtrack.referenceSurface().getSharedPtr(),
                                                                  outtrack.parameters(), outtrack.covariance(), outtrack.particleHypothesis()}});

    if (Verbosity() > 2)
    {
      std::cout << "Fitted parameters for track" << std::endl;
      std::cout << " position : " << outtrack.referenceSurface().localToGlobal(m_transient_geocontext, Acts::Vector2(outtrack.loc0(), outtrack.loc1()), Acts::Vector3(1, 1, 1)).transpose()

                << std::endl;
      int otcharge = outtrack.qOverP() > 0 ? 1 : -1;
      std::cout << "charge: " << otcharge << std::endl;
      std::cout << " momentum : " << outtrack.momentum().transpose()
                << std::endl;
      std::cout << "For trackTip == " << outtrack.tipIndex() << std::endl;
    }

    /// Get position, momentum from the Acts output. Update the values of
    /// the proto track
    PHTimer updateTrackTimer("UpdateTrackTimer");
    updateTrackTimer.stop();
    updateTrackTimer.restart();
    updateSvtxTrack(trackTips, indexedParams, tracks, track);

    if (m_commissioning)
    {
      if (track->get_silicon_seed() && track->get_tpc_seed())
      {
        m_alignStates.fillAlignmentStateMap(tracks, trackTips,
                                            track, measurements);
      }
    }

    updateTrackTimer.stop();
    auto updateTime = updateTrackTimer.get_accumulated_time();

    if (Verbosity() > 1)
      std::cout << "PHActsTrkFitter update SvtxTrack time "
                << updateTime << std::endl;

    if (m_timeAnalysis)
    {
      h_updateTime->Fill(updateTime);
    }

    Trajectory trajectory(tracks.trackStateContainer(),
                          trackTips, indexedParams);

    m_trajectories->insert(std::make_pair(track->get_id(), trajectory));

    if (m_actsEvaluator)
    {
      m_evaluator->evaluateTrackFit(tracks, trackTips, indexedParams, track,
                                    seed, measurements);
    }

    return true;
  }

  return false;
}

ActsTrackFittingAlgorithm::TrackFitterResult PHActsTrkFitter::fitTrack(
    const std::vector<Acts::SourceLink>& sourceLinks,
    const ActsTrackFittingAlgorithm::TrackParameters& seed,
    const ActsTrackFittingAlgorithm::GeneralFitterOptions& kfOptions,
    const SurfacePtrVec& surfSequence,
    const CalibratorAdapter& calibrator,
    ActsTrackFittingAlgorithm::TrackContainer& tracks)
{
  if (m_fitSiliconMMs)
  {
    return (*m_fitCfg.dFit)(sourceLinks, seed, kfOptions,
                            surfSequence, calibrator, tracks);
  }
  else
  {
    return (*m_fitCfg.fit)(sourceLinks, seed, kfOptions,
                           calibrator, tracks);
  }
}

SourceLinkVec PHActsTrkFitter::getSurfaceVector(const SourceLinkVec& sourceLinks,
                                                SurfacePtrVec& surfaces) const
{
  SourceLinkVec siliconMMSls;

  //   if(Verbosity() > 1)
  //     std::cout << "Sorting " << sourceLinks.size() << " SLs" << std::endl;

  for (const auto& sl : sourceLinks)
  {
    const ActsSourceLink asl = sl.get<ActsSourceLink>();
    if (Verbosity() > 1)
    {
      std::cout << "SL available on : " << asl.geometryId() << std::endl;
    }

    const auto surf = m_tGeometry->geometry().tGeometry->findSurface(asl.geometryId());
    // skip TPC surfaces
    if (m_tGeometry->maps().isTpcSurface(surf)) continue;

    // also skip micromegas surfaces if not used
    if (m_tGeometry->maps().isMicromegasSurface(surf) && !m_useMicromegas) continue;

    // update vectors
    siliconMMSls.push_back(sl);
    surfaces.push_back(surf);
  }

  /// Surfaces need to be sorted in order, i.e. from smallest to
  /// largest radius extending from target surface
  /// Add a check to ensure this
  if (!surfaces.empty())
  {
    checkSurfaceVec(surfaces);
  }

  if (Verbosity() > 1)
  {
    for (const auto& surf : surfaces)
    {
      std::cout << "Surface vector : " << surf->geometryId() << std::endl;
    }
  }

  return siliconMMSls;
}

void PHActsTrkFitter::checkSurfaceVec(SurfacePtrVec& surfaces) const
{
  for (int i = 0; i < surfaces.size() - 1; i++)
  {
    const auto& surface = surfaces.at(i);
    const auto thisVolume = surface->geometryId().volume();
    const auto thisLayer = surface->geometryId().layer();

    const auto nextSurface = surfaces.at(i + 1);
    const auto nextVolume = nextSurface->geometryId().volume();
    const auto nextLayer = nextSurface->geometryId().layer();

    /// Implement a check to ensure surfaces are sorted
    if (nextVolume == thisVolume)
    {
      if (nextLayer < thisLayer)
      {
        std::cout
            << "PHActsTrkFitter::checkSurfaceVec - "
            << "Surface not in order... removing surface"
            << surface->geometryId() << std::endl;

        surfaces.erase(surfaces.begin() + i);

        /// Subtract one so we don't skip a surface
        --i;
        continue;
      }
    }
    else
    {
      if (nextVolume < thisVolume)
      {
        std::cout
            << "PHActsTrkFitter::checkSurfaceVec - "
            << "Volume not in order... removing surface"
            << surface->geometryId() << std::endl;

        surfaces.erase(surfaces.begin() + i);

        /// Subtract one so we don't skip a surface
        --i;
        continue;
      }
    }
  }
}

void PHActsTrkFitter::updateSvtxTrack(std::vector<Acts::MultiTrajectoryTraits::IndexType>& tips,
                                      Trajectory::IndexedParameters& paramsMap,
                                      ActsTrackFittingAlgorithm::TrackContainer& tracks,
                                      SvtxTrack* track)
{
  const auto& mj = tracks.trackStateContainer();

  /// only one track tip in the track fit Trajectory
  auto& trackTip = tips.front();

  if (Verbosity() > 2)
  {
    std::cout << "Identify (proto) track before updating with acts results " << std::endl;
    track->identify();
  }

  if (!m_fitSiliconMMs)
  {
    track->clear_states();
  }

  // create a state at pathlength = 0.0
  // This state holds the track parameters, which will be updated below
  float pathlength = 0.0;
  SvtxTrackState_v1 out(pathlength);
  out.set_x(0.0);
  out.set_y(0.0);
  out.set_z(0.0);
  track->insert_state(&out);

  auto trajState =
      Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);

  const auto& params = paramsMap.find(trackTip)->second;

  /// Acts default unit is mm. So convert to cm
  track->set_x(params.position(m_transient_geocontext)(0) / Acts::UnitConstants::cm);
  track->set_y(params.position(m_transient_geocontext)(1) / Acts::UnitConstants::cm);
  track->set_z(params.position(m_transient_geocontext)(2) / Acts::UnitConstants::cm);

  track->set_px(params.momentum()(0));
  track->set_py(params.momentum()(1));
  track->set_pz(params.momentum()(2));

  track->set_charge(params.charge());
  track->set_chisq(trajState.chi2Sum);
  track->set_ndf(trajState.NDF);

  ActsTransformations rotater;
  rotater.setVerbosity(Verbosity());

  if (params.covariance())
  {
    auto rotatedCov = rotater.rotateActsCovToSvtxTrack(params);

    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        track->set_error(i, j, rotatedCov(i, j));
      }
    }
  }

  // Also need to update the state list and cluster ID list for all measurements associated with the acts track
  // loop over acts track states, copy over to SvtxTrackStates, and add to SvtxTrack
  PHTimer trackStateTimer("TrackStateTimer");
  trackStateTimer.stop();
  trackStateTimer.restart();

  if (m_fillSvtxTrackStates)
  {
    rotater.fillSvtxTrackStates(mj, trackTip, track,
				m_transient_geocontext);
  }

  trackStateTimer.stop();
  auto stateTime = trackStateTimer.get_accumulated_time();

  if (Verbosity() > 1)
    std::cout << "PHActsTrkFitter update SvtxTrackStates time "
              << stateTime << std::endl;

  if (m_timeAnalysis)
  {
    h_stateTime->Fill(stateTime);
  }

  if (Verbosity() > 2)
  {
    std::cout << " Identify fitted track after updating track states:"
              << std::endl;
    track->identify();
  }

  return;
}

Acts::BoundSquareMatrix PHActsTrkFitter::setDefaultCovariance() const
{
  Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Zero();

  /// Acts cares about the track covariance as it helps the KF
  /// know whether or not to trust the initial track seed or not.
  /// We reset it here to some loose values as it helps Acts improve
  /// the fitting.
  /// If the covariance is too loose, it won't be able to propagate,
  /// but if it is too tight, it will just "believe" the track seed over
  /// the hit data

  /// If we are using distortions, then we need to blow up the covariance
  /// a bit since the seed was created with distorted TPC clusters
  if (m_fitSiliconMMs)
  {
    cov << 1000 * Acts::UnitConstants::um, 0., 0., 0., 0., 0.,
        0., 1000 * Acts::UnitConstants::um, 0., 0., 0., 0.,
        0., 0., 0.1, 0., 0., 0.,
        0., 0., 0., 0.1, 0., 0.,
        0., 0., 0., 0., 0.005, 0.,
        0., 0., 0., 0., 0., 1.;
  }
  else
  {
    // cppcheck-suppress duplicateAssignExpression
    double sigmaD0 = 50 * Acts::UnitConstants::um;
    double sigmaZ0 = 50 * Acts::UnitConstants::um;
    // cppcheck-suppress duplicateAssignExpression
    double sigmaPhi = 1 * Acts::UnitConstants::degree;
    double sigmaTheta = 1 * Acts::UnitConstants::degree;
    double sigmaT = 1. * Acts::UnitConstants::ns;

    cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = sigmaD0 * sigmaD0;
    cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = sigmaZ0 * sigmaZ0;
    cov(Acts::eBoundTime, Acts::eBoundTime) = sigmaT * sigmaT;
    cov(Acts::eBoundPhi, Acts::eBoundPhi) = sigmaPhi * sigmaPhi;
    cov(Acts::eBoundTheta, Acts::eBoundTheta) = sigmaTheta * sigmaTheta;
    /// Acts takes this value very seriously - tuned to be in a "sweet spot"
    cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = 0.0001;
  }

  return cov;
}

void PHActsTrkFitter::printTrackSeed(const ActsTrackFittingAlgorithm::TrackParameters& seed) const
{
  std::cout
      << PHWHERE
      << " Processing proto track:"
      << std::endl;

  std::cout
       << "position: " << seed.position(m_transient_geocontext).transpose()
      << std::endl
      << "momentum: " << seed.momentum().transpose()
      << std::endl;

  std::cout << "charge : " << seed.charge() << std::endl;
  std::cout << "absolutemom : " << seed.absoluteMomentum() << std::endl;
}

int PHActsTrkFitter::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsTrkFitter::createNodes");
  }

  PHNodeIterator dstIter(topNode);
  PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  if (m_fitSiliconMMs)
  {
    m_directedTrackMap = findNode::getClass<SvtxTrackMap>(topNode,
                                                          "SvtxSiliconMMTrackMap");
    if (!m_directedTrackMap)
    {
      /// Copy this trackmap, then use it for the rest of processing
      m_directedTrackMap = new SvtxTrackMap_v2;

      PHIODataNode<PHObject>* trackNode =
          new PHIODataNode<PHObject>(m_directedTrackMap, "SvtxSiliconMMTrackMap", "PHObject");
      svtxNode->addNode(trackNode);
    }
  }

  m_trajectories = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsTrajectories");
  if (!m_trajectories)
  {
    m_trajectories = new std::map<const unsigned int, Trajectory>;
    auto node =
        new PHDataNode<std::map<const unsigned int, Trajectory>>(m_trajectories, "ActsTrajectories");
    svtxNode->addNode(node);
  }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);

  if (!m_trackMap)
  {
    m_trackMap = new SvtxTrackMap_v2;
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_trackMap, _track_map_name, "PHObject");
    svtxNode->addNode(node);
  }

  m_alignmentStateMap = findNode::getClass<SvtxAlignmentStateMap>(topNode, "SvtxAlignmentStateMap");
  if (!m_alignmentStateMap)
  {
    m_alignmentStateMap = new SvtxAlignmentStateMap_v1;
    auto node = new PHDataNode<SvtxAlignmentStateMap>(m_alignmentStateMap, "SvtxAlignmentStateMap", "PHObject");
    svtxNode->addNode(node);
  }

  if (m_actsEvaluator)
  {
    m_seedTracks = findNode::getClass<SvtxTrackMap>(topNode, _seed_track_map_name);

    if (!m_seedTracks)
    {
      m_seedTracks = new SvtxTrackMap_v2;

      PHIODataNode<PHObject>* seedNode =
          new PHIODataNode<PHObject>(m_seedTracks, _seed_track_map_name, "PHObject");
      svtxNode->addNode(seedNode);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::getNodes(PHCompositeNode* topNode)
{

  m_alignmentTransformationMap = findNode::getClass<alignmentTransformationContainer>(topNode, "alignmentTransformationContainer");
  if(!m_alignmentTransformationMap)
    {
      std::cout << PHWHERE << "alignmentTransformationContainer not on node tree. Bailing"
                << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_alignmentTransformationMapTransient = findNode::getClass<alignmentTransformationContainer>(topNode, "alignmentTransformationContainerTransient");
  if(!m_alignmentTransformationMapTransient)
    {
      std::cout << PHWHERE << "alignmentTransformationContainerTransient not on node tree. Bailing"
                << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }


  m_tpcSeeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!m_tpcSeeds)
  {
    std::cout << PHWHERE << "TpcTrackSeedContainer not on node tree. Bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_siliconSeeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!m_siliconSeeds)
  {
    std::cout << PHWHERE << "SiliconTrackSeedContainer not on node tree. Bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterContainer)
  {
    std::cout << PHWHERE
              << "No trkr cluster container, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << "ActsGeometry not on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_seedMap = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");
  if (!m_seedMap)
  {
    std::cout << "No Svtx seed map on node tree. Exiting."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc distortion corrections
  _dcc_static = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerStatic");
  if (_dcc_static)
  {
    std::cout << PHWHERE << "  found static TPC distortion correction container" << std::endl;
  }
  _dcc_average = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerAverage");
  if (_dcc_average)
  {
    std::cout << PHWHERE << "  found average TPC distortion correction container" << std::endl;
  }
  _dcc_fluctuation = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerFluctuation");
  if (_dcc_fluctuation)
  {
    std::cout << PHWHERE << "  found fluctuation TPC distortion correction container" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
