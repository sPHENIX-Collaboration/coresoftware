/*!
 *  \file		PHActsTrkFitter.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */


#include "PHActsTrkFitter.h"

#include "ActsPropagator.h"
#include "MakeSourceLinks.h"

#include <tpc/TpcDistortionCorrectionContainer.h>

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
//#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxTrackState_v3.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <g4detectors/PHG4TpcGeomContainer.h>

#include <micromegas/MicromegasDefs.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>
#include <Acts/EventData/SourceLink.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

namespace
{
  // check vector validity
  inline bool is_valid(const Acts::Vector3& vec)
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

  // configure alignStates
  m_alignStates.loadNodes(topNode);
  m_alignStates.verbosity(Verbosity());
  m_alignStates.fieldMap(m_fieldMap);

  // detect const field
  std::istringstream stringline(m_fieldMap);
  stringline >> fieldstrength;
  if (!stringline.fail())  // it is a float
  {
    m_ConstField = true;
  }

  auto level = Acts::Logging::FATAL;
  if (Verbosity() > 5)
  {
    level = Acts::Logging::VERBOSE;
  }

  m_fitCfg.fit = ActsTrackFittingAlgorithm::makeKalmanFitterFunction(
      m_tGeometry->geometry().tGeometry,
      m_tGeometry->geometry().magField,
      true, true, 0.0, Acts::FreeToBoundCorrection(), *Acts::getDefaultLogger("Kalman", level));

  m_fitCfg.dFit = ActsTrackFittingAlgorithm::makeDirectedKalmanFitterFunction(
      m_tGeometry->geometry().tGeometry,
      m_tGeometry->geometry().magField);

  MaterialSurfaceSelector selector;
  if (m_fitSiliconMMs || m_directNavigation)
  {
    m_tGeometry->geometry().tGeometry->visitSurfaces(selector,false);
    //std::cout<<"selector.surfaces.size() "<<selector.surfaces.size()<<std::endl;
    m_materialSurfaces = selector.surfaces;
  }

  m_outlierFinder.verbosity = Verbosity();
  std::map<long unsigned int, float> chi2Cuts;
  chi2Cuts.insert(std::make_pair(10, 4));
  chi2Cuts.insert(std::make_pair(12, 4));
  chi2Cuts.insert(std::make_pair(14, 9));
  chi2Cuts.insert(std::make_pair(16, 4));
  m_outlierFinder.chi2Cuts = chi2Cuts;
  if (m_useOutlierFinder)
  {
    m_outlierFinder.m_tGeometry = m_tGeometry;
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
    if(m_actsEvaluator && !m_simActsEvaluator)
    {
      m_evaluator->isData();
    }
    m_evaluator->verbosity(Verbosity());
  }

  _tpccellgeo = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  if (!_tpccellgeo)
    {
      std::cout << PHWHERE << " unable to find DST node TPCGEOMCONTAINER" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

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
  m_nBadFits = 0;
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
    {
      logLevel = Acts::Logging::VERBOSE;
    }
  }

  // in case the track map already exist in the file, we want to replace it
  m_trackMap->Reset();

  loopTracks(logLevel);

  eventTimer.stop();
  auto eventTime = eventTimer.get_accumulated_time();

  if (Verbosity() > 1)
  {
    std::cout << "PHActsTrkFitter total event time "
              << eventTime << std::endl;
  }

  if (m_timeAnalysis)
  {
    h_eventTime->Fill(eventTime);
  }

  if (Verbosity() > 1)
  {
    std::cout << "PHActsTrkFitter::process_event finished"
              << std::endl;
  }

  // put this in the output file
  if (Verbosity() > 0)
  {
    std::cout << "The Acts track fitter had " << m_nBadFits
              << " fits return an error" << std::endl;
    std::cout << " seed map size " << m_seedMap->size() << std::endl;

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
  if(m_useOutlierFinder)
  {
    m_outlierFinder.Write();
  }
  if (Verbosity() > 0)
  {
    std::cout << "Finished PHActsTrkFitter" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrkFitter::loopTracks(Acts::Logging::Level logLevel)
{
  auto logger = Acts::getDefaultLogger("PHActsTrkFitter", logLevel);

  for (auto track : *m_seedMap)
  {
    if (!track)
    {
      continue;
    }

    unsigned int tpcid = track->get_tpc_seed_index();
    unsigned int siid = track->get_silicon_seed_index();

    // capture the input crossing value, and set crossing parameters
    //==============================
    short silicon_crossing =  SHRT_MAX;
    auto siseed = m_siliconSeeds->get(siid);
    if(siseed)
      {
	silicon_crossing = siseed->get_crossing();
      }
    short crossing = silicon_crossing;
    short int crossing_estimate = crossing;

    if(m_enable_crossing_estimate)
      {
	crossing_estimate = track->get_crossing_estimate();  // geometric crossing estimate from matcher
     }
    //===============================


    // must have silicon seed with valid crossing if we are doing a SC calibration fit
    if (m_fitSiliconMMs)
      {
	if( (siid == std::numeric_limits<unsigned int>::max()) || (silicon_crossing == SHRT_MAX))
	  {
	    continue;
	  }
      }

    // do not skip TPC only tracks, just set crossing to the nominal zero
    if(!siseed)
      {
	crossing = 0;
      }

    if (Verbosity() > 1)
    {
      if(siseed)
	{
	  std::cout << "tpc and si id " << tpcid << ", " << siid << " silicon_crossing " << silicon_crossing
		    << " crossing " << crossing << " crossing estimate " << crossing_estimate << std::endl;
	}
    }

    auto tpcseed = m_tpcSeeds->get(tpcid);

    /// Need to also check that the tpc seed wasn't removed by the ghost finder
    if (!tpcseed)
    {
      std::cout << "no tpc seed" << std::endl;
      continue;
    }

    if (Verbosity() > 0)
    {
      if (siseed)
      {
        const auto si_position = TrackSeedHelper::get_xyz(siseed);
        const auto tpc_position = TrackSeedHelper::get_xyz(tpcseed);
        std::cout << "    silicon seed position is (x,y,z) = " << si_position.x() << "  " << si_position.y() << "  " << si_position.z() << std::endl;
        std::cout << "    tpc seed position is (x,y,z) = " << tpc_position.x() << "  " << tpc_position.y() << "  " << tpc_position.z() << std::endl;
      }
    }

    PHTimer trackTimer("TrackTimer");
    trackTimer.stop();
    trackTimer.restart();

    if (Verbosity() > 1 && siseed)
    {
      std::cout << " m_pp_mode " << m_pp_mode << " m_enable_crossing_estimate " << m_enable_crossing_estimate
        << " INTT crossing " << crossing << " crossing_estimate " << crossing_estimate << std::endl;
    }

    short int this_crossing = crossing;
    bool use_estimate = false;
    short int nvary = 0;
    std::vector<float> chisq_ndf;
    std::vector<SvtxTrack_v4> svtx_vec;

    if(m_pp_mode)
      {
	if (m_enable_crossing_estimate && crossing == SHRT_MAX)
	  {
	    // this only happens if there is a silicon seed but no assigned INTT crossing, and only in pp_mode
	    // If there is no INTT crossing, start with the crossing_estimate value, vary up and down, fit, and choose the best chisq/ndf
	    use_estimate = true;
	    nvary = max_bunch_search;
	    if (Verbosity() > 1)
	      {
		std::cout << " No INTT crossing: use crossing_estimate " << crossing_estimate << " with nvary " << nvary << std::endl;
	      }
	  }
	else
	  {
	    // use INTT crossing
	    crossing_estimate = crossing;
	  }
      }
    else
      {
	// non pp mode, we want only crossing zero, veto others
	if(siseed && silicon_crossing != 0)
	  {
	    crossing = 0;
	    //continue;
	  }
	crossing_estimate = crossing;
      }

    // Fit this track assuming either:
    //    crossing = INTT value, if it exists (uses nvary = 0)
    //    crossing = crossing_estimate +/- max_bunch_search, if no INTT value exists and m_enable_crossing_estimate flag is set.

    for (short int ivary = -nvary; ivary <= nvary; ++ivary)
    {
      this_crossing = crossing_estimate + ivary;

      if (Verbosity() > 1)
      {
        std::cout << "   nvary " << nvary << " trial fit with ivary " << ivary << " this_crossing = " << this_crossing << std::endl;
      }

      ActsTrackFittingAlgorithm::MeasurementContainer measurements;

      SourceLinkVec sourceLinks;

      MakeSourceLinks makeSourceLinks;
      makeSourceLinks.initialize(_tpccellgeo);
      makeSourceLinks.setVerbosity(Verbosity());
      makeSourceLinks.set_pp_mode(m_pp_mode);
      for(const auto& layer : m_ignoreLayer)
      {
        makeSourceLinks.ignoreLayer(layer);
      }
      // loop over modifiedTransformSet and replace transient elements modified for the previous track with the default transforms
      // does nothing if m_transient_id_set is empty
      makeSourceLinks.resetTransientTransformMap(
        m_alignmentTransformationMapTransient,
        m_transient_id_set,
        m_tGeometry);

      if (m_use_clustermover)
      {
        // make source links using cluster mover after making distortion correction
        if (siseed && !m_ignoreSilicon)
        {
          // silicon source links
          sourceLinks = makeSourceLinks.getSourceLinksClusterMover(
            siseed,
            measurements,
            m_clusterContainer,
            m_tGeometry,
            m_globalPositionWrapper,
            this_crossing);
        }

        // tpc source links
        const auto tpcSourceLinks = makeSourceLinks.getSourceLinksClusterMover(
          tpcseed,
          measurements,
          m_clusterContainer,
          m_tGeometry,
          m_globalPositionWrapper,
          this_crossing);

        // add tpc sourcelinks to silicon source links
        sourceLinks.insert(sourceLinks.end(), tpcSourceLinks.begin(), tpcSourceLinks.end());

      } else {

        // make source links using transient transforms for distortion corrections
        if(Verbosity() > 1)
        { std::cout << "Calling getSourceLinks for si seed, siid " << siid << " and tpcid " << tpcid << std::endl; }

        if (siseed && !m_ignoreSilicon)
        {
          // silicon source links
          sourceLinks = makeSourceLinks.getSourceLinks(
            siseed,
            measurements,
            m_clusterContainer,
            m_tGeometry,
            m_globalPositionWrapper,
            m_alignmentTransformationMapTransient,
            m_transient_id_set,
            this_crossing);
        }

        if(Verbosity() > 1)
        { std::cout << "Calling getSourceLinks for tpc seed, siid " << siid << " and tpcid " << tpcid << std::endl; }

        // tpc source links
        const auto tpcSourceLinks = makeSourceLinks.getSourceLinks(
          tpcseed,
          measurements,
          m_clusterContainer,
          m_tGeometry,
          m_globalPositionWrapper,
          m_alignmentTransformationMapTransient,
          m_transient_id_set,
          this_crossing);

        // add tpc sourcelinks to silicon source links
        sourceLinks.insert(sourceLinks.end(), tpcSourceLinks.begin(), tpcSourceLinks.end());
      }

      // copy transient map for this track into transient geoContext
      m_transient_geocontext = m_alignmentTransformationMapTransient;

      // position comes from the silicon seed, unless there is no silicon seed
      Acts::Vector3 position(0, 0, 0);
      if (siseed)
      {
        position = TrackSeedHelper::get_xyz(siseed)*Acts::UnitConstants::cm;
      }
      if(!siseed || !is_valid(position) || m_ignoreSilicon)
      {
        position = TrackSeedHelper::get_xyz(tpcseed)*Acts::UnitConstants::cm;
      }
      if (!is_valid(position))
      {
       if(Verbosity() > 4)
        {
          std::cout << "Invalid position of " << position.transpose() << std::endl;
        }
        continue;
      }

      if (sourceLinks.empty())
      {
        continue;
      }

      /// If using directed navigation, collect surface list to navigate
      SurfacePtrVec surfaces_tmp;
      SurfacePtrVec surfaces;
      if (m_fitSiliconMMs || m_directNavigation)
      {
        sourceLinks = getSurfaceVector(sourceLinks, surfaces_tmp);

        // skip if there is no surfaces
        if (surfaces_tmp.empty())
        {
          continue;
        }

        for (const auto& surface_apr : m_materialSurfaces)
        {
          if(m_forceSiOnlyFit)
          {
            if(surface_apr->geometryId().volume() >12)
            {
              continue;
            }
          }
          bool pop_flag = false;
          if(surface_apr->geometryId().approach() == 1)
          {
            surfaces.push_back(surface_apr);
          }
          else
          {
            pop_flag = true;
            for (const auto& surface_sns: surfaces_tmp)
            {
              if (surface_apr->geometryId().volume() == surface_sns->geometryId().volume())
              {
                if ( surface_apr->geometryId().layer()==surface_sns->geometryId().layer())
                {
                  pop_flag = false;
                  surfaces.push_back(surface_sns);
                }
              }
            }
            if (!pop_flag)
            {
              surfaces.push_back(surface_apr);
            }
            else
            {
              surfaces.pop_back();
              pop_flag = false;
            }
            if (surface_apr->geometryId().volume() == 12&& surface_apr->geometryId().layer()==8)
            {
              for (const auto& surface_sns: surfaces_tmp)
              {
                if (14 == surface_sns->geometryId().volume())
                {
                  surfaces.push_back(surface_sns);
                }
              }
            }
          }
        }
        checkSurfaceVec(surfaces);
        if (Verbosity() > 1)
        {
          for (const auto& surf : surfaces)
          {
            std::cout << "Surface vector : " << surf->geometryId() << std::endl;
          }
        }

        if (m_fitSiliconMMs)
        {
          // make sure micromegas are in the tracks, if required
          if (m_useMicromegas &&
            std::none_of(surfaces.begin(), surfaces.end(), [this](const auto& surface)
          { return m_tGeometry->maps().isMicromegasSurface(surface); }))
        {
          continue;
        }
      }
    }

      float px = std::numeric_limits<float>::quiet_NaN();
      float py = std::numeric_limits<float>::quiet_NaN();
      float pz = std::numeric_limits<float>::quiet_NaN();

      // get phi and theta from the silicon seed, momentum from the TPC seed
      float seedphi = 0;
      float seedtheta = 0;
      float seedeta = 0;
      if(siseed)
      {
        seedphi = siseed->get_phi();
        seedtheta = siseed->get_theta();
        seedeta = siseed->get_eta();
      }
      else
      {
        seedphi = tpcseed->get_phi();
        seedtheta = tpcseed->get_theta();
        seedeta = tpcseed->get_eta();
      }

      float seedpt = tpcseed->get_pt();

      if (m_ConstField)
      {
        float pt = fabs(1. / tpcseed->get_qOverR()) * (0.3 / 100) * fieldstrength;
        float phi = seedphi;
        float eta = seedeta;
        float theta = seedtheta;
        px = pt * std::cos(phi);
        py = pt * std::sin(phi);
        pz = pt * std::cosh(eta) * std::cos(theta);
      } else {
        px = seedpt * std::cos(seedphi);
        py = seedpt * std::sin(seedphi);
        pz = seedpt * std::cosh(seedeta) * std::cos(seedtheta);
      }

      Acts::Vector3 momentum(px, py, pz);
      if (!is_valid(momentum))
      {
        if(Verbosity() > 4)
        {
          std::cout << "Invalid momentum of " << momentum.transpose() << std::endl;
        }
        continue;
      }

      auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>( position);

      Acts::Vector4 actsFourPos(position(0), position(1), position(2), 10 * Acts::UnitConstants::ns);
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

      auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
      auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
      ActsTrackFittingAlgorithm::TrackContainer tracks(trackContainer, trackStateContainer);

      if(Verbosity() > 1)
      {  std::cout << "Calling fitTrack for track with siid " << siid << " tpcid " << tpcid << " crossing " << crossing << std::endl; }

      auto result = fitTrack(sourceLinks, seed, kfOptions, surfaces, calibrator, tracks);
      fitTimer.stop();

      if (Verbosity() > 1)
      {
        const auto fitTime = fitTimer.get_accumulated_time();
        std::cout << "PHActsTrkFitter Acts fit time " << fitTime << std::endl;
      }

      /// Check that the track fit result did not return an error
      if (result.ok())
      {
        if (use_estimate)  // trial variation case
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
            if (Verbosity() > 1)
            {
              std::cout << "   tpcid " << tpcid << " siid " << siid << " ivary " << ivary << " this_crossing " << this_crossing << " chi2ndf " << chi2ndf << std::endl;
            }
          }

          if (ivary != nvary)
          {
            if(Verbosity() > 3)
            {
              std::cout << "Skipping track fit for trial variation" << std::endl;
            }
            continue;
          }

          // if we are here this is the last crossing iteration, evaluate the results
          if (Verbosity() > 1)
          {
            std::cout << "Finished with trial fits, chisq_ndf size is " << chisq_ndf.size() << " chisq_ndf values are:" << std::endl;
          }
          float best_chisq = 1000.0;
          short int best_ivary = 0;
          for (unsigned int i = 0; i < chisq_ndf.size(); ++i)
          {
            if (chisq_ndf[i] < best_chisq)
            {
              best_chisq = chisq_ndf[i];
              best_ivary = i;
            }
            if (Verbosity() > 1)
            {
              std::cout << "  trial " << i << " chisq_ndf " << chisq_ndf[i] << " best_chisq " << best_chisq << " best_ivary " << best_ivary << std::endl;
            }
          }
          unsigned int trid = m_trackMap->size();
          svtx_vec[best_ivary].set_id(trid);

          m_trackMap->insertWithKey(&svtx_vec[best_ivary], trid);
        }
        else  // case where INTT crossing is known
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

              // insert in dedicated map
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
        }    // end case where INTT crossing is known




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
    }    // end ivary loop

    trackTimer.stop();
    auto trackTime = trackTimer.get_accumulated_time();

    if (Verbosity() > 1)
    {
      std::cout << "PHActsTrkFitter total single track time " << trackTime << std::endl;
    }
  }

  return;
}

bool PHActsTrkFitter::getTrackFitResult(
  const FitResult& fitOutput,
  TrackSeed* seed, SvtxTrack* track,
  const ActsTrackFittingAlgorithm::TrackContainer& tracks,
  const ActsTrackFittingAlgorithm::MeasurementContainer& measurements)
{
  /// Make a trajectory state for storage, which conforms to Acts track fit
  /// analysis tool
  std::vector<Acts::MultiTrajectoryTraits::IndexType> trackTips;
  trackTips.reserve(1);
  const auto& outtrack = fitOutput.value();
  if (outtrack.hasReferenceSurface())
  {
    trackTips.emplace_back(outtrack.tipIndex());
    Trajectory::IndexedParameters indexedParams;

    // retrieve track parameters from fit result
    Acts::BoundTrackParameters parameters = ActsExamples::TrackParameters(outtrack.referenceSurface().getSharedPtr(),
      outtrack.parameters(), outtrack.covariance(), outtrack.particleHypothesis());

    indexedParams.emplace(
      outtrack.tipIndex(),
      ActsExamples::TrackParameters{outtrack.referenceSurface().getSharedPtr(),
      outtrack.parameters(), outtrack.covariance(), outtrack.particleHypothesis()});

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
    {
      std::cout << "PHActsTrkFitter update SvtxTrack time "
                << updateTime << std::endl;
    }

    if (m_timeAnalysis)
    {
      h_updateTime->Fill(updateTime);
    }

    Trajectory trajectory(tracks.trackStateContainer(),
                          trackTips, indexedParams);

    if (m_actsEvaluator)
    {
      m_evaluator->evaluateTrackFit(tracks, trackTips, indexedParams, track,
                                    seed, measurements);
    }

    return true;
  }

  return false;
}

//__________________________________________________________________________________
ActsTrackFittingAlgorithm::TrackFitterResult PHActsTrkFitter::fitTrack(
    const std::vector<Acts::SourceLink>& sourceLinks,
    const ActsTrackFittingAlgorithm::TrackParameters& seed,
    const ActsTrackFittingAlgorithm::GeneralFitterOptions& kfOptions,
    const SurfacePtrVec& surfSequence,
    const CalibratorAdapter& calibrator,
    ActsTrackFittingAlgorithm::TrackContainer& tracks)
{
  if (m_fitSiliconMMs || m_directNavigation)
  {
    return (*m_fitCfg.dFit)(sourceLinks, seed, kfOptions, surfSequence, calibrator, tracks);
  } else {
    return (*m_fitCfg.fit)(sourceLinks, seed, kfOptions, calibrator, tracks);
  }
}

//__________________________________________________________________________________
SourceLinkVec PHActsTrkFitter::getSurfaceVector(const SourceLinkVec& sourceLinks, SurfacePtrVec& surfaces) const
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
    if (m_fitSiliconMMs)
    {
      // skip TPC surfaces
      if (m_tGeometry->maps().isTpcSurface(surf))
      {
        continue;
      }

      // also skip micromegas surfaces if not used
      if (m_tGeometry->maps().isMicromegasSurface(surf) && !m_useMicromegas)
      {
        continue;
      }
    }

    if(m_forceSiOnlyFit)
    {
      if(m_tGeometry->maps().isMicromegasSurface(surf)||m_tGeometry->maps().isTpcSurface(surf))
      {
        continue;
      }
    }

    // update vectors
    siliconMMSls.push_back(sl);
    surfaces.push_back(surf);
  }

  if (Verbosity() > 10)
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
  for (unsigned int i = 0; i < surfaces.size() - 1; i++)
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

void PHActsTrkFitter::updateSvtxTrack(
  const std::vector<Acts::MultiTrajectoryTraits::IndexType>& tips,
  const Trajectory::IndexedParameters& paramsMap,
  const ActsTrackFittingAlgorithm::TrackContainer& tracks,
  SvtxTrack* track)
{
  const auto& mj = tracks.trackStateContainer();

  /// only one track tip in the track fit Trajectory
  const auto& trackTip = tips.front();

  if (Verbosity() > 2)
  {
    std::cout << "Identify (proto) track before updating with acts results " << std::endl;
    track->identify();
  }

  if (!m_fitSiliconMMs && !m_forceSiOnlyFit)
  {
    track->clear_states();
  }

  // create a state at pathlength = 0.0
  // This state holds the track parameters, which will be updated below
  float pathlength = 0.0;
  //  SvtxTrackState_v1 out(pathlength);
  SvtxTrackState_v3 out(pathlength);
  out.set_localX(0.0);
  out.set_localY(0.0);
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

  ActsTransformations transformer;
  transformer.setVerbosity(Verbosity());

  if (params.covariance())
  {
    auto rotatedCov = transformer.rotateActsCovToSvtxTrack(params);

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
  { transformer.fillSvtxTrackStates(mj, trackTip, track, m_transient_geocontext); }

  // in using silicon mm fit also extrapolate track parameters to all TPC surfaces with clusters
  // get all tpc clusters
  auto seed = track->get_tpc_seed();
  if( m_fitSiliconMMs && seed )
  {

    // acts propagator
    ActsPropagator propagator(m_tGeometry);

    // loop over cluster keys associated to TPC seed
    for( auto key_iter = seed->begin_cluster_keys(); key_iter != seed->end_cluster_keys(); ++key_iter )
    {
      const auto& cluskey = *key_iter;

      // make sure cluster is from TPC
      const auto detId = TrkrDefs::getTrkrId(cluskey);
      if (detId != TrkrDefs::tpcId)
      { continue; }

      // get layer, propagate
      const auto layer = TrkrDefs::getLayer(cluskey);
      auto result = propagator.propagateTrack(params, layer);
      if( !result.ok() ) continue;

      // get path length and extrapolated parameters
      auto& [pathLength, trackStateParams] = result.value();
      pathLength /= Acts::UnitConstants::cm;

      // create track state and add to track
      transformer.addTrackState(track, cluskey, pathLength, trackStateParams, m_transient_geocontext);
    }
  }

  trackStateTimer.stop();
  auto stateTime = trackStateTimer.get_accumulated_time();

  if (Verbosity() > 1)
  {
    std::cout << "PHActsTrkFitter update SvtxTrackStates time "
              << stateTime << std::endl;
  }

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
    double sigmaPhi = 1. * Acts::UnitConstants::degree;
    double sigmaTheta = 1. * Acts::UnitConstants::degree;
    double sigmaT = 1. * Acts::UnitConstants::ns;

    cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = sigmaD0 * sigmaD0;
    cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = sigmaZ0 * sigmaZ0;
    cov(Acts::eBoundTime, Acts::eBoundTime) = sigmaT * sigmaT;
    cov(Acts::eBoundPhi, Acts::eBoundPhi) = sigmaPhi * sigmaPhi;
    cov(Acts::eBoundTheta, Acts::eBoundTheta) = sigmaTheta * sigmaTheta;
    /// Acts takes this value very seriously - tuned to be in a "sweet spot"
    //    cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = 0.0001;
    cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = 0.025;
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

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);

  if (!m_trackMap)
  {
    m_trackMap = new SvtxTrackMap_v2;
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_trackMap, _track_map_name, "PHObject");
    svtxNode->addNode(node);
  }

  m_alignmentStateMap = findNode::getClass<SvtxAlignmentStateMap>(topNode, _svtx_alignment_state_map_name);
  if (!m_alignmentStateMap)
  {
    m_alignmentStateMap = new SvtxAlignmentStateMap_v1;
    auto node = new PHDataNode<SvtxAlignmentStateMap>(m_alignmentStateMap, _svtx_alignment_state_map_name, "PHObject");
    svtxNode->addNode(node);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::getNodes(PHCompositeNode* topNode)
{
  m_alignmentTransformationMap = findNode::getClass<alignmentTransformationContainer>(topNode, "alignmentTransformationContainer");
  if (!m_alignmentTransformationMap)
  {
    std::cout << PHWHERE << "alignmentTransformationContainer not on node tree. Bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_alignmentTransformationMapTransient = findNode::getClass<alignmentTransformationContainer>(topNode, "alignmentTransformationContainerTransient");
  if (!m_alignmentTransformationMapTransient)
  {
    std::cout << PHWHERE << "alignmentTransformationContainerTransient not on node tree. Bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc seeds
  m_tpcSeeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!m_tpcSeeds)
  {
    std::cout << PHWHERE << "TpcTrackSeedContainer not on node tree. Bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // silicon seeds
  m_siliconSeeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!m_siliconSeeds)
  {
    std::cout << PHWHERE << "SiliconTrackSeedContainer not on node tree. Bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // clusters
  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  if (!m_clusterContainer)
  {
    std::cout << PHWHERE
              << "No trkr cluster container, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << "ActsGeometry not on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // track seeds
  m_seedMap = findNode::getClass<TrackSeedContainer>(topNode, _svtx_seed_map_name);
  if (!m_seedMap)
  {
    std::cout << "No Svtx seed map on node tree. Exiting."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc global position wrapper
  m_globalPositionWrapper.loadNodes(topNode);
  m_globalPositionWrapper.set_suppressCrossing(true);

  return Fun4AllReturnCodes::EVENT_OK;
}
