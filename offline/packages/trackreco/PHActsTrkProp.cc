
#include "PHActsTrkProp.h"
#include "PHActsTracks.h"

/// Fun4All includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

/// Tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <Acts/EventData/ChargePolicy.hpp>
#include <Acts/EventData/SingleCurvilinearTrackParameters.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/EventData/NeutralParameters.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Propagator/AbortList.hpp>
#include <Acts/Propagator/ActionList.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/Units.hpp>

#include <ACTFW/Framework/BareAlgorithm.hpp>
#include <ACTFW/Framework/ProcessCode.hpp>
#include <ACTFW/Framework/RandomNumbers.hpp>
#include <ACTFW/Framework/WhiteBoard.hpp>
#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/Framework/AlgorithmContext.hpp>

#include <cmath>
#include <iostream>
#include <vector>
#include <utility>
#include <TMatrixDSym.h>

PHActsTrkProp::PHActsTrkProp(const std::string& name)
  : PHTrackPropagating(name)
  , m_event(0)
  , m_actsGeometry(nullptr)
  , m_minTrackPt(0.15)

{
  Verbosity(0);
}

 PHActsTrkProp::~PHActsTrkProp()
{
}

int PHActsTrkProp::Setup(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkProp::Process()
{
  m_event++;

  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
    std::cout << "Start PHActsTrkfitter::process_event" << std::endl;
  }

  PerigeeSurface surface = Acts::Surface::makeShared<Acts::PerigeeSurface>
    (Acts::Vector3D(0., 0., 0.));

  
  for(std::vector<ActsTrack>::iterator trackIter = m_actsTracks->begin();
      trackIter != m_actsTracks->end();
      ++trackIter)
	
    {
      ActsTrack track = *trackIter;
     
      const FW::TrackParameters actsTrack = track.trackParams;
      
      PropagationOutput pOutput = propagate(actsTrack);

    }


  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkProp::End(PHCompositeNode* topNode)
{
  if (Verbosity() > 10)
  {
    std::cout << "Finished PHActsTrkProp" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

 
PropagationOutput PHActsTrkProp::propagate(FW::TrackParameters parameters)
{
  PropagationOutput pOutput;
  
  // The step length logger for testing & end of world aborter
  using MaterialInteractor = Acts::MaterialInteractor;
  using SteppingLogger     = Acts::detail::SteppingLogger;
  using DebugOutput        = Acts::detail::DebugOutputActor;
  using EndOfWorld         = Acts::detail::EndOfWorldReached;
  
  // Action list and abort list
  using ActionList
    = Acts::ActionList<SteppingLogger, MaterialInteractor, DebugOutput>;
  using AbortList         = Acts::AbortList<EndOfWorld>;
  using PropagatorOptions = Acts::PropagatorOptions<ActionList, AbortList>;
  
  PropagatorOptions options(m_actsGeometry->geoContext, 
			    m_actsGeometry->magFieldContext);
  options.pathLimit = std::numeric_limits<double>::max();
  options.debug     = true;
  
  // Activate loop protection at some pt value
  options.loopProtection
    = (Acts::VectorHelpers::perp(parameters.momentum())
       < m_minTrackPt);

  // Switch the material interaction on/off & eventually into logging mode
  // Should all of these switches be configurable from e.g. constructor?
  auto& mInteractor = options.actionList.get<MaterialInteractor>();
  mInteractor.multipleScattering = true;
  mInteractor.energyLoss         = true;
  mInteractor.recordInteractions = true;

  // Set a maximum step size
  options.maxStepSize = 3. * Acts::UnitConstants::mm;
  
  /*
  // Propagate using the propagator
  /// Will need to setup some acts class that creates the propagator
  /// similar to the track fitter
  const auto& result
    = m_cfg.propagator.propagate(startParameters, options).value();

  auto steppingResults = result.template get<SteppingLogger::result_type>();
  
  // Set the stepping result
  pOutput.first = std::move(steppingResults.steps);
  // Also set the material recording result - if configured
  if (mInteractor.recordInteractions) {
    auto materialResult
      = result.template get<MaterialInteractor::result_type>();
    pOutput.second = std::move(materialResult);
  }
  
  if(Verbosity() > 1)
    {
      auto& debugResult = result.template get<DebugOutput::result_type>();
      std::cout << "Acts Propagator Debug: " << debugResult.debugString 
		<< std::endl;
    }
  */
  return pOutput;
 
}



Acts::BoundSymMatrix PHActsTrkProp::getActsCovMatrix(const SvtxTrack *track)
{
  Acts::BoundSymMatrix matrix = Acts::BoundSymMatrix::Zero();
  const double px = track->get_px();
  const double py = track->get_py();
  const double pz = track->get_pz();
  const double p = sqrt(px * px + py * py + pz * pz);

  // Get the track seed covariance matrix
  // These are the variances, so the std devs are sqrt(seed_cov[i][j])
  TMatrixDSym seed_cov(6);
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      seed_cov[i][j] = track->get_error(i, j);
    }
  }

  const double sigmap = sqrt(px * px * seed_cov[3][3] + py * py * seed_cov[4][4] + pz * pz * seed_cov[5][5]) / p;

  // Need to convert seed_cov from x,y,z,px,py,pz basis to Acts basis of
  // x,y,phi/theta of p, qoverp, time
  double phi = track->get_phi();

  const double pxfracerr = seed_cov[3][3] / (px * px);
  const double pyfracerr = seed_cov[4][4] / (py * py);
  const double phiPrefactor = fabs(py) / (fabs(px) * (1 + (py / px) * (py / px)));
  const double sigmaPhi = phi * phiPrefactor * sqrt(pxfracerr + pyfracerr);
  const double theta = acos(pz / p);
  const double thetaPrefactor = ((fabs(pz)) / (p * sqrt(1 - (pz / p) * (pz / p))));
  const double sigmaTheta = thetaPrefactor * sqrt(sigmap * sigmap / (p * p) + seed_cov[5][5] / (pz * pz));
  const double sigmaQOverP = sigmap / (p * p);

  // Just set to 0 for now?
  const double sigmaTime = 0;

  if (Verbosity() > 10)
  {
    std::cout << "Track (px,py,pz,p) = (" << px << "," << py
              << "," << pz << "," << p << ")" << std::endl;
    std::cout << "Track covariance matrix: " << std::endl;

    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        std::cout << seed_cov[i][j] << ", ";
      }
      std::cout << std::endl;
    }
    std::cout << "Corresponding uncertainty calculations: " << std::endl;
    std::cout << "perr: " << sigmap << std::endl;
    std::cout << "phi: " << phi << std::endl;
    std::cout << "pxfracerr: " << pxfracerr << std::endl;
    std::cout << "pyfracerr: " << pyfracerr << std::endl;
    std::cout << "phiPrefactor: " << phiPrefactor << std::endl;
    std::cout << "sigmaPhi: " << sigmaPhi << std::endl;
    std::cout << "theta: " << theta << std::endl;
    std::cout << "thetaPrefactor: " << thetaPrefactor << std::endl;
    std::cout << "sigmaTheta: " << sigmaTheta << std::endl;
    std::cout << "sigmaQOverP: " << sigmaQOverP << std::endl;
  }

  /// Seed covariances are already variances, so don't need to square them
  matrix(Acts::eLOC_0, Acts::eLOC_0) = seed_cov[0][0];
  matrix(Acts::eLOC_1, Acts::eLOC_1) = seed_cov[1][1];
  matrix(Acts::ePHI, Acts::ePHI) = sigmaPhi * sigmaPhi;
  matrix(Acts::eTHETA, Acts::eTHETA) = sigmaTheta * sigmaTheta;
  matrix(Acts::eQOP, Acts::eQOP) = sigmaQOverP * sigmaQOverP;
  matrix(Acts::eT, Acts::eT) = sigmaTime;

  return matrix;
}


int PHActsTrkProp::createNodes(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkProp::getNodes(PHCompositeNode* topNode)
{

  m_actsGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  if (!m_actsGeometry)
  {
    std::cout << "ActsGeometry not on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_actsTracks = findNode::getClass<std::vector<ActsTrack>>(topNode, "ActsTracks");
  
  if (!m_actsTracks)
    {
      std::cout << "ActsTracks not on node tree. Exiting."

		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }


  return Fun4AllReturnCodes::EVENT_OK;
}
