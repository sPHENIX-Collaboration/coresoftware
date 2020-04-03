
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

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/Fitting/TrkrClusterFittingAlgorithm.hpp>
#include <ACTFW/Framework/AlgorithmContext.hpp>

#include <cmath>
#include <iostream>
#include <vector>
#include <utility>
#include <TMatrixDSym.h>

PHActsTrkProp::PHActsTrkProp(const std::string& name)
  : PHTrackPropagating(name)
  , m_event(0)
  , m_fitCfgOptions(nullptr)

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


  for(SvtxTrackMap::Iter trackIter = _track_map->begin();
      trackIter != _track_map->end(); ++trackIter)
    {
      const SvtxTrack *track = trackIter->second;

      if(!track) 
	continue;

      if(Verbosity() > 1)
	{
	  std::cout << "Found SvtxTrack " << trackIter->first << std::endl;
	  track->identify();
	}

      const double x = track->get_x();
      const double y = track->get_y();
      const double d0 = sqrt(x*x + y*y);
      const double z0 = track->get_z();
      const double phi = track->get_phi();
      const double theta = 2.*atan(exp(-1*track->get_eta()));
      const int charge = track->get_charge();
      const double qop = (double)charge / track->get_p();
      /// Just set to 0 for now?
      const double time = 0;

      Acts::BoundSymMatrix cov = getActsCovMatrix(track);

      Acts::BoundVector pars;
      pars << d0, z0, phi, theta, qop, time;

      Acts::Vector3D sPosition(0., 0., 0.);
      Acts::Vector3D sMomentum(0., 0., 0.);
      Acts::BoundParameters startParameters(
          m_fitCfgOptions->geoContext, 
	  std::move(cov), std::move(pars), surface);
      sPosition = startParameters.position();
      sMomentum = startParameters.momentum();


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




Acts::BoundSymMatrix PHActsTrkProp::getActsCovMatrix(const SvtxTrack *track)
{
  Acts::BoundSymMatrix matrix = Acts::BoundSymMatrix::Zero();
  const double px = track->get_px();
  const double py = track->get_py();
  const double pz = track->get_pz();
  const double p = sqrt(px * px + py * py + pz * pz);
  const double x = track->get_x();
  const double y = track->get_y();

  const double d = sqrt(x*x + y*y);
  
  // Get the track seed covariance matrix
  // These are the variances, so the std devs are sqrt(seed_cov[i][j])
  TMatrixDSym cov(6);
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      cov[i][j] = track->get_error(i, j);
    }
  }
  const double sigmad = (1./d) * sqrt( x*x*cov[0][0] + y*y*cov[1][1]);
  const double sigmap = sqrt(px * px * cov[3][3] + py * py * cov[4][4] + pz * pz * cov[5][5]) / p;

  // Need to convert cov from x,y,z,px,py,pz basis to Acts basis of
  // x,y,phi/theta of p, qoverp, time
  double phi = track->get_phi();

  const double pxfracerr = cov[3][3] / (px * px);
  const double pyfracerr = cov[4][4] / (py * py);
  const double phiPrefactor = fabs(py) / (fabs(px) * (1 + (py / px) * (py / px)));
  const double sigmaPhi = phi * phiPrefactor * sqrt(pxfracerr + pyfracerr);
  const double theta = acos(pz / p);
  const double thetaPrefactor = ((fabs(pz)) / (p * sqrt(1 - (pz / p) * (pz / p))));
  const double sigmaTheta = thetaPrefactor * sqrt(sigmap * sigmap / (p * p) + cov[5][5] / (pz * pz));
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
        std::cout << cov[i][j] << ", ";
      }
      std::cout << std::endl;
    }
    std::cout << "Corresponding uncertainty calculations: " << std::endl;
    std::cout << "d: " << d << " +/- " << sigmad << std::endl;
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
  /// Rest are std devs, so need to square them
  matrix(Acts::eLOC_0, Acts::eLOC_0) = sigmad * sigmad;
  matrix(Acts::eLOC_1, Acts::eLOC_1) = cov[2][2]; // z
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

  m_fitCfgOptions = findNode::getClass<FitCfgOptions>(topNode, "ActsFitCfg");

  if (!m_fitCfgOptions)
  {
    std::cout << "Acts FitCfgOptions not on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
