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

#include <Acts/EventData/ChargePolicy.hpp>
#include <Acts/EventData/SingleCurvilinearTrackParameters.hpp>

/// Tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

/// std (and the like) includes
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

#include <TMatrixDSym.h>

PHActsTracks::PHActsTracks(const std::string &name)
  : SubsysReco(name)
  , m_actsProtoTracks(nullptr)
  , m_trackMap(nullptr)
  , m_hitIdClusKey(nullptr)
  , m_sourceLinks(nullptr)
{
  Verbosity(0);
}

int PHActsTracks::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 10)
  {
    std::cout << "Finished PHActsTracks" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTracks::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTracks::InitRun(PHCompositeNode *topNode)
{
  createNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTracks::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "Start process_event in PHActsTracks" << std::endl;
  }

  /// Check to get the nodes needed
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  /// Vector to hold source links for a particular track
  std::vector<SourceLink> trackSourceLinks;
  std::vector<FW::TrackParameters> trackSeeds;

  for (SvtxTrackMap::Iter trackIter = m_trackMap->begin();
       trackIter != m_trackMap->end(); ++trackIter)
  {
    const SvtxTrack *track = trackIter->second;

    if (!track)
      continue;

    if (Verbosity() > 1)
    {
      std::cout << "found SvtxTrack " << trackIter->first << std::endl;
      track->identify();
    }

    /// Get the necessary parameters and values for the TrackParameters
    const Acts::BoundSymMatrix seedCov = getActsCovMatrix(track);
    const Acts::Vector3D seedPos(track->get_x(),
                                 track->get_y(),
                                 track->get_z());
    const Acts::Vector3D seedMom(track->get_px(),
                                 track->get_py(),
                                 track->get_pz());

    // just set to 0 for now?
    const double trackTime = 0;
    const int trackQ = track->get_charge();

    const FW::TrackParameters trackSeed(seedCov, seedPos,
                                        seedMom, trackQ, trackTime);

    /// Start fresh for this track
    trackSourceLinks.clear();
    for (SvtxTrack::ConstClusterKeyIter clusIter = track->begin_cluster_keys();
         clusIter != track->end_cluster_keys();
         ++clusIter)
    {
      const TrkrDefs::cluskey key = *clusIter;

      const unsigned int hitId = m_hitIdClusKey->find(key)->second;

      if (Verbosity() > 0)
      {
        std::cout << "cluskey " << key
                  << " has hitid " << hitId
                  << std::endl;
      }
      trackSourceLinks.push_back(m_sourceLinks->find(hitId)->second);
    }

    if (Verbosity() > 0)
    {
      for (unsigned int i = 0; i < trackSourceLinks.size(); ++i)
      {
        std::cout << "proto_track readback: hitid " << trackSourceLinks.at(i).hitID() << std::endl;
      }
    }

    ActsTrack actsTrack(trackSeed, trackSourceLinks);
    m_actsProtoTracks->push_back(actsTrack);
  }

  if (Verbosity() > 20)
    std::cout << "Finished PHActsTrack::process_event" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::BoundSymMatrix PHActsTracks::getActsCovMatrix(const SvtxTrack *track)
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

void PHActsTracks::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsTracks::createNodes");
  }

  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_actsProtoTracks = findNode::getClass<std::vector<ActsTrack>>(topNode, "ActsProtoTracks");

  if (!m_actsProtoTracks)
  {
    m_actsProtoTracks = new std::vector<ActsTrack>;

    PHDataNode<std::vector<ActsTrack>> *protoTrackNode =
        new PHDataNode<std::vector<ActsTrack>>(m_actsProtoTracks, "ActsProtoTracks");

    svtxNode->addNode(protoTrackNode);
  }

  return;
}

int PHActsTracks::getNodes(PHCompositeNode *topNode)
{
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  if (!m_trackMap)
  {
    std::cout << PHWHERE << "SvtxTrackMap not found on node tree. Exiting."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_sourceLinks = findNode::getClass<std::map<unsigned int, SourceLink>>(topNode, "TrkrClusterSourceLinks");

  if (!m_sourceLinks)
  {
    std::cout << PHWHERE << "TrkrClusterSourceLinks node not found on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_hitIdClusKey = findNode::getClass<std::map<TrkrDefs::cluskey, unsigned int>>(topNode, "HitIDClusIDActsMap");

  if (!m_hitIdClusKey)
  {
    std::cout << PHWHERE << "HitID cluster key map not found on node tree. Exiting. "
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
