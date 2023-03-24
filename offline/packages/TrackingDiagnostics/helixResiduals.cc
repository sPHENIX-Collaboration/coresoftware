#include "helixResiduals.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <trackbase/TrkrDefs.h>

#include <TFile.h>
#include <TNtuple.h>

helixResiduals::helixResiduals(const std::string &name)
  : SubsysReco(name)
{
  _fitter = new HelicalFitter;
}

int helixResiduals::InitRun(PHCompositeNode *topNode)
{
  _fitter->InitRun(topNode);
  const char *cfilepath = filepath.c_str();
  fout = new TFile(cfilepath, "recreate");
  ntp_residuals = new TNtuple("ntp_residuals", "Seed Residuals", "seed_id:layer:dphi:dz:x:y:z:pt:px:py:pz:crossing:isSilicon:isTpc");

  getNodes(topNode);

  return 0;
}

int helixResiduals::End(PHCompositeNode * /**topNode*/)
{
  fout->cd();
  ntp_residuals->Write();
  fout->Close();

  delete _fitter;

  return 0;
}

int helixResiduals::getNodes(PHCompositeNode *topNode)
{
  _tracks = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_tracks)
  {
    std::cerr << PHWHERE << "No SvtxTrackMap on node tree, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _tpc_seeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!_tpc_seeds)
  {
    std::cerr << PHWHERE << " ERROR: Can't find TpcTrackSeedContainer " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _si_seeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!_si_seeds)
  {
    std::cerr << PHWHERE << " ERROR: Can't find SiliconTrackSeedContainer" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _clusters = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_clusters)
  {
    std::cerr << PHWHERE << " ERROR: Can't find TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tGeometry)
  {
    std::cerr << PHWHERE << "No acts tracking geometry, can't proceed" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int helixResiduals::process_event(PHCompositeNode * /*topNode*/)
{
  for (auto tpcseed_iter = _tpc_seeds->begin(); tpcseed_iter != _tpc_seeds->end(); ++tpcseed_iter)
  {
    int id = _tpc_seeds->index(tpcseed_iter);
    TrackSeed *tpcseed = _tpc_seeds->get(id);
    if (!tpcseed)
    {
      continue;
    }
    std::cout << "processing tpc seed " << id << std::endl;
    fill_residuals(tpcseed, id, true);
  }
  for (auto siseed_iter = _si_seeds->begin(); siseed_iter != _si_seeds->end(); ++siseed_iter)
  {
    int id = _si_seeds->index(siseed_iter);
    TrackSeed *siseed = _si_seeds->get(id);
    if (!siseed)
    {
      continue;
    }
    std::cout << "processing si seed " << id << std::endl;
    fill_residuals(siseed, id, false);
  }
  std::cout << "done" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

void helixResiduals::fill_residuals(TrackSeed *seed, int seed_id, bool isTpc)
{
  if (seed->size_cluster_keys() == 0)
  {
    return;
  }

  std::vector<Acts::Vector3> clusterPositions;
  std::vector<TrkrDefs::cluskey> clusterKeys;
  _fitter->getTrackletClusters(seed, clusterPositions, clusterKeys);
  std::vector<float> fitparams = _fitter->fitClusters(clusterPositions, clusterKeys);

  float pt = seed->get_pt();
  float px = seed->get_px(_clusters, tGeometry);
  float py = seed->get_py(_clusters, tGeometry);
  float pz = seed->get_pz();
  unsigned int crossing = seed->get_crossing();

  for (size_t i = 0; i < clusterPositions.size(); i++)
  {
    unsigned int layer = TrkrDefs::getLayer(clusterKeys[i]);
    Acts::Vector3 position = clusterPositions[i];
    Acts::Vector3 pca = _fitter->get_helix_pca(fitparams, position);
    float cluster_phi = atan2(position(1), position(0));
    float pca_phi = atan2(pca(1), pca(0));
    float dphi = cluster_phi - pca_phi;
    if (dphi > M_PI)
    {
      dphi = 2 * M_PI - dphi;
    }
    if (dphi < -M_PI)
    {
      dphi = 2 * M_PI + dphi;
    }
    float dz = position(2) - pca(2);

    ntp_residuals->Fill(seed_id, layer, dphi, dz, position(0), position(1), position(2), pt, px, py, pz, crossing, !isTpc, isTpc);
  }
}
