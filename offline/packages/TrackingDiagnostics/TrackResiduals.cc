
#include "TrackResiduals.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/TrackSeed.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

namespace
{
  template <class T>
  inline T square(const T& t)
  {
    return t * t;
  }
  template <class T>
  inline T r(const T& x, const T& y)
  {
    return std::sqrt(square(x) + square(y));
  }

  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track)
  {
    std::vector<TrkrDefs::cluskey> out;
    for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
    {
      if (seed)
      {
        std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
      }
    }

    return out;
  }
}  // namespace

//____________________________________________________________________________..
TrackResiduals::TrackResiduals(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
TrackResiduals::~TrackResiduals()
{
}

//____________________________________________________________________________..
int TrackResiduals::Init(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrackResiduals::InitRun(PHCompositeNode*)
{
  m_outfile = new TFile(m_outfileName.c_str(), "RECREATE");
  createBranches();

  return Fun4AllReturnCodes::EVENT_OK;
}
void TrackResiduals::clearClusterStateVectors()
{
  m_cluslx.clear();
  m_cluslz.clear();
  m_cluselx.clear();
  m_cluselz.clear();
  m_clusgx.clear();
  m_clusgy.clear();
  m_clusgz.clear();
  m_cluslayer.clear();
  m_clussize.clear();

  m_statelx.clear();
  m_statelz.clear();
  m_stateelx.clear();
  m_stateelz.clear();
  m_stategx.clear();
  m_stategy.clear();
  m_stategz.clear();
  m_statepx.clear();
  m_statepy.clear();
  m_statepz.clear();
}
//____________________________________________________________________________..
int TrackResiduals::process_event(PHCompositeNode* topNode)
{
  auto trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  auto clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  auto geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  auto vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!trackmap or !clustermap or !geometry or !vertexmap)
  {
    std::cout << "Missing node, can't continue" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  ActsTransformations transformer;

  if (Verbosity() > 1)
  {
    std::cout << "Track map size is " << trackmap->size() << std::endl;
  }

  for (const auto& [key, track] : *trackmap)
  {
    if (!track)
    {
      continue;
    }
    m_trackid = key;
    m_crossing = track->get_crossing();
    m_px = track->get_px();
    m_py = track->get_py();
    m_pz = track->get_pz();
    m_pt = std::sqrt(square(m_px) + square(m_py));
    m_eta = atanh(m_pz / std::sqrt(square(m_pt) + square(m_pz)));
    m_phi = atan2(m_py, m_px);
    float CVxx = track->get_error(3, 3);
    float CVxy = track->get_error(3, 4);
    float CVyy = track->get_error(4, 4);
    m_deltapt = std::sqrt((CVxx * square(m_px) + 2 * CVxy * m_px * m_py + CVyy * square(m_py)) / (square(m_px) + square(m_py)));

    m_charge = track->get_charge();
    m_quality = track->get_quality();
    m_chisq = track->get_chisq();
    m_ndf = track->get_ndf();
    m_nmaps = 0;
    m_nintt = 0;
    m_ntpc = 0;
    m_nmms = 0;
    m_vertexid = track->get_vertex_id();
    auto vertex = vertexmap->find(m_vertexid)->second;
    if (vertex)
    {
      m_vx = vertex->get_x();
      m_vy = vertex->get_y();
      m_vz = vertex->get_z();
    }

    m_pcax = track->get_x();
    m_pcay = track->get_y();
    m_pcaz = track->get_z();

    clearClusterStateVectors();
    if (Verbosity() > 1)
    {
      std::cout << "Track " << key << " has cluster/states"
                << std::endl;
    }

    for (const auto& ckey : get_cluster_keys(track))
    {
      TrkrCluster* cluster = clustermap->findCluster(ckey);
      switch (TrkrDefs::getTrkrId(ckey))
      {
      case TrkrDefs::mvtxId:
        m_nmaps++;
        break;
      case TrkrDefs::inttId:
        m_nintt++;
        break;
      case TrkrDefs::tpcId:
        m_ntpc++;
        break;
      case TrkrDefs::micromegasId:
        m_nmms++;
        break;
      }

      Acts::Vector3 clusglob = geometry->getGlobalPosition(ckey, cluster);

      auto matched_state = track->begin_states();
      float drmin = -1;
      float clusr = r(clusglob.x(), clusglob.y());

      for (auto state_iter = track->begin_states();
           state_iter != track->end_states();
           ++state_iter)
      {
        SvtxTrackState* state = state_iter->second;
        float stater = r(state->get_x(), state->get_y());
        float thisdr = std::abs(clusr - stater);
        if (drmin < 0 or thisdr < drmin)
        {
          matched_state = state_iter;
          drmin = thisdr;
        }
        else
        {
          break;
        }
      }

      SvtxTrackState* state = matched_state->second;

      //! have cluster and state, fill vectors
      m_cluslx.push_back(cluster->getLocalX());
      float clusz = cluster->getLocalY();
      if(TrkrDefs::getTrkrId(ckey) == TrkrDefs::TrkrId::tpcId)
	{
	  clusz = convertTimeToZ(geometry,ckey, cluster);
	}
      m_cluslz.push_back(clusz);
      m_cluselx.push_back(cluster->getRPhiError());
      m_cluselz.push_back(cluster->getZError());
      m_clusgx.push_back(clusglob.x());
      m_clusgy.push_back(clusglob.y());
      m_clusgz.push_back(clusglob.z());
      m_cluslayer.push_back(TrkrDefs::getLayer(ckey));
      m_clussize.push_back(cluster->getPhiSize() + cluster->getZSize());

      if (Verbosity() > 1)
      {
        std::cout << "Track state/clus in layer "
                  << TrkrDefs::getLayer(ckey) << std::endl;
      }

      auto surf = geometry->maps().getSurface(ckey, cluster);
      Acts::Vector3 stateglob(state->get_x(), state->get_y(), state->get_z());
      Acts::Vector2 stateloc;
      auto norm = surf->normal(geometry->geometry().getGeoContext());
      auto result = surf->globalToLocal(geometry->geometry().getGeoContext(),
                                        stateglob * Acts::UnitConstants::cm,
                                        norm);
      if (result.ok())
      {
        stateloc = result.value() / Acts::UnitConstants::cm;
      }
      else
      {
        //! manual transform for tpc
        Acts::Vector3 loct = surf->transform(geometry->geometry().getGeoContext()).inverse() * (stateglob * Acts::UnitConstants::cm);
        loct /= Acts::UnitConstants::cm;
        stateloc(0) = loct(0);
        stateloc(1) = loct(1);
      }

      const Acts::BoundSymMatrix actscov =
          transformer.rotateSvtxTrackCovToActs(state);

      m_statelx.push_back(stateloc(0));
      m_statelz.push_back(stateloc(1));
      m_stateelx.push_back(std::sqrt(actscov(Acts::eBoundLoc0, Acts::eBoundLoc0)) / Acts::UnitConstants::cm);
      m_stateelz.push_back(std::sqrt(actscov(Acts::eBoundLoc1, Acts::eBoundLoc1)) / Acts::UnitConstants::cm);
      m_stategx.push_back(state->get_x());
      m_stategy.push_back(state->get_y());
      m_stategz.push_back(state->get_z());
      m_statepx.push_back(state->get_px());
      m_statepy.push_back(state->get_py());
      m_statepz.push_back(state->get_pz());
      m_statepl.push_back(state->get_pathlength());
    }

    m_nhits = m_nmaps + m_nintt + m_ntpc + m_nmms;

    m_tree->Fill();
  }

  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}

float TrackResiduals::convertTimeToZ(ActsGeometry* geometry, TrkrDefs::cluskey cluster_key, TrkrCluster *cluster)
{
  // must convert local Y from cluster average time of arival to local cluster z position
  double drift_velocity = geometry->get_drift_velocity();
  double zdriftlength = cluster->getLocalY() * drift_velocity;
  double surfCenterZ = 52.89; // 52.89 is where G4 thinks the surface center is
  double zloc = surfCenterZ - zdriftlength;   // converts z drift length to local z position in the TPC in north
  unsigned int side = TpcDefs::getSide(cluster_key);
  if(side == 0) zloc = -zloc;
  float z = zloc;  // in cm
 
  return z; 
}

//____________________________________________________________________________..
int TrackResiduals::End(PHCompositeNode*)
{
  m_outfile->cd();
  m_tree->Write();
  m_outfile->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackResiduals::createBranches()
{
  m_tree = new TTree("residualtree", "A tree with track, cluster, and state info");
  m_tree->Branch("trackid", &m_trackid, "m_trackid/I");
  m_tree->Branch("crossing", &m_crossing, "m_crossing/I");
  m_tree->Branch("px", &m_px, "m_px/F");
  m_tree->Branch("py", &m_py, "m_py/F");
  m_tree->Branch("pz", &m_pz, "m_pz/F");
  m_tree->Branch("pt", &m_pt, "m_pt/F");
  m_tree->Branch("eta", &m_eta, "m_eta/F");
  m_tree->Branch("phi", &m_phi, "m_phi/F");
  m_tree->Branch("deltapt", &m_deltapt, "m_deltapt/F");
  m_tree->Branch("charge", &m_charge, "m_charge/I");
  m_tree->Branch("quality", &m_quality, "m_quality/F");
  m_tree->Branch("ndf", &m_ndf, "m_ndf/F");
  m_tree->Branch("nhits", &m_nhits, "m_nhits/I");
  m_tree->Branch("nmaps", &m_nmaps, "m_nmaps/I");
  m_tree->Branch("nintt", &m_nintt, "m_nintt/I");
  m_tree->Branch("ntpc", &m_ntpc, "m_ntpc/I");
  m_tree->Branch("nmms", &m_nmms, "m_nmms/I");
  m_tree->Branch("vertexid", &m_vertexid, "m_vertexid/I");
  m_tree->Branch("vx", &m_vx, "m_vx/F");
  m_tree->Branch("vy", &m_vy, "m_vy/F");
  m_tree->Branch("vz", &m_vz, "m_vz/F");
  m_tree->Branch("pcax", &m_pcax, "m_pcax/F");
  m_tree->Branch("pcay", &m_pcay, "m_pcay/F");
  m_tree->Branch("pcaz", &m_pcaz, "m_pcaz/F");

  m_tree->Branch("cluslx", &m_cluslx);
  m_tree->Branch("cluslz", &m_cluslz);
  m_tree->Branch("cluselx", &m_cluselx);
  m_tree->Branch("cluselz", &m_cluselz);
  m_tree->Branch("clusgx", &m_clusgx);
  m_tree->Branch("clusgy", &m_clusgy);
  m_tree->Branch("clusgz", &m_clusgz);
  m_tree->Branch("cluslayer", &m_cluslayer);
  m_tree->Branch("clussize", &m_clussize);

  m_tree->Branch("statelx", &m_statelx);
  m_tree->Branch("statelz", &m_statelz);
  m_tree->Branch("stateelx", &m_stateelx);
  m_tree->Branch("stateelz", &m_stateelz);
  m_tree->Branch("stategx", &m_stategx);
  m_tree->Branch("stategy", &m_stategy);
  m_tree->Branch("stategz", &m_stategz);
  m_tree->Branch("statepx", &m_statepx);
  m_tree->Branch("statepy", &m_statepy);
  m_tree->Branch("statepz", &m_statepz);
  m_tree->Branch("statepl", &m_statepl);
}
