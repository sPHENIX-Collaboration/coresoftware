#include "KFParticle_truthAndDetTools.h"

#include <intt/InttDefs.h>

#include <mvtx/MvtxDefs.h>

#include <tpc/TpcDefs.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTruthEval.h>
#include <g4eval/SvtxVertexEval.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phool/getClass.h>

#include <KFParticle.h>

#include <TTree.h>

std::map<std::string, int> Use =
    {
        {"MVTX", 1},
        {"INTT", 1},
        {"TPC", 1},
        {"EMCAL", 0},
        {"OHCAL", 0},
        {"IHCAL", 0}};

KFParticle_truthAndDetTools::KFParticle_truthAndDetTools()
  : m_svtx_evalstack(nullptr)
{
}  //Constructor

KFParticle_truthAndDetTools::~KFParticle_truthAndDetTools() {}  //Destructor

SvtxTrack *KFParticle_truthAndDetTools::getTrack(unsigned int track_id, SvtxTrackMap *trackmap)
{
  SvtxTrack *matched_track = NULL;

  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter)
  {
    if (iter->first == track_id) matched_track = iter->second;
  }

  return matched_track;
}

SvtxVertex *KFParticle_truthAndDetTools::getVertex(unsigned int vertex_id, SvtxVertexMap *vertexmap)
{
  SvtxVertex *matched_vertex = NULL;

  for (SvtxVertexMap::Iter iter = vertexmap->begin();
       iter != vertexmap->end();
       ++iter)
  {
    if (iter->first == vertex_id) matched_vertex = iter->second;
  }

  return matched_vertex;
}

void KFParticle_truthAndDetTools::initializeTruthBranches(TTree *m_tree, int daughter_id, std::string daughter_number, bool m_constrain_to_vertex_truthMatch)
{
  m_tree->Branch(TString(daughter_number) + "_true_vertex_x", &m_true_daughter_vertex_x[daughter_id], TString(daughter_number) + "_true_vertex_x/F");
  m_tree->Branch(TString(daughter_number) + "_true_vertex_y", &m_true_daughter_vertex_y[daughter_id], TString(daughter_number) + "_true_vertex_y/F");
  m_tree->Branch(TString(daughter_number) + "_true_vertex_z", &m_true_daughter_vertex_z[daughter_id], TString(daughter_number) + "_true_vertex_z/F");
  if (m_constrain_to_vertex_truthMatch) m_tree->Branch(TString(daughter_number) + "_true_IP", &m_true_daughter_ip[daughter_id], TString(daughter_number) + "_true_IP/F");
  if (m_constrain_to_vertex_truthMatch) m_tree->Branch(TString(daughter_number) + "_true_IP_xy", &m_true_daughter_ip_xy[daughter_id], TString(daughter_number) + "_true_IP_xy/F");
  m_tree->Branch(TString(daughter_number) + "_true_px", &m_true_daughter_px[daughter_id], TString(daughter_number) + "_true_px/F");
  m_tree->Branch(TString(daughter_number) + "_true_py", &m_true_daughter_py[daughter_id], TString(daughter_number) + "_true_py/F");
  m_tree->Branch(TString(daughter_number) + "_true_pz", &m_true_daughter_pz[daughter_id], TString(daughter_number) + "_true_pz/F");
  m_tree->Branch(TString(daughter_number) + "_true_p", &m_true_daughter_p[daughter_id], TString(daughter_number) + "_true_p/F");
  m_tree->Branch(TString(daughter_number) + "_true_pT", &m_true_daughter_pt[daughter_id], TString(daughter_number) + "_true_pT/F");
  m_tree->Branch(TString(daughter_number) + "_true_ID", &m_true_daughter_id[daughter_id], TString(daughter_number) + "_true_ID/I");
  if (m_constrain_to_vertex_truthMatch)
  {
    m_tree->Branch(TString(daughter_number) + "_true_PV_x", &m_true_daughter_pv_x[daughter_id], TString(daughter_number) + "_true_pv_x/F");
    m_tree->Branch(TString(daughter_number) + "_true_PV_y", &m_true_daughter_pv_y[daughter_id], TString(daughter_number) + "_true_pv_y/F");
    m_tree->Branch(TString(daughter_number) + "_true_PV_z", &m_true_daughter_pv_z[daughter_id], TString(daughter_number) + "_true_pv_x/F");
  }
}

void KFParticle_truthAndDetTools::fillTruthBranch(PHCompositeNode *topNode, TTree *m_tree, KFParticle daughter, int daughter_id, KFParticle vertex, bool m_constrain_to_vertex_truthMatch)
{
  float true_px, true_py, true_pz, true_p, true_pt;

  if (!m_svtx_evalstack)
  {
    m_svtx_evalstack = new SvtxEvalStack(topNode);
    //trackeval = m_svtx_evalstack->get_track_eval();
    clustereval = m_svtx_evalstack->get_cluster_eval();
    trutheval = m_svtx_evalstack->get_truth_eval();
    vertexeval = m_svtx_evalstack->get_vertex_eval();
  }
  //m_svtx_evalstack->next_event(topNode);

  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode*>(nodeIter.findFirst("SvtxTrackMap"));
  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  }
  else
  {
    std::cout << "KFParticle truth matching: SvtxTrackMap does not exist" << std::endl;
  }
  findNode = dynamic_cast<PHNode*>(nodeIter.findFirst("SvtxVertexMap"));
  if (findNode)
  {
    dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  }
  else
  {
    std::cout << "KFParticle truth matching: SvtxVertexMap does not exist" << std::endl;
  }

  m_svtx_evalstack->next_event(topNode);

  track = getTrack(daughter.Id(), dst_trackmap);

  TrkrDefs::cluskey clusKey = *track->begin_cluster_keys();
  g4particle = clustereval->max_truth_particle_by_cluster_energy(clusKey);

  bool isParticleValid = g4particle == nullptr ? 0 : 1;
  
  true_px = isParticleValid ? (Float_t) g4particle->get_px() : 0.;
  true_py = isParticleValid ? (Float_t) g4particle->get_py() : 0.;
  true_pz = isParticleValid ? (Float_t) g4particle->get_pz() : 0.;
  true_p = sqrt(pow(true_px, 2) + pow(true_py, 2) + pow(true_pz, 2));
  true_pt = sqrt(pow(true_px, 2) + pow(true_py, 2));

  m_true_daughter_px[daughter_id] = true_px;
  m_true_daughter_py[daughter_id] = true_py;
  m_true_daughter_pz[daughter_id] = true_pz;
  m_true_daughter_p[daughter_id] = true_p;
  m_true_daughter_pt[daughter_id] = true_pt;
  m_true_daughter_id[daughter_id] = isParticleValid ? g4particle->get_pid() : 0;

  if (isParticleValid)
  {
    g4vertex_point = trutheval->get_vertex(g4particle);
  }

  m_true_daughter_vertex_x[daughter_id] = isParticleValid ? g4vertex_point->get_x() : 0.;
  m_true_daughter_vertex_y[daughter_id] = isParticleValid ? g4vertex_point->get_y() : 0.;
  m_true_daughter_vertex_z[daughter_id] = isParticleValid ? g4vertex_point->get_z() : 0.;

  if (m_constrain_to_vertex_truthMatch)  
  {
    //Calculate true DCA
    SvtxVertex* recoVertex = getVertex(vertex.Id(), dst_vertexmap);
    PHG4VtxPoint* truePoint = vertexeval->max_truth_point_by_ntracks(recoVertex);

    KFParticle trueKFParticleVertex;

    float f_vertexParameters[6] = {0};

    if (truePoint == NULL)
    {
      std::cout << "KFParticle truth matching: This event has no PHG4VtxPoint information!\n";
      std::cout << "Your truth track DCA will be measured wrt a reconstructed vertex!" << std::endl; 

      f_vertexParameters[0] = recoVertex->get_x(); 
      f_vertexParameters[1] = recoVertex->get_y();
      f_vertexParameters[2] = recoVertex->get_z();
    }
    else
    {
      f_vertexParameters[0] = truePoint->get_x();
      f_vertexParameters[1] = truePoint->get_y();
      f_vertexParameters[2] = truePoint->get_z();
    }

    float f_vertexCovariance[21] = {0};

    trueKFParticleVertex.Create(f_vertexParameters, f_vertexCovariance, 0, -1);

    KFParticle trueKFParticle;

    float f_trackParameters[6] = {m_true_daughter_vertex_x[daughter_id],
                                  m_true_daughter_vertex_y[daughter_id],
                                  m_true_daughter_vertex_z[daughter_id],
                                  true_px,
                                  true_py,
                                  true_pz};

    float f_trackCovariance[21] = {0};

    trueKFParticle.Create(f_trackParameters, f_trackCovariance, 1, -1);

    m_true_daughter_ip[daughter_id] = trueKFParticle.GetDistanceFromVertex(trueKFParticleVertex);
    m_true_daughter_ip_xy[daughter_id] = trueKFParticle.GetDistanceFromVertexXY(trueKFParticleVertex);

    m_true_daughter_pv_x[daughter_id] = truePoint == NULL ? -99. : truePoint->get_x();
    m_true_daughter_pv_y[daughter_id] = truePoint == NULL ? -99. : truePoint->get_y();
    m_true_daughter_pv_z[daughter_id] = truePoint == NULL ? -99. : truePoint->get_z();

  }
}

void KFParticle_truthAndDetTools::initializeDetectorBranches(TTree *m_tree, int daughter_id, std::string daughter_number)
{
  m_tree->Branch(TString(daughter_number) + "_local_x", &detector_local_x[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_local_y", &detector_local_y[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_local_z", &detector_local_z[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_layer", &detector_layer[daughter_id]);

  for (auto const &subdetector : Use)
  {
    if (subdetector.second) initializeSubDetectorBranches(m_tree, subdetector.first, daughter_id, daughter_number);
  }
}

void KFParticle_truthAndDetTools::initializeSubDetectorBranches(TTree *m_tree, std::string detectorName, int daughter_id, std::string daughter_number)
{
  if (detectorName == "MVTX")
  {
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_staveID", &mvtx_staveID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_chipID", &mvtx_chipID[daughter_id]);
  }
  if (detectorName == "INTT")
  {
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_ladderZID", &intt_ladderZID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_ladderPhiID", &intt_ladderPhiID[daughter_id]);
  }
  if (detectorName == "TPC")
  {
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_sectorID", &tpc_sectorID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_side", &tpc_side[daughter_id]);
  }
}

void KFParticle_truthAndDetTools::fillDetectorBranch(PHCompositeNode *topNode,
                                                     TTree *m_tree, KFParticle daughter, int daughter_id)
{
  PHNodeIterator nodeIter(topNode);

  PHNode *findNode = dynamic_cast<PHNode*>(nodeIter.findFirst("SvtxTrackMap"));
  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  }
  else
  {
    std::cout << "KFParticle detector info: SvtxTrackMap does not exist" << std::endl;
  }

  findNode = dynamic_cast<PHNode*>(nodeIter.findFirst("TRKR_CLUSTER"));
  if (findNode)
  {
    dst_clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  }
  else
  {
    std::cout << "KFParticle detector info: TRKR_CLUSTER does not exist" << std::endl;
  }

  track = getTrack(daughter.Id(), dst_trackmap);

  for (SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys();
       iter != track->end_cluster_keys();
       ++iter)
  {
    TrkrDefs::cluskey clusKey = *iter;
    TrkrCluster *cluster = dst_clustermap->findCluster(clusKey);
    const unsigned int trkrId = TrkrDefs::getTrkrId(clusKey);

    detector_local_x[daughter_id].push_back(cluster->getX());
    detector_local_y[daughter_id].push_back(cluster->getY());
    detector_local_z[daughter_id].push_back(cluster->getZ());
    detector_layer[daughter_id].push_back(TrkrDefs::getLayer(clusKey));
    unsigned int staveId, chipId, ladderZId, ladderPhiId, sectorId, side;
    staveId = chipId = ladderZId = ladderPhiId = sectorId = side = -99;

    if (Use["MVTX"] && trkrId == TrkrDefs::mvtxId)
    {
      staveId = MvtxDefs::getStaveId(clusKey);
      chipId = MvtxDefs::getChipId(clusKey);
    }
    else if (Use["INTT"] && trkrId == TrkrDefs::inttId)
    {
      ladderZId = InttDefs::getLadderZId(clusKey);
      ladderPhiId = InttDefs::getLadderPhiId(clusKey);
    }
    else if (Use["TPC"] && trkrId == TrkrDefs::tpcId)
    {
      sectorId = TpcDefs::getSectorId(clusKey);
      side = TpcDefs::getSide(clusKey);
    }

    mvtx_staveID[daughter_id].push_back(staveId);
    mvtx_chipID[daughter_id].push_back(chipId);
    intt_ladderZID[daughter_id].push_back(ladderZId);
    intt_ladderPhiID[daughter_id].push_back(ladderPhiId);
    tpc_sectorID[daughter_id].push_back(sectorId);
    tpc_side[daughter_id].push_back(side);
  }
}
