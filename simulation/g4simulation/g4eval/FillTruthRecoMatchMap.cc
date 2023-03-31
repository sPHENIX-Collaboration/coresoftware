#include "FillTruthRecoMatchMap.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/phool.h>  // for PHWHERE
#include <phool/getClass.h>
#include <g4tracking/EmbRecoMatch.h>
#include <g4tracking/EmbRecoMatchContainer.h>
#include <trackbase_historic/PHG4ParticleSvtxMap.h>
#include <trackbase_historic/PHG4ParticleSvtxMap_v1.h>
#include <trackbase_historic/SvtxPHG4ParticleMap_v1.h>
#include <iostream>

using std::cout;
using std::endl;
//____________________________________________________________________________..
FillTruthRecoMatchMap::FillTruthRecoMatchMap(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
FillTruthRecoMatchMap::~FillTruthRecoMatchMap()
{
}

//____________________________________________________________________________..
int FillTruthRecoMatchMap::Init(PHCompositeNode *topNode)
{
  if (Verbosity()>1) topNode->print();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int FillTruthRecoMatchMap::InitRun(PHCompositeNode *topNode)
{
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int FillTruthRecoMatchMap::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cout << PHWHERE << " DST node is missing, quitting" << std::endl;
    std::cerr << PHWHERE << " DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in FillTruthRecoMatchMap::createNodes");
  }

  m_EmbRecoMatchContainer = findNode::getClass<EmbRecoMatchContainer>(topNode,"TRKR_EMBRECOMATCHCONTAINER");
  if (!m_EmbRecoMatchContainer) {
    std::cout << PHWHERE << " Cannot find node TRKR_EMBRECOMATCHCONTAINER on node tree; quitting " << std::endl;
    std::cerr << PHWHERE << " Cannot find node TRKR_EMBRECOMATCHCONTAINER on node tree; quitting " << std::endl;
    throw std::runtime_error(" Cannot find node TRKR_EMBRECOMATCHCONTAINER on node tree; quitting");
  }

  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_PHG4ParticleSvtxMap = findNode::getClass<PHG4ParticleSvtxMap>(topNode, "PHG4ParticleToRecoMap");
  if (!m_PHG4ParticleSvtxMap)
  {
    m_PHG4ParticleSvtxMap = new PHG4ParticleSvtxMap_v1;
    PHIODataNode<PHObject> *truthNode =
        new PHIODataNode<PHObject>(m_PHG4ParticleSvtxMap, "PHG4ParticleToRecoMap", "PHObject");
    svtxNode->addNode(truthNode);
  }

  m_SvtxPHG4ParticleMap = findNode::getClass<SvtxPHG4ParticleMap>(topNode, "RecoToPHG4ParticleMap");
  if (!m_SvtxPHG4ParticleMap)
  {
    m_SvtxPHG4ParticleMap = new SvtxPHG4ParticleMap_v1;
    PHIODataNode<PHObject> *recoNode =
        new PHIODataNode<PHObject>(m_SvtxPHG4ParticleMap, "RecoToPHG4ParticleMap", "PHObject");
    svtxNode->addNode(recoNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int FillTruthRecoMatchMap::process_event(PHCompositeNode * /*topNode*/)
{
  if (Verbosity()>5) cout << " FillTruthRecoMatchMap::process_event() " << endl;
  // make maps to all the matches for both truth and reco and fill the 
  /* auto matches = m_EmbRecoMatchContainer->getMatches(); */

  /* map<unsigned short, set<int>> truth_entries; */
  /* map<unsigned short, set<int>> m_reco_entries; */
  /* int cnt {0}; */
  /* for (auto& match : m_EmbRecoMatchContainer->getMatches()) { */
    
  /* } */


  SvtxPHG4ParticleMap::Map map_RtoT {}; // RtoT = "Reco to Truth"
  PHG4ParticleSvtxMap::Map map_TtoR {}; // TtoR = "Truth to Reco"

  for (auto& match : m_EmbRecoMatchContainer->getMatches()) {
    // every match is a unique match from truth to embedded particles
    const int            gtrackID = match->idTruthTrack()   ;
    const unsigned int   id_reco  = match->idRecoTrack()    ;
    const unsigned short n_match  = match->nClustersMatched() ;
    const unsigned short n_truth  = match->nClustersTruth() ;
    const unsigned short n_reco   = match->nClustersReco()  ;

    // fill the map_TtoR
    /* PHG4ParticleSvtxMap::WeightedRecoTrackMap*  entry_TtoR; */
    if (map_TtoR.find(gtrackID) == map_TtoR.end()) {
      map_TtoR[gtrackID] = PHG4ParticleSvtxMap::WeightedRecoTrackMap{};
    }
    auto& entry_TtoR = map_TtoR[gtrackID];
    float weight_TtoR = (float)n_match + (float)n_truth/100.;
    if (entry_TtoR.find(weight_TtoR) == entry_TtoR.end()) {
      entry_TtoR[weight_TtoR] = { id_reco }; // i.e. std::set<unsigned int> { id_reco };
    } else {
      entry_TtoR[weight_TtoR].insert(id_reco);
    }

    // fill the map_RtoT
    /* SvtxPHG4ParticleMap::WeightedTruthTrackMap*  entry_RtoT; */
    if (map_RtoT.find(id_reco) == map_RtoT.end()) {
      map_RtoT[id_reco] = SvtxPHG4ParticleMap::WeightedTruthTrackMap{};
    } 
    auto& entry_RtoT = map_RtoT[id_reco];
    float weight_RtoT = (float)n_match + (float)n_reco/100.;
    if (entry_RtoT.find(weight_RtoT) == entry_RtoT.end()) {
        entry_RtoT[weight_RtoT] = { gtrackID }; // i.e. std::set<int> { gtrackID }
    } else {
      entry_RtoT[weight_RtoT].insert(gtrackID);
    }

    if (Verbosity() > 20) {
      printf("EmbRecoMatch: gtrackID(%2i) id_reco(%2i) nclusters:match(%i),gtrack(%2i),reco(%2i)\n",
          gtrackID, (int)id_reco, (int)n_match, (int)n_truth, (int)n_reco);
      printf("   -> in SvtxPHG4ParticleMap {id_reco->{weight->id_true}} = {%2i->{%5.2f->%2i}}\n", 
          (int)id_reco,  weight_RtoT, (int)gtrackID);
      printf("   -> in PHG4ParticleSvtxMap {id_true->{weight->id_reco}} = {%2i->{%5.2f->%2i}}\n", 
         gtrackID, weight_TtoR, (int)id_reco );
    }
  }

  // fill the output maps
  for (auto& entry : map_TtoR) {
    m_PHG4ParticleSvtxMap->insert(entry.first, entry.second);
  }
  for (auto& entry : map_RtoT) {
    m_SvtxPHG4ParticleMap->insert(entry.first, entry.second);
  }

  if (Verbosity()>15) {
    cout << PHWHERE << " Print out: " << endl << " Contents of SvtxPHG4ParticleMap (node \"RecoToPHG4ParticleMap\")" << endl
                    << " and PHG4ParticleSvtxMap (node \"PHG4ParticleToRecoMap\")" << endl;

  std::cout << " --BEGIN-- Contents of SvtxPHG4ParticleMap: " << std::endl;
    for (auto iter = m_SvtxPHG4ParticleMap->begin(); iter != m_SvtxPHG4ParticleMap->end(); ++iter) {
      printf("    { %2i ", (int)iter->first); // id_reco
      auto n_matches = iter->second.size();
      long unsigned int cnt_matches = 0;
      for (auto matches : iter->second) {
        if (cnt_matches==0) printf("-> { %5.2f -> ", matches.first);
        else                printf("         -> { %5.2f -> ", matches.first);
        auto size = matches.second.size();
        long unsigned int i=0;
        for (auto id_true : matches.second) {
          if (i<size-1) printf("%2i, ",id_true);
          else          printf("%2i }",id_true);
          ++i;
        } // end matches cnt
        if (cnt_matches < (n_matches-1)) printf(",\n");
        else                             printf("}\n");
        ++cnt_matches;
      }
    }
    std::cout << " --END-- Contents of SvtxPHG4ParticleMap: " << std::endl << std::endl;

  std::cout << " --BEGIN-- Contents of PHG4ParticleToRecoMap: " << std::endl;
    for (auto iter = m_PHG4ParticleSvtxMap->begin(); iter != m_PHG4ParticleSvtxMap->end(); ++iter) {
      printf("    { %2i ", iter->first); // id_true
      auto n_matches = iter->second.size();
      long unsigned int cnt_matches = 0;
      for (auto matches : iter->second) {
        if (cnt_matches==0) printf("-> { %5.2f -> ", matches.first);
        else                printf("         -> { %5.2f -> ", matches.first);
        auto size = matches.second.size();
        long unsigned int i=0;
        for (auto id_reco : matches.second) {
          if (i<size-1) printf("%2i, ",id_reco);
          else          printf("%2i }",id_reco);
          ++i;
        } // end matches cnt
        if (cnt_matches < (n_matches-1)) printf(",\n");
        else                             printf("}\n");
        ++cnt_matches;
      }
    }
    std::cout << " --END-- Contents of SvtxPHG4ParticleMap: " << std::endl << std::endl;


  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int FillTruthRecoMatchMap::End(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
