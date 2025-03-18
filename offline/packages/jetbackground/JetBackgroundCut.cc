#include "JetBackgroundCut.h"
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <ffarawobjects/Gl1Packetv2.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <globalvertex/GlobalVertexMapv1.h>
#include <jetbase/Jet.h>
#include <jetbase/JetContainerv1.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phparameter/PHParameters.h>
#include <cmath>

//____________________________________________________________________________..
JetBackgroundCut::JetBackgroundCut(const std::string &jetNodeName, const std::string &name, const int debug, const bool doAbort, GlobalVertex::VTXTYPE vtxtype, int sysvar)
  : SubsysReco(name)
  , _name(name)
  , _jetNodeName(jetNodeName)
  , _vtxtype(vtxtype)
  , _cutParams(name)
{
  _debug = debug;
  _doAbort = doAbort;
  _sysvar = sysvar;
  SetDefaultParams();
}

//____________________________________________________________________________..
JetBackgroundCut::~JetBackgroundCut() = default;

//____________________________________________________________________________..
int JetBackgroundCut::Init(PHCompositeNode *topNode)
{
  CreateNodeTree(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

void JetBackgroundCut::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  if (!parNode)
  {
    std::cout << "No RUN node found; cannot create PHParameters for storing cut results. Aborting run!";
  }

  _cutParams.SaveToNodeTree(parNode, "JetCutParams");
}

//____________________________________________________________________________..
int JetBackgroundCut::process_event(PHCompositeNode *topNode)
{
  TowerInfoContainer *towersEM = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  TowerInfoContainer *towersOH = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  JetContainer *jets = findNode::getClass<JetContainerv1>(topNode, _jetNodeName);
  GlobalVertexMap *gvtxmap = findNode::getClass<GlobalVertexMapv1>(topNode, "GlobalVertexMap");

  RawTowerGeomContainer *geom[2];
  geom[0] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  geom[1] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  float zvtx = NAN;
  float maxJetET = 0;
  float maxJetPhi = NAN;
  float subJetET = 0;
  float subJetPhi = NAN;
  float frcem = 0;
  float frcoh = 0;
  float dPhi = NAN;

  if (!towersEM || !towersOH || !geom[0] || !geom[1] || !gvtxmap)
  {
    if (_debug > 0 && !_missingInfoWarningPrinted)
    {
      std::cout << "Missing critical info; abort event. Further warnings will be suppressed. AddressOf towersEM/towersOH/geomIH/geomOH/gvtxmap : " << towersEM << "/" << towersOH << "/" << geom[0] << "/" << geom[1] << "/" << gvtxmap << std::endl;
    }
    _missingInfoWarningPrinted = true;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (gvtxmap)
  {
    if (gvtxmap->empty())
    {
      if (_debug > 0)
      {
        std::cout << "gvtxmap empty - aborting event." << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    GlobalVertex *gvtx = gvtxmap->begin()->second;
    if (gvtx)
    {
      auto startIter = gvtx->find_vertexes(_vtxtype);
      auto endIter = gvtx->end_vertexes();
      for (auto iter = startIter; iter != endIter; ++iter)
      {
        const auto &[type, vertexVec] = *iter;
        if (type != _vtxtype)
        {
          continue;
        }
        for (const auto *vertex : vertexVec)
        {
          if (!vertex)
          {
            continue;
          }
          zvtx = vertex->get_z();
        }
      }
    }
    else
    {
      if (_debug > 0)
      {
        std::cout << "gvtx is NULL! Aborting event." << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if (std::isnan(zvtx))
  {
    if (_debug > 0)
    {
      std::cout << "zvtx is NAN after attempting to grab it. ABORT EVENT!" << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (_debug > 1)
  {
    std::cout << "Getting jets: " << std::endl;
  }

  if (jets)
  {
    int tocheck = jets->size();
    if (_debug > 2)
    {
      std::cout << "Found " << tocheck << " jets to check..." << std::endl;
    }
    for (int i = 0; i < tocheck; ++i)
    {
      float jetET = 0;
      float jetPhi = NAN;
      float jetEta = NAN;
      Jet *jet = jets->get_jet(i);
      if (jet)
      {
        jetEta = jet->get_eta();
        jetET = jet->get_e() / std::cosh(jetEta);
        jetPhi = jet->get_phi();
      }
      else
      {
        continue;
      }
      if (jetET < 8)
      {
        continue;
      }
      if (_debug > 2)
      {
        std::cout << "found a good jet!" << std::endl;
      }
      if (jetET > maxJetET)
      {
        if (maxJetET)
        {
          subJetET = maxJetET;
          subJetPhi = maxJetPhi;
        }
        maxJetET = jetET;
        maxJetPhi = jetPhi;
      }
      else if (jetET > subJetET)
      {
        subJetET = jetET;
        subJetPhi = jetPhi;
        continue;
      }
      else
      {
        continue;
      }
      frcem = 0;
      frcoh = 0;
      for (auto comp : jet->get_comp_vec())
      {
        unsigned int channel = comp.second;
        TowerInfo *tower;
        if (comp.first == 13 || comp.first == 28 || comp.first == 25)
        {
          tower = towersEM->get_tower_at_channel(channel);
          int key = towersEM->encode_key(channel);
          const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towersEM->getTowerEtaBin(key), towersEM->getTowerPhiBin(key));
          RawTowerGeom *tower_geom = geom[0]->get_tower_geometry(geomkey);
          float radius = 93.5;
          float ihEta = tower_geom->get_eta();
          float emZ = radius / (std::tan(2 * std::atan(std::exp(-ihEta))));
          float newz = emZ - zvtx;
          float newTheta = std::atan2(radius, newz);
          float towerEta = -log(std::tan(0.5 * newTheta));
          frcem += tower->get_energy() / std::cosh(towerEta);
        }
        if (comp.first == 7 || comp.first == 27)
        {
          tower = towersOH->get_tower_at_channel(channel);
          int key = towersOH->encode_key(channel);
          const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, towersOH->getTowerEtaBin(key), towersOH->getTowerPhiBin(key));
          RawTowerGeom *tower_geom = geom[1]->get_tower_geometry(geomkey);
          float radius = tower_geom->get_center_radius();
          float newz = tower_geom->get_center_z() - zvtx;
          float newTheta = std::atan2(radius, newz);
          float towerEta = -log(std::tan(0.5 * newTheta));
          frcoh += tower->get_energy() / std::cosh(towerEta);
        }
      }

      frcem /= maxJetET;
      frcoh /= maxJetET;
    }
  }
  else
  {
    if (_debug > 0)
    {
      std::cout << "No jet node!" << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  bool isDijet = false;
  if (subJetET > 8)
  {
    isDijet = true;
    dPhi = std::abs(maxJetPhi - subJetPhi);
    if (dPhi > M_PI)
    {
      dPhi = 2 * M_PI - dPhi;
    }
  }
  bool dPhiCut = failsdPhiCut(dPhi, isDijet);

  bool failsLoEm = failsLoEmFracETCut(frcem, maxJetET, dPhiCut, isDijet);
  bool failsHiEm = failsHiEmFracETCut(frcem, maxJetET, dPhiCut, isDijet);
  bool failsIhCut = failsIhFracCut(frcem, frcoh);

  bool failsAnyCut = failsLoEm || failsHiEm || failsIhCut;

  if (failsAnyCut && _doAbort)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  PHNodeIterator iter(topNode);
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  _cutParams.set_int_param("failsLoEmJetCut", failsLoEm);
  _cutParams.set_int_param("failsHiEmJetCut", failsHiEm);
  _cutParams.set_int_param("failsIhJetCut", failsIhCut);
  _cutParams.set_int_param("failsAnyJetCut", failsAnyCut);
  _cutParams.set_int_param("isDijet", isDijet);
  _cutParams.set_double_param("frcem", frcem);
  _cutParams.set_double_param("frcoh", frcoh);
  _cutParams.set_double_param("maxJetET", maxJetET);
  _cutParams.set_double_param("dPhi", dPhi);
  _cutParams.UpdateNodeTree(parNode, "JetCutParams");

  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________________________________..
int JetBackgroundCut::ResetEvent(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "JetBackgroundCut::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetBackgroundCut::End(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "JetBackgroundCut::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetBackgroundCut::Reset(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "JetBackgroundCut::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void JetBackgroundCut::Print(const std::string &what) const
{
  std::cout << "JetBackgroundCut::Print(const std::string &what) const Printing info for " << what << std::endl;
}
