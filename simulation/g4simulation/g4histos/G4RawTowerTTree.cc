#include "G4RawTowerTTree.h"

#include "G4RootRawTower.h"
#include "G4RootRawTowerContainer.h"

#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <TH1.h>
#include <TSystem.h>

#include <iostream>  // for operator<<, endl, basic_...

G4RawTowerTTree::G4RawTowerTTree(const std::string &name)
  : SubsysReco(name)
{
}

int G4RawTowerTTree::Init(PHCompositeNode *topNode)
{
  if (_detector.empty())
  {
    std::cout << "Detector not set via Detector(<name>) method" << std::endl;
    std::cout << "(it is the name appended to the G4TOWER_<name> nodename)" << std::endl;
    std::cout << "you do not want to run like this, exiting now" << std::endl;
    gSystem->Exit(1);
  }
  hm = new Fun4AllHistoManager("TOWERHIST");
  etot_hist = new TH1F("etot", "total deposited energy", 200, 0, 20);
  hm->registerHisto(etot_hist);
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  G4RootRawTowerContainer *towers = new G4RootRawTowerContainer();
  PHIODataNode<PHObject> *node = new PHIODataNode<PHObject>(towers, _outnodename, "PHObject");
  dstNode->addNode(node);
  evtno = 0;
  return 0;
}

int G4RawTowerTTree::process_event(PHCompositeNode *topNode)
{
  evtno++;
  G4RootRawTowerContainer *towers = findNode::getClass<G4RootRawTowerContainer>(topNode, _outnodename);
  RawTowerGeomContainer *rawtowergeom = findNode::getClass<RawTowerGeomContainer>(topNode, _towergeomnodename);

  // RawTowerContainer *g4towers = findNode::getClass<RawTowerContainer>(topNode, _towernodename);
  TowerInfoContainer *g4towers = findNode::getClass<TowerInfoContainer>(topNode, _towernodename);
  if (!g4towers)
  {
    std::cout << "could not find " << _towernodename << std::endl;
    gSystem->Exit(1);
  }

  double etot = 0;

  unsigned int nchannels = g4towers->size();
  for (unsigned int channel = 0; channel < nchannels; channel++)
  {
    TowerInfo *intower = g4towers->get_tower_at_channel(channel);
    if (savetowers)
    {
      unsigned int towerkey = g4towers->encode_key(channel);
      int ieta = g4towers->getTowerEtaBin(towerkey);
      int iphi = g4towers->getTowerPhiBin(towerkey);

      G4RootRawTower roottwr(rawtowergeom->get_etacenter(ieta), rawtowergeom->get_phicenter(iphi), intower->get_energy());
      towers->AddG4RootRawTower(roottwr);
    }
    etot += intower->get_energy();
  }
  etot_hist->Fill(etot);
  towers->set_etotal(etot);
  towers->set_event(evtno);
  return 0;
}

int G4RawTowerTTree::End(PHCompositeNode * /*topNode*/)
{
  hm->dumpHistos(_histofilename);
  delete hm;
  return 0;
}

void G4RawTowerTTree::Detector(const std::string &det)
{
  _detector = det;
  _outnodename = "G4RootRawTower_" + det;
  _towernodename = "TOWERINFO_CALIB_" + det;
  _towergeomnodename = "TOWERGEOM_" + det;
  if (_histofilename.empty())
  {
    _histofilename = "RawTowerHistos_" + det + ".root";
  }
}
