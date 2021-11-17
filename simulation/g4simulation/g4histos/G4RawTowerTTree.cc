#include "G4RawTowerTTree.h"

#include "G4RootRawTower.h"
#include "G4RootRawTowerContainer.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>              // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>              // for PHIODataNode
#include <phool/PHNodeIterator.h>            // for PHNodeIterator
#include <phool/PHObject.h>                  // for PHObject
#include <phool/getClass.h>

#include <TH1.h>
#include <TSystem.h>

#include <iostream>                          // for operator<<, endl, basic_...
#include <map>                               // for _Rb_tree_const_iterator
#include <utility>                           // for pair

using namespace std;

G4RawTowerTTree::G4RawTowerTTree(const std::string &name)
  : SubsysReco(name)
  , savetowers(1)
  , evtno(0)
  , hm(nullptr)
  , etot_hist(nullptr)
{
}

int G4RawTowerTTree::Init(PHCompositeNode *topNode)
{
  if (_detector.empty())
  {
    cout << "Detector not set via Detector(<name>) method" << endl;
    cout << "(it is the name appended to the G4TOWER_<name> nodename)" << endl;
    cout << "you do not want to run like this, exiting now" << endl;
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

  RawTowerContainer *g4towers = findNode::getClass<RawTowerContainer>(topNode, _towernodename);
  if (!g4towers)
  {
    cout << "could not find " << _towernodename << endl;
    gSystem->Exit(1);
  }
  RawTowerContainer::ConstRange tower_range = g4towers->getTowers();

  double etot = 0;
  for (RawTowerContainer::ConstIterator tower_iter = tower_range.first; tower_iter != tower_range.second; tower_iter++)
  {
    RawTower *intower = tower_iter->second;
    if (savetowers)
    {
      G4RootRawTower roottwr(rawtowergeom->get_etacenter(intower->get_bineta()), rawtowergeom->get_phicenter(intower->get_binphi()), intower->get_energy());
      towers->AddG4RootRawTower(roottwr);
    }
    etot += intower->get_energy();
  }
  etot_hist->Fill(etot);
  towers->set_etotal(etot);
  towers->set_event(evtno);
  return 0;
}

int G4RawTowerTTree::End(PHCompositeNode */*topNode*/)
{
  hm->dumpHistos(_histofilename);
  delete hm;
  return 0;
}

void G4RawTowerTTree::Detector(const std::string &det)
{
  _detector = det;
  _outnodename = "G4RootRawTower_" + det;
  _towernodename = "TOWER_CALIB_" + det;
  _towergeomnodename = "TOWERGEOM_" + det;
  if (!_histofilename.size())
  {
    _histofilename = "RawTowerHistos_" + det + ".root";
  }
}
