#include "G4RawTowerTTree.h"

#include "G4RootRawTower.h"
#include "G4RootRawTowerContainer.h"

#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainerv1.h>


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

  // RawTowerContainer *g4towers = findNode::getClass<RawTowerContainer>(topNode, _towernodename);
  TowerInfoContainer *g4towers = findNode::getClass<TowerInfoContainerv1>(topNode, _towernodename);
  if (!g4towers)
  {
    cout << "could not find " << _towernodename << endl;
    gSystem->Exit(1);
  }

  double etot = 0;
  
  unsigned int nchannels = g4towers->size();
  for (unsigned int channel = 0; channel < nchannels;channel++)
    {
      TowerInfo *intower =g4towers->get_tower_at_channel(channel);
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
  _towernodename = "TOWERINFO_CALIB_" + det;
  _towergeomnodename = "TOWERGEOM_" + det;
  if (!_histofilename.size())
  {
    _histofilename = "RawTowerHistos_" + det + ".root";
  }
}
