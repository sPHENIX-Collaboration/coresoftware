#include "G4ScintillatorTowerTTree.h"
#include "G4RootScintillatorTowerContainer.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainerv1.h>


#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>                // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>                // for PHIODataNode
#include <phool/PHNode.h>                      // for PHNode
#include <phool/PHNodeIterator.h>              // for PHNodeIterator
#include <phool/PHObject.h>                    // for PHObject
#include <phool/getClass.h>

#include <TH1.h>
#include <TSystem.h>

#include <iostream>                            // for operator<<, endl, basi...
#include <map>                                 // for _Rb_tree_const_iterator
#include <utility>                             // for pair

using namespace std;

G4ScintillatorTowerTTree::G4ScintillatorTowerTTree(const std::string &name)
  : SubsysReco(name)
  , savetowers(1)
  , evtno(0)
  , hm(nullptr)
  , etot_hist(nullptr)
{
}

int G4ScintillatorTowerTTree::Init(PHCompositeNode *topNode)
{
  if (!_detector.size())
  {
    cout << "Detector not set via Detector(<name>) method" << endl;
    cout << "(it is the name appended to the G4TOWER_<name> nodename)" << endl;
    cout << "you do not want to run like this, exiting now" << endl;
    gSystem->Exit(1);
  }
  hm = new Fun4AllHistoManager("SCINTILLATORTOWERHIST");
  etot_hist = new TH1F("etot", "total deposited energy", 200, 0, 20);
  hm->registerHisto(etot_hist);
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  G4RootScintillatorTowerContainer *towers = new G4RootScintillatorTowerContainer();
  PHIODataNode<PHObject> *node = new PHIODataNode<PHObject>(towers, _outnodename, "PHObject");
  dstNode->addNode(node);
  evtno = 0;
  return 0;
}

int G4ScintillatorTowerTTree::process_event(PHCompositeNode *topNode)
{
  evtno++;
  G4RootScintillatorTowerContainer *towers = findNode::getClass<G4RootScintillatorTowerContainer>(topNode, _outnodename);

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
	  double towerenergy = intower->get_energy();
	  towers->AddTower(towerenergy,ieta,iphi);
	}
      etot += intower->get_energy();
    }
  etot_hist->Fill(etot);
  towers->set_etotal(etot);
  towers->set_event(evtno);
  return 0;
}

int G4ScintillatorTowerTTree::End(PHCompositeNode */*topNode*/)
{
  hm->dumpHistos(_histofilename);
  delete hm;
  return 0;
}

void G4ScintillatorTowerTTree::Detector(const std::string &det)
{
  _detector = det;
  _outnodename = "G4RootScintillatorTower_" + det;
  _towernodename = "TOWERINFO_" + det;
  if (!_histofilename.size())
  {
    _histofilename = "ScintillatorTowerHistos_" + det + ".root";
  }
}
