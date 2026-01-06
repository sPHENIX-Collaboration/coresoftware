#include "G4ScintillatorSlatTTree.h"

#include "G4RootScintillatorSlatContainer.h"

#include <g4detectors/PHG4ScintillatorSlat.h>
#include <g4detectors/PHG4ScintillatorSlatContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <TH1.h>
#include <TSystem.h>

#include <iostream>  // for operator<<, endl
#include <map>       // for _Rb_tree_cons...
#include <utility>   // for pair

G4ScintillatorSlatTTree::G4ScintillatorSlatTTree(const std::string &name)
  : SubsysReco(name)
{
}

int G4ScintillatorSlatTTree::Init(PHCompositeNode *topNode)
{
  if (_detector.empty())
  {
    std::cout << "Detector not set via Detector(<name>) method" << std::endl;
    std::cout << "(it is the name appended to the G4CELL_<name> nodename)" << std::endl;
    std::cout << "you do not want to run like this, exiting now" << std::endl;
    gSystem->Exit(1);
  }
  hm = new Fun4AllHistoManager("SCINTILLATORSLATHIST");
  etot_hist = new TH1F("etot", "total deposited energy", 200, 0, 20);
  hm->registerHisto(etot_hist);
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  G4RootScintillatorSlatContainer *slats = new G4RootScintillatorSlatContainer();
  PHIODataNode<PHObject> *node = new PHIODataNode<PHObject>(slats, _outnodename, "PHObject");
  dstNode->addNode(node);
  evtno = 0;
  return 0;
}

int G4ScintillatorSlatTTree::process_event(PHCompositeNode *topNode)
{
  evtno++;
  G4RootScintillatorSlatContainer *slats = findNode::getClass<G4RootScintillatorSlatContainer>(topNode, _outnodename);

  PHG4ScintillatorSlatContainer *g4slats = findNode::getClass<PHG4ScintillatorSlatContainer>(topNode, _slatnodename);
  if (!g4slats)
  {
    std::cout << "could not find " << _slatnodename << std::endl;
    gSystem->Exit(1);
  }
  PHG4ScintillatorSlatContainer::ConstRange slat_range = g4slats->getScintillatorSlats();

  double etot = 0;
  for (PHG4ScintillatorSlatContainer::ConstIterator slat_iter = slat_range.first; slat_iter != slat_range.second; slat_iter++)
  {
    PHG4ScintillatorSlat *inslat = slat_iter->second;
    if (saveslats)
    {
      slats->AddSlat(*inslat);
    }
    etot += inslat->get_edep();
  }
  etot_hist->Fill(etot);
  slats->set_etotal(etot);
  slats->set_event(evtno);
  return 0;
}

int G4ScintillatorSlatTTree::End(PHCompositeNode * /*topNode*/)
{
  hm->dumpHistos(_histofilename);
  delete hm;
  return 0;
}

void G4ScintillatorSlatTTree::Detector(const std::string &det)
{
  _detector = det;
  _outnodename = "G4RootScintillatorSlat_" + det;
  _slatnodename = "G4CELL_" + det;
  if (_histofilename.empty())
  {
    _histofilename = "ScintillatorSlatHistos_" + det + ".root";
  }
}
