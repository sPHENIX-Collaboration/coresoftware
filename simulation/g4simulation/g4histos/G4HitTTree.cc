#include "G4HitTTree.h"

#include "G4RootHitContainer.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>           // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>           // for PHIODataNode
#include <phool/PHNodeIterator.h>         // for PHNodeIterator
#include <phool/PHObject.h>               // for PHObject
#include <phool/getClass.h>

#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>

#include <iostream>                       // for operator<<, endl, basic_ost...
#include <map>                            // for _Rb_tree_const_iterator
#include <utility>                        // for pair

using namespace std;

G4HitTTree::G4HitTTree(const std::string &name)
  : SubsysReco(name)
  , savehits(1)
  , evtno(0)
  , hm(nullptr)
  , etot_hist(nullptr)
  , eion_etot_hist(nullptr)
{
  BlackHoleName("BH_1");  // initialize this to what we have in our common sims
}

int G4HitTTree::Init(PHCompositeNode *topNode)
{
  if (!_detector.size())
  {
    cout << "Detector not set via Detector(<name>) method" << endl;
    cout << "(it is the name appended to the G4HIT_<name> nodename)" << endl;
    cout << "you do not want to run like this, exiting now" << endl;
    gSystem->Exit(1);
  }
  hm = new Fun4AllHistoManager("HITHIST");
  etot_hist = new TH1F("etot", "total deposited energy", 200, 0, 20);
  hm->registerHisto(etot_hist);
  eion_etot_hist = new TH2F("eion_etot", "ionization energy vs total energy", 200, 0, 20, 200, 0, 20);
  hm->registerHisto(eion_etot_hist);

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  G4RootHitContainer *hits = new G4RootHitContainer();
  PHIODataNode<PHObject> *node = new PHIODataNode<PHObject>(hits, _outnodename, "PHObject");
  dstNode->addNode(node);
  evtno = 0;
  return 0;
}

int G4HitTTree::process_event(PHCompositeNode *topNode)
{
  evtno++;
  G4RootHitContainer *hits = findNode::getClass<G4RootHitContainer>(topNode, _outnodename);
  PHG4HitContainer *g4hits = findNode::getClass<PHG4HitContainer>(topNode, _hitnodename);
  double etot = 0;
  double eion = 0;
  if (g4hits)
  {
    PHG4HitContainer::ConstRange hit_range = g4hits->getHits();
    //shower_z->Reset();
    //  cout << "Number of Hits: " << g4hits->size() << endl;
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
    {
      PHG4Hit *inhit = hit_iter->second;
      if (savehits)
      {
        hits->AddHit(inhit);
      }
      etot += inhit->get_edep();
      eion += inhit->get_eion();
    }
    etot_hist->Fill(etot);
    eion_etot_hist->Fill(etot, eion);
  }
  g4hits = findNode::getClass<PHG4HitContainer>(topNode, _absorbernodename);
  if (g4hits)
  {
    PHG4HitContainer::ConstRange hit_range = g4hits->getHits();
    //shower_z->Reset();
    //  cout << "Number of Hits: " << g4hits->size() << endl;
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
    {
      PHG4Hit *inhit = hit_iter->second;
      if (savehits)
      {
        hits->AddHit(inhit);
      }
    }
  }
  double eleak = 0;
  g4hits = findNode::getClass<PHG4HitContainer>(topNode, _blackholenodename);
  if (g4hits)
  {
    PHG4HitContainer::ConstRange hit_range = g4hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
    {
      if (savehits)
      {
        //	      PHG4Hit *g4h = hits->AddHit( *(hit_iter->second));
      }
      eleak += hit_iter->second->get_edep();
    }
  }
  hits->set_etotal(etot);
  hits->set_eion(eion);
  hits->set_leakage(eleak);
  hits->set_event(evtno);
  return 0;
}

int G4HitTTree::End(PHCompositeNode */*topNode*/)
{
  hm->dumpHistos("HitHistos.root");
  delete hm;
  return 0;
}

void G4HitTTree::Detector(const std::string &det)
{
  _detector = det;
  _outnodename = "G4RootHit_" + det;
  _hitnodename = "G4HIT_" + det;
  _absorbernodename = "G4HIT_ABSORBER_" + det;
}

void G4HitTTree::BlackHoleName(const std::string &bh)
{
  _blackholenodename = "G4HIT_" + bh;
}
