#include "G4EdepNtuple.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>           // for SubsysReco

#include <phool/getClass.h>

#include <TFile.h>
#include <TNtuple.h>

#include <sstream>
#include <utility>                        // for pair

using namespace std;

G4EdepNtuple::G4EdepNtuple(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , nblocks(0)
  , hm(nullptr)
  , _filename(filename)
  , ntup(nullptr)
  , outfile(nullptr)
{
}

G4EdepNtuple::~G4EdepNtuple()
{
  delete hm;
}

int G4EdepNtuple::Init(PHCompositeNode *)
{
  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(_filename.c_str(), "RECREATE");
  ntup = new TNtuple("edepntup", "G4Edeps", "detid:layer:edep");
  return 0;
}

int G4EdepNtuple::process_event(PHCompositeNode *topNode)
{
  ostringstream nodename;
  set<string>::const_iterator iter;
  map<int, double> layer_edep_map;
  map<int, double>::const_iterator edepiter;
  for (iter = _node_postfix.begin(); iter != _node_postfix.end(); ++iter)
  {
    layer_edep_map.clear();
    int detid = (_detid.find(*iter))->second;
    nodename.str("");
    nodename << "G4HIT_" << *iter;
    PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
    if (hits)
    {
      double esum = 0;
      //          double numhits = hits->size();
      //          nhits[i]->Fill(numhits);
      PHG4HitContainer::ConstRange hit_range = hits->getHits();
      for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)

      {
        layer_edep_map[hit_iter->second->get_layer()] += hit_iter->second->get_edep();
        esum += hit_iter->second->get_edep();
      }
      for (edepiter = layer_edep_map.begin(); edepiter != layer_edep_map.end(); ++edepiter)
      {
        ntup->Fill(detid, edepiter->first, edepiter->second);
      }
      ntup->Fill(detid, -1, esum);  // fill sum over all layers for each detector
    }
  }
  return 0;
}

int G4EdepNtuple::End(PHCompositeNode */*topNode*/)
{
  outfile->cd();
  ntup->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(_filename, "UPDATE");
  return 0;
}

void G4EdepNtuple::AddNode(const std::string &name, const int detid)
{
  _node_postfix.insert(name);
  _detid[name] = detid;
  return;
}
