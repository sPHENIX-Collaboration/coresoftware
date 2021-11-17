#include "G4HitNtuple.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>           // for SubsysReco

#include <phool/getClass.h>

#include <TFile.h>
#include <TH1.h>
#include <TNtuple.h>

#include <sstream>
#include <utility>                        // for pair

using namespace std;

G4HitNtuple::G4HitNtuple(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , nblocks(0)
  , hm(nullptr)
  , _filename(filename)
  , ntup(nullptr)
  , outfile(nullptr)
{
}

G4HitNtuple::~G4HitNtuple()
{
  //  delete ntup;
  delete hm;
}

int G4HitNtuple::Init(PHCompositeNode *)
{
  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(_filename.c_str(), "RECREATE");
  ntup = new TNtuple("hitntup", "G4Hits", "detid:lyr:slat:x0:y0:z0:x1:y1:z1:edep");
  //  ntup->SetDirectory(0);
  TH1 *h1 = new TH1F("edep1GeV", "edep 0-1GeV", 1000, 0, 1);
  eloss.push_back(h1);
  h1 = new TH1F("edep100GeV", "edep 0-100GeV", 1000, 0, 100);
  eloss.push_back(h1);
  return 0;
}

int G4HitNtuple::process_event(PHCompositeNode *topNode)
{
  ostringstream nodename;
  set<string>::const_iterator iter;
  vector<TH1 *>::const_iterator eiter;
  for (iter = _node_postfix.begin(); iter != _node_postfix.end(); ++iter)
  {
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
        esum += hit_iter->second->get_edep();
        ntup->Fill(detid,
                   hit_iter->second->get_layer(),
                   hit_iter->second->get_scint_id(),
                   hit_iter->second->get_x(0),
                   hit_iter->second->get_y(0),
                   hit_iter->second->get_z(0),
                   hit_iter->second->get_x(1),
                   hit_iter->second->get_y(1),
                   hit_iter->second->get_z(1),
                   hit_iter->second->get_edep());
      }
      for (eiter = eloss.begin(); eiter != eloss.end(); ++eiter)
      {
        (*eiter)->Fill(esum);
      }
    }
  }
  return 0;
}

int G4HitNtuple::End(PHCompositeNode */*topNode*/)
{
  outfile->cd();
  ntup->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(_filename, "UPDATE");
  return 0;
}

void G4HitNtuple::AddNode(const std::string &name, const int detid)
{
  _node_postfix.insert(name);
  _detid[name] = detid;
  return;
}
