#include "G4TowerNtuple.h"

#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainerv1.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>

#include <TFile.h>
#include <TH1.h>
#include <TNtuple.h>

#include <cmath>     // for isfinite
#include <iostream>  // for operator<<, basic_ostream
#include <utility>   // for pair, make_pair

G4TowerNtuple::G4TowerNtuple(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , m_filename(filename)
{
}

G4TowerNtuple::~G4TowerNtuple()
{
  delete hm;
}

int G4TowerNtuple::Init(PHCompositeNode * /*unused*/)
{
  delete hm; // make cppcheck happy
  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(m_filename.c_str(), "RECREATE");
  ntup = new TNtuple("towerntup", "G4Towers", "detid:phi:eta:energy");
  //  ntup->SetDirectory(0);
  TH1 *h1 = new TH1F("energy1GeV", "energy 0-1GeV", 1000, 0, 1);
  eloss.push_back(h1);
  h1 = new TH1F("energy100GeV", "energy 0-100GeV", 1000, 0, 100);
  eloss.push_back(h1);
  return 0;
}

int G4TowerNtuple::process_event(PHCompositeNode *topNode)
{
  std::string nodename;
  std::string geonodename;
  for (const auto &iter : m_node_postfix)
  {
    int detid = (m_detid.find(iter))->second;
    nodename = "TOWERINFO_" + m_tower_type[iter];
    geonodename = "TOWERGEOM_" + iter;
    RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, geonodename);
    if (!towergeom)
    {
      std::cout << "no geometry node " << geonodename << " for " << iter << std::endl;
      continue;
    }

    TowerInfoContainer *towers = findNode::getClass<TowerInfoContainer>(topNode, nodename);
    if (towers)
    {
      double esum = 0;
      unsigned int nchannels = towers->size();
      for (unsigned int channel = 0; channel < nchannels; channel++)
      {
        TowerInfo *tower = towers->get_tower_at_channel(channel);
        double energy = tower->get_energy();
        if (!std::isfinite(energy))
        {
          std::cout << "invalid energy: " << energy << std::endl;
        }
        esum += energy;
        unsigned int towerkey = towers->encode_key(channel);
        int etabin = towers->getTowerEtaBin(towerkey);
        int phibin = towers->getTowerPhiBin(towerkey);

        // to search the map fewer times, cache the geom object until the layer changes
        double phi = towergeom->get_phicenter(phibin);
        double eta = towergeom->get_etacenter(etabin);
        ntup->Fill(detid,
                   phi,
                   eta,
                   energy);
      }
      for (auto *eiter : eloss)
      {
        eiter->Fill(esum);
      }
    }
  }
  return 0;
}

int G4TowerNtuple::End(PHCompositeNode * /*topNode*/)
{
  outfile->cd();
  ntup->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(m_filename, "UPDATE");
  return 0;
}

void G4TowerNtuple::AddNode(const std::string &name, const std::string &twrtype, const int detid)
{
  std::string twrname;
  twrname = twrtype + "_" + name;
  m_node_postfix.insert(name);
  m_tower_type.insert(make_pair(name, twrname));
  m_detid[name] = detid;
  return;
}
