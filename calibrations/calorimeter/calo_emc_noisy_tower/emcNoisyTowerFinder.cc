#include "emcNoisyTowerFinder.h"

// Tower stuff
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

// Fun4All
#include <ffaobjects/EventHeader.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

// ROOT
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLine.h>
#include <TPad.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <boost/format.hpp>

//________________________________
emcNoisyTowerFinder::emcNoisyTowerFinder(const std::string &name, const std::string &outputName)
  : SubsysReco(name)
  , Outfile(outputName)
{
  // initialize tree with SG vs Kurary fiber information

  std::cout << "emcNoisyTowerFinder::emcNoisyTowerFinder Calling ctor" << std::endl;
}
//__________________________________
emcNoisyTowerFinder::~emcNoisyTowerFinder()
{
  std::cout << "emcNoisyTowerFinder::~emcNoisyTowerFinder() Calling dtor" << std::endl;
}

//_____________________________
int emcNoisyTowerFinder::Init(PHCompositeNode * /*topNode*/)
{
  std::cout << "emcNoisyTowerFinder::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  foutput = new TFile(Outfile.c_str(), "recreate");

  h_hits_eta_phi_adc = new TH2F("h_hits_eta_phi_adc", "", Neta, 0, Neta, Nphi, 0, Nphi);
  pr_hits_eta_phi_adc = new TProfile2D("pr_hits_eta_phi_adc", "", Neta, 0, Neta, Nphi, 0, Nphi);

  h_hits_eta_phi_gev = new TH2F("h_hits_eta_phi_gev", "", Neta, 0, Neta, Nphi, 0, Nphi);
  pr_hits_eta_phi_gev = new TProfile2D("pr_hits_eta_phi_gev", "", Neta, 0, Neta, Nphi, 0, Nphi);

  std::string default_time_independent_calib = "cemc_pi0_twrSlope_v1_default";
  m_fieldname = "Femc_datadriven_qm1_correction";

  std::string calibdir = CDBInterface::instance()->getUrl(default_time_independent_calib);
  if (calibdir.empty())
  {
    std::cout << "CaloTowerCalib::::InitRun No EMCal Calibration NOT even a default" << std::endl;
    exit(1);
  }
  cdbttree = new CDBTTree(calibdir);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________
int emcNoisyTowerFinder::InitRun(PHCompositeNode * /*topNode*/)
{
  std::cout << "emcNoisyTowerFinder::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________
int emcNoisyTowerFinder::process_event(PHCompositeNode *topNode)
{
  // Get TowerInfoContainer
  TowerInfoContainer *towers;
  towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_CEMC");
  if (!towers)
  {
    std::cout << PHWHERE << "emcNoisyTowerFinder::process_event Could not find node TOWERS_CEMC" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  // iterate through all towers, incrementing their Frequency arrays if they record a hit

  int tower_range = towers->size();
  for (int j = 0; j < tower_range; j++)
  {
    TowerInfo *tower = towers->get_tower_at_channel(j);
    float energy = tower->get_energy();
    unsigned int towerkey = towers->encode_key(j);
    int ieta = towers->getTowerEtaBin(towerkey);
    int iphi = towers->getTowerPhiBin(towerkey);

    float calibconst = cdbttree->GetFloatValue(towerkey, m_fieldname);
    float calib_energy = calibconst * energy;

    if (energy > energy_threshold_adc)
    {
      h_hits_eta_phi_adc->Fill(ieta + 1, iphi + 1);
      pr_hits_eta_phi_adc->Fill(ieta + 1, iphi + 1, 1);
    }
    else
    {
      pr_hits_eta_phi_adc->Fill(ieta + 1, iphi + 1, 0);
    }

    if (calib_energy > energy_threshold_gev)
    {
      h_hits_eta_phi_gev->Fill(ieta + 1, iphi + 1);
      pr_hits_eta_phi_gev->Fill(ieta + 1, iphi + 1, 1);
    }
    else
    {
      pr_hits_eta_phi_gev->Fill(ieta + 1, iphi + 1, 0);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//__________________________
int emcNoisyTowerFinder::ResetEvent(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
//__________________________

//__________________________
void emcNoisyTowerFinder::FindHot(std::string &infilename, std::string &outfilename, const std::string &inHist)
{
  TFile *fin = new TFile(infilename.c_str());
  if (!fin)
  {
    std::cout << "emcNoisyTowerFinder::FindHot: input file not found " << infilename.c_str() << std::endl;
    return;
  }
  h_hits_eta_phi_adc = (TH2F *) fin->Get(inHist.c_str());
    if (!h_hits_eta_phi_adc)
  {
    std::cout << "emcNoisyTowerFinder::FindHot: input hist not found " << inHist.c_str() << std::endl;
    return;
  }

  TH2F *h_hits = h_hits_eta_phi_adc;

  TFile *fout = new TFile(outfilename.c_str(), "recreate");

  TH2F *h_hot = new TH2F("h_hot", "", Neta, 0, Neta, Nphi, 0, Nphi);
  TH1F *h1_hits[Neta];
  TH1F *h1_hits2[Neta];

  float max = h_hits->GetBinContent(h_hits->GetMaximumBin());
  float min = h_hits->GetMinimum();
  if (min == 0)
  {
    min = 1;
  }

  for (int ie = 0; ie < Neta; ie++)
  {
    h1_hits[ie] = new TH1F((std::string("h1_hits") + std::to_string(ie)).c_str(), "", 200, min, max);
    h1_hits2[ie] = new TH1F((std::string("h1_hits2") + std::to_string(ie)).c_str(), "", 200, min, max);
    for (int iphi = 0; iphi < 256; iphi++)
    {
      h1_hits[ie]->Fill(h_hits->GetBinContent(ie + 1, iphi + 1));
    }

    float mean = h1_hits[ie]->GetMean();
    float std = h1_hits[ie]->GetStdDev();
    for (int iphi = 0; iphi < Nphi; iphi++)
    {
      float val = h_hits->GetBinContent(ie + 1, iphi + 1);
      if (std::fabs(mean - val) / std > 3)
      {
        continue;
      }
      h1_hits2[ie]->Fill(val);
    }
    mean = h1_hits[ie]->GetMean();
    std = h1_hits[ie]->GetStdDev();
    for (int iphi = 0; iphi < Nphi; iphi++)
    {
      h_hot->SetBinContent(ie + 1, iphi + 1, 0);
      float val = h_hits->GetBinContent(ie + 1, iphi + 1);

      if ((val - mean) / std > sigma_bad_thresh)  // hot tower
      {
        h_hot->SetBinContent(ie + 1, iphi + 1, 2);
        h_hits->SetBinContent(ie + 1, iphi + 1, 0);
      }
      if ((val - mean) / std < -1 * sigma_bad_thresh || val < mean * percent_cold_thresh)  // cold tower
      {
        h_hot->SetBinContent(ie + 1, iphi + 1, 3);
        h_hits->SetBinContent(ie + 1, iphi + 1, 0);
      }
      if (val == 0)  // dead tower
      {
        h_hot->SetBinContent(ie + 1, iphi + 1, 1);
      }
    }
  }
  fout->Write();
  fout->Close();
}

//____________________________________________________________________________..
int emcNoisyTowerFinder::Reset(PHCompositeNode * /*topNode*/)
{
  std::cout << "emcNoisyTowerFinder::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________
void emcNoisyTowerFinder::Print(const std::string &what) const
{
  std::cout << "emcNoisyTowerFinder::Print(const std::string &what) const Printing info for " << what << std::endl;
}
//____________________________________________________
int emcNoisyTowerFinder::End(PHCompositeNode * /*topNode*/)
{
  foutput->Write();
  // fchannels -> cd();
  // fchannels->Close();
  // delete fchannels;
  // fchannels=NULL;

  std::cout << "emcNoisyTowerFinder::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________
int emcNoisyTowerFinder::EndRun(int runnumber)
{
  std::cout << "emcNoisyTowerFinder::EndRun: this is the end of run: " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________
