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

#include <fstream>
#include <iostream>
#include <sstream>


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
  //TH2F *h_hits_eta_phi_adc = nullptr;
  bool isListFile = (infilename.substr(infilename.rfind('.') + 1) == "txt" || infilename.substr(infilename.rfind('.') + 1) == "list");

  if (isListFile)
  {
    std::ifstream fileList(infilename.c_str());
    std::string line;
    while (std::getline(fileList, line))
    {
      TFile *fin = new TFile(line.c_str());
      if (!fin || fin->IsZombie())
      {
        std::cout << "emcNoisyTowerFinder::FindHot: input file not found or is corrupted " << line.c_str() << std::endl;
        continue;
      }
      TH2F *tempHist = (TH2F *) fin->Get(inHist.c_str());
      if (!tempHist)
      {
        std::cout << "emcNoisyTowerFinder::FindHot: input hist not found in file " << line.c_str() << std::endl;
        delete fin;
        continue;
      }
      if (!h_hits_eta_phi_adc)
      {
        h_hits_eta_phi_adc = (TH2F *) tempHist->Clone();
        h_hits_eta_phi_adc->SetDirectory(nullptr);  // Detach from the file to keep it in memory
      }
      else
      {
        h_hits_eta_phi_adc->Add(tempHist);
      }
      delete fin;
    }
    fileList.close();
  }
  else
  {
    TFile *fin = new TFile(infilename.c_str());
    if (!fin || fin->IsZombie())
    {
      std::cout << "emcNoisyTowerFinder::FindHot: input file not found or is corrupted " << infilename.c_str() << std::endl;
      return;
    }
    h_hits_eta_phi_adc = (TH2F *) fin->Get(inHist.c_str());
    if (!h_hits_eta_phi_adc)
    {
      std::cout << "emcNoisyTowerFinder::FindHot: input hist not found " << inHist.c_str() << std::endl;
      delete fin;
      return;
    }
    h_hits_eta_phi_adc->SetDirectory(nullptr);  // Detach from the file to keep it in memory
    delete fin;
  }

  if (!h_hits_eta_phi_adc)
  {
    std::cout << "emcNoisyTowerFinder::FindHot: no valid histogram found" << std::endl;
    return;
  }
  TH2F *h_hits = h_hits_eta_phi_adc;

  TFile *fout = new TFile(outfilename.c_str(), "recreate");

  TH2F *h_hot = new TH2F("h_hot", "", Neta, 0, Neta, Nphi, 0, Nphi);
  TH2F *h_hitClean;  // = new TH2F("h_hitClean", "", Neta, 0, Neta, Nphi, 0, Nphi);
  TH2F *h_heatSigma = new TH2F("h_heatSigma", "", Neta, 0, Neta, Nphi, 0, Nphi);
  TH1F *h_perMedian = new TH1F("h_perMedian", "", 500, 0, 5);
  std::vector<TH1 *> h1_hits;
  h1_hits.resize(Neta);
  std::vector<TH1 *> h1_hits2;
  h1_hits2.resize(Neta);
  h_hits->Write();
  h_hitClean = (TH2F *) h_hits->Clone("h_hitClean");

  float max = h_hits->GetBinContent(h_hits->GetMaximumBin());
  float min = h_hits->GetMinimum();
  if (min == 0)
  {
    min = 1;
  }

  for (int ie = 0; ie < Neta; ie++)
  {
    h1_hits[ie] = new TH1F((std::string("h1_hits") + std::to_string(ie)).c_str(), "", 200, min, max);
    h1_hits2[ie] = new TH1F((std::string("h1_hits2_") + std::to_string(ie)).c_str(), "", 200, min, max);
    std::vector<float> vals;
    for (int iphi = 0; iphi < Nphi; iphi++)
    {
      float val = h_hits->GetBinContent(ie + 1, iphi + 1);
      h1_hits[ie]->Fill(val);
      vals.push_back(val);
    }

    float median = findMedian(vals);
    float mean = h1_hits[ie]->GetMean();
    float std = h1_hits[ie]->GetStdDev();
    for (int iphi = 0; iphi < Nphi; iphi++)
    {
      float val = h_hits->GetBinContent(ie + 1, iphi + 1);
      if (val / median > 5 || val / median < 0.1)
      {
        continue;
      }
      h1_hits2[ie]->Fill(val);
    }
    mean = h1_hits2[ie]->GetMean();
    std = h1_hits2[ie]->GetStdDev();
    for (int iphi = 0; iphi < Nphi; iphi++)
    {
      h_hot->SetBinContent(ie + 1, iphi + 1, 0);
      float val = h_hits->GetBinContent(ie + 1, iphi + 1);
      h_perMedian->Fill(val / median);
      float sigma = (val - mean) / std;
      h_heatSigma->SetBinContent(ie + 1, iphi + 1, sigma);
      if (sigma > sigma_bad_thresh)  // hot tower
      {
        h_hot->SetBinContent(ie + 1, iphi + 1, 2);
        h_hitClean->SetBinContent(ie + 1, iphi + 1, 0);
      }
      if ((sigma < -1 * sigma_bad_thresh || val < mean * percent_cold_thresh) && val != 0)  // cold tower
      {
        h_hot->SetBinContent(ie + 1, iphi + 1, 3);
        h_hitClean->SetBinContent(ie + 1, iphi + 1, 0);
      }
      if (val == 0)  // dead tower
      {
        h_hot->SetBinContent(ie + 1, iphi + 1, 1);
      }
    }
  }
  h_hitClean->Write();

  fout->Write();

  ////////////////////////////////////////////
  // make cdb tree
  size_t pos = outfilename.find_last_of('.');
  std::string f_cdbout_name = outfilename;
  f_cdbout_name.insert(pos, "cdb");

  CDBTTree *cdbttree_out = new CDBTTree(f_cdbout_name.c_str());
  std::string m_fieldname_out = "status";
  for (int i = 0; i < Neta; i++)
  {
    for (int j = 0; j < Nphi; j++)
    {
      unsigned int key = TowerInfoDefs::encode_emcal(i, j);
      if (Neta == 24)
      {
        key = TowerInfoDefs::encode_hcal(i, j);
      }
      int val = h_hot->GetBinContent(i + 1, j + 1);
      float sigma = h_heatSigma->GetBinContent(i + 1, j + 1);
      cdbttree_out->SetIntValue(key, m_fieldname_out, val);
      cdbttree_out->SetFloatValue(key, "CEMC_sigma", sigma);
    }
  }
  cdbttree_out->Commit();
  cdbttree_out->WriteCDBTTree();
  delete cdbttree_out;

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

float emcNoisyTowerFinder::findMedian(const std::vector<float> &arr)
{
  if (arr.empty())
  {
    return std::numeric_limits<float>::quiet_NaN();  // Return NaN if the array is empty
  }
  std::vector<float> sortedArr = arr;
  std::sort(sortedArr.begin(), sortedArr.end());

  size_t n = sortedArr.size();
  if (n % 2 == 0)
  {
    return (sortedArr[n / 2 - 1] + sortedArr[n / 2]) / 2.0;
  }
  else
  {
    return sortedArr[n / 2];
  }
}
