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

// ROOT
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <boost/format.hpp>

//________________________________
emcNoisyTowerFinder::emcNoisyTowerFinder(const std::string &name, const std::string &outputName, const std::string &cdbtreename_in, float adccut_sg_in, float adccut_k_in, float sigmas_lo_in, float sigmas_hi_in, float SG_f_in, float Kur_f_in, float region_f_in)
  : SubsysReco(name)
  , T(nullptr)
  , out(nullptr)
  , Outfile(outputName)
  , cdbtreename(cdbtreename_in)
  , adccut_sg(adccut_sg_in)  // The minimum ADC required for a tower with Saint Gobain fibers to register a hit
  , adccut_k(adccut_k_in)    // Minimum ADC for Kurary fibers to register a hit
  , sigmas_lo(sigmas_lo_in)  // # of standard deviations from the mode for a cold tower
  , sigmas_hi(sigmas_hi_in)  // # of standard deviations from the mode for a hot tower
  , SG_f(SG_f_in)            // Fiducial cut (artificial maximum) for Saint Gobain towers
  , Kur_f(Kur_f_in)          // Fiducial cut for Kurary
  , region_f(region_f_in)    // Fiducial cut for Sectors and IBs
{
  // initialize tree with SG vs Kurary fiber information
  fchannels = TFile::Open("channels.root", "read");

  channels = (TTree *) fchannels->Get("myTree");
  channels->SetBranchAddress("fiber_type", &fiber_type);

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
  TowerInfoContainer *emcTowerContainer;
  emcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_CEMC");
  if (!emcTowerContainer)
  {
    std::cout << PHWHERE << "emcNoisyTowerFinder::process_event Could not find node TOWERS_CEMC" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  // iterate through all towers, incrementing their Frequency arrays if they record a hit

  // bool goodevent = false;
  int tower_range = emcTowerContainer->size();
  for (int j = 0; j < tower_range; j++)
  {
    double energy = emcTowerContainer->get_tower_at_channel(j)->get_energy();
    channels->GetEntry(j);
    if (fabs(energy) > 0)
    {
      // keeping track of how many times a tower actually saw
      // some reasonable amount of energy and normalizing the response by
      // that as opposed to the raw number of events
      goodevents[j]++;
      goodeventsIB[j / nTowersIB]++;  // add these to the header file
      goodeventsSec[j / nTowersSec]++;
    }
    if (((fiber_type == 0) && (energy > adccut_sg)) ||
        ((fiber_type == 1) && (energy > adccut_k)))
    {
      towerF[j]++;
      // goodevent = true; //counter of events with nonzero EMCal energy
      sectorF[j / nTowersSec]++;
      ibF[j / nTowersIB]++;
    }

    towerE[j] += energy;
    sectorE[j / nSectors] += energy;
    ibE[j / nTowersIB] += energy;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
//__________________________
int emcNoisyTowerFinder::ResetEvent(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
//__________________________
void emcNoisyTowerFinder::FillHistograms(const int /*runnumber*/, const int /*segment*/)
{
  const int nBins = 101;

  out = new TFile(Outfile.c_str(), "RECREATE");

  // NOLINTNEXTLINE(bugprone-integer-division)
  Fspec = new TH2F("Fspec", "Fspec", nTowers + 1, -0.5, nTowers + 0.5, nBins, 0. - 1 / (2 * nBins), 1 + 1 / (2 * nBins));
  Fspec_SG = new TH2F("Fspec_SG", "Fspec_SG", nTowers + 1, -0.5, nTowers + 0.5, nBins, 0. - SG_f / (2 * nBins), SG_f + SG_f / (2 * nBins));
  Fspec_K = new TH2F("Fspec_K", "Fspec_K", nTowers + 1, -0.5, nTowers + 0.5, nBins, 0. - Kur_f / (2 * nBins), Kur_f + Kur_f / (2 * nBins));
  Fspec_sector = new TH2F("Fspec_sector", "Fspec_sector", nSectors + 1, -0.5, nSectors + 0.5, nBins, 0. - region_f / (2 * nBins), region_f + region_f / (2 * nBins));
  Fspec_IB = new TH2F("Fspec_IB", "Fspec_IB", nIB + 1, -0.5, nIB + 0.5, nBins, 0. - region_f / (2 * nBins), region_f + region_f / (2 * nBins));

  // NOLINTNEXTLINE(bugprone-integer-division)
  Fspeci = new TH2F("Fspeci", "Fspec_init", nTowers + 1, -0.5, nTowers + 0.5, nBins, 0. - 1 / (2 * nBins), 1);
  // NOLINTNEXTLINE(bugprone-integer-division)
  Fspeci_SG = new TH2F("Fspeci_SG", "Fspeci_SG", nTowers + 1, -0.5, nTowers + 0.5, nBins, 0. - 1 / (2 * nBins), 1);
  // NOLINTNEXTLINE(bugprone-integer-division)
  Fspeci_K = new TH2F("Fspeci_K", "Fspeci_K", nTowers + 1, -0.5, nTowers + 0.5, nBins, 0. - 1 / (2 * nBins), 1);
  // NOLINTNEXTLINE(bugprone-integer-division)
  Fspeci_sector = new TH2F("Fspeci_sector", "Fspeci_sector", nSectors + 1, -0.5, nSectors + 0.5, nBins, 0. - 1 / (2 * nBins), 1);
  // NOLINTNEXTLINE(bugprone-integer-division)
  Fspeci_IB = new TH2F("Fspeci_IB", "Fspeci_IB", nIB + 1, -0.5, nIB + 0.5, nBins, 0. - 1 / (2 * nBins), 1);

  // NOLINTNEXTLINE(bugprone-integer-division)
  Espec = new TH2F("Espec", "Espec", nTowers + 1, -0.5, nTowers - 0.5, nBins, 0. - 1 / (2 * nBins), 1);
  // NOLINTNEXTLINE(bugprone-integer-division)
  Espec_SG = new TH2F("Espec_SG", "Espec_SG", nTowers + 1, -0.5, nTowers - 0.5, nBins, 0. - 1 / (2 * nBins), 1);
  // NOLINTNEXTLINE(bugprone-integer-division)
  Espec_K = new TH2F("Espec_K", "Espec_K", nTowers + 1, -0.5, nTowers - 0.5, nBins, 0. - 1 / (2 * nBins), 1);
  // NOLINTNEXTLINE(bugprone-integer-division)
  Espec_sector = new TH2F("Espec_sector", "Espec_sector", nSectors, -0.5, nSectors - 0.5, nBins, 0. - 1 / (2 * nBins), 1);
  // NOLINTNEXTLINE(bugprone-integer-division)
  Espec_IB = new TH2F("Espec_IB", "Espec_IB", nTowers + 1, -0.5, nTowers - 0.5, nBins, 0. - 1 / (2 * nBins), 1);

  // hEventCounter = new TH1F("goodEventCounter","goodEventCounter",1,-0.5,0.5);
  // hEventCounter -> SetBinContent(1,goodevents);

  // fill histograms

  for (int i = 0; i < nTowers; i++)
  {
    Fspec->Fill(i, towerF[i] / goodevents[i]);
    Espec->Fill(i, towerE[i] / goodevents[i]);
    Fspeci->Fill(i, towerF[i] / goodevents[i]);

    channels->GetEntry(i);
    if (fiber_type == 0)
    {
      Fspec_SG->Fill(i, towerF[i] / goodevents[i]);
      Espec_SG->Fill(i, towerE[i] / goodevents[i]);
      Fspeci_SG->Fill(i, towerF[i] / goodevents[i]);
    }
    else
    {
      Fspec_K->Fill(i, towerF[i] / goodevents[i]);
      Espec_K->Fill(i, towerE[i] / goodevents[i]);
      Fspeci_K->Fill(i, towerF[i] / goodevents[i]);
    }
  }
  for (int j = 0; j < nIB; j++)
  {
    // int nTowerEvents = 1;
    // for(int evt = j*nTowersIB; evt < (j+1)*nTowersIB; j++)nTowerEvents+=goodevents[evt];
    //  Fspec_IB->Fill(j,1.0*ibF[j]/nTowersIB/(goodeventsIB[j]/nTowerEvents));
    //  Espec_IB->Fill(j,1.0*ibE[j]/nTowersIB/(goodeventsIB[j]/nTowerEvents));
    //  Fspeci_IB->Fill(j,1.0*ibF[j]/nTowersIB/(goodeventsIB[j]/nTowerEvents));
    // NOLINTNEXTLINE(bugprone-integer-division)
    Fspec_IB->Fill(j, 1.0 * ibF[j] / nTowersIB / (goodeventsIB[j] / nTowersIB));
    // NOLINTNEXTLINE(bugprone-integer-division)
    Espec_IB->Fill(j, 1.0 * ibE[j] / nTowersIB / (goodeventsIB[j] / nTowersIB));
    // NOLINTNEXTLINE(bugprone-integer-division)
    Fspeci_IB->Fill(j, 1.0 * ibF[j] / nTowersIB / (goodeventsIB[j] / nTowersIB));
  }
  for (int k = 0; k < nSectors; k++)
  {
    // int nTowerEvents = 1;
    // for(int evt = k*nTowersSec; evt < (k+1)*nTowersSec; k++)nTowerEvents+=goodevents[evt];
    //  Fspec_sector->Fill(k,1.0*sectorF[k]/nTowersSec/(goodeventsSec[k]/nTowerEvents));
    //  Espec_sector->Fill(k,1.0*sectorE[k]/nTowersSec/(goodeventsSec[k]/nTowerEvents));
    //  Fspeci_sector->Fill(k,1.0*sectorF[k]/nTowersSec/(goodeventsSec[k]/nTowerEvents));
    // NOLINTNEXTLINE(bugprone-integer-division)
    Fspec_sector->Fill(k, 1.0 * sectorF[k] / nTowersSec / (goodeventsSec[k] / nTowersSec));
    // NOLINTNEXTLINE(bugprone-integer-division)
    Espec_sector->Fill(k, 1.0 * sectorE[k] / nTowersSec / (goodeventsSec[k] / nTowersSec));
    // NOLINTNEXTLINE(bugprone-integer-division)
    Fspeci_sector->Fill(k, 1.0 * sectorF[k] / nTowersSec / (goodeventsSec[k] / nTowersSec));
  }

  // kill zeros:

  // // Fspec -> SetBinContent(1,0);
  // // Fspec_SG -> SetBinContent(1,0);
  // // Fspec_K -> SetBinContent(1,0);
  // // Fspec_IB -> SetBinContent(1,0);
  // // Fspec_sector -> SetBinContent(1,0);

  out->cd();

  // hEventCounter -> Write();
  // delete hEventCounter;

  Fspec->Write();
  delete Fspec;

  Fspec_K->Write();
  delete Fspec_K;

  Fspec_SG->Write();
  delete Fspec_SG;

  Fspec_sector->Write();
  delete Fspec_sector;

  Fspec_IB->Write();
  delete Fspec_IB;

  Espec->Write();
  delete Espec;

  Espec_SG->Write();
  delete Espec_SG;

  Espec_K->Write();
  delete Espec_K;

  Espec_sector->Write();
  delete Espec_sector;

  Espec_IB->Write();
  delete Espec_IB;

  Fspeci->Write();
  delete Fspeci;

  Fspeci_K->Write();
  delete Fspeci_K;

  Fspeci_SG->Write();
  delete Fspeci_SG;

  Fspeci_sector->Write();
  delete Fspeci_sector;

  Fspeci_IB->Write();
  delete Fspeci_IB;

  out->Close();
  delete out;
}
//__________________________
void emcNoisyTowerFinder::CalculateCutOffs(const int runnumber)
{
  TFile *histsIn = new TFile(boost::str(boost::format("output/%d/DST_CALOR-%08d_badTowerMapTree.root") % runnumber % runnumber).c_str());

  float cutoffFreq_SG;
  float cutoffFreq_K;
  float cutoffFreq;

  // float cutoffFreq_sector;
  float cutoffFreq_IB;

  float cutoffFreq_SG_lo;
  float cutoffFreq_K_lo;

  // float cutoffFreq_sector_lo;
  float cutoffFreq_IB_lo;

  // hEventCounter = (TH1F*)histsIn -> Get("goodEventCounter");
  // goodevents = hEventCounter -> GetBinContent(1);

  // calculate initial cut off values
  TH1F *dummyProj = nullptr;
  gStyle->SetOptStat(0);
  Fspec_SG = (TH2F *) histsIn->Get("Fspec_SG");
  dummyProj = (TH1F *) Fspec_SG->ProjectionY("dummy");
  dummyProj->SetBinContent(1, 0);
  cutoffFreq_SG = dummyProj->GetStdDev() * sigmas_hi + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());
  cutoffFreq_SG_lo = -1 * dummyProj->GetStdDev() * sigmas_lo + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());

  Fspec_K = (TH2F *) histsIn->Get("Fspec_K");
  dummyProj = (TH1F *) Fspec_K->ProjectionY("dummy");
  dummyProj->SetBinContent(1, 0);

  cutoffFreq_K = dummyProj->GetStdDev() * sigmas_hi + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());
  cutoffFreq_K_lo = -1 * Fspec_K->GetStdDev() * sigmas_lo + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());

  Fspec = (TH2F *) histsIn->Get("Fspec");
  dummyProj = (TH1F *) Fspec->ProjectionY("dummy");
  dummyProj->SetBinContent(1, 0);
  cutoffFreq = dummyProj->GetStdDev() * sigmas_hi + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());

  Fspec_IB = (TH2F *) histsIn->Get("Fspec_IB");
  dummyProj = (TH1F *) Fspec_IB->ProjectionY("dummy");
  dummyProj->SetBinContent(1, 0);

  cutoffFreq_IB = dummyProj->GetStdDev() * sigmas_hi + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());
  cutoffFreq_IB_lo = -1 * Fspec_IB->GetStdDev() * sigmas_lo + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());

  // Fspec_sector = (TH2F*)histsIn -> Get("Fspec_sector");
  // dummyProj = (TH1F*)Fspec_sector -> ProjectionY("dummy");
  // dummyProj -> SetBinContent(1,0);

  // cutoffFreq_sector = dummyProj->GetStdDev()*sigmas_hi + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());
  // cutoffFreq_sector_lo = -1*dummyProj->GetStdDev()*sigmas_lo + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());

  // get the rest of the histograms we nede
  Fspeci_K = (TH2F *) histsIn->Get("Fspeci_K");

  Fspeci_SG = (TH2F *) histsIn->Get("Fspeci_SG");
  Fspeci_SG->SetTitle(";Tower ID; Hit Frequency");

  Fspeci = (TH2F *) histsIn->Get("Fspeci");
  Fspeci_IB = (TH2F *) histsIn->Get("Fspeci_IB");
  Fspeci_IB->SetTitle(";IB ID; Hit Frequency");

  Fspeci_sector = (TH2F *) histsIn->Get("Fspeci_sector");

  // initial hot tower calculation
  for (int i = 0; i < nTowers; i++)
  {
    channels->GetEntry(i);

    TH1F *fSpecProj_SG = (TH1F *) Fspeci_SG->ProjectionY("dummySG", i + 1, i + 1);
    TH1F *fSpecProj_K = (TH1F *) Fspeci_K->ProjectionY("dummyK", i + 1, i + 1);
    TH1F *fSpecProj = (TH1F *) Fspeci->ProjectionY("dummySpec", i + 1, i + 1);

    if ((fiber_type == 0 && fSpecProj_SG->GetBinCenter(fSpecProj_SG->GetMaximumBin()) > cutoffFreq_SG) ||
        (fiber_type == 1 && fSpecProj_K->GetBinCenter(fSpecProj_K->GetMaximumBin()) > cutoffFreq_K))
    {
      hottowers[i]++;
    }
    else if (fSpecProj->GetBinCenter(fSpecProj->GetMaximumBin()) == 0)
    {
      deadtowers[i]++;
    }
    else if ((fiber_type == 0 && fSpecProj_SG->GetBinCenter(fSpecProj_SG->GetMaximumBin()) < cutoffFreq_SG_lo && fSpecProj_SG->GetBinCenter(fSpecProj_SG->GetMaximumBin()) > 0) ||
             (fiber_type == 1 && fSpecProj_K->GetBinCenter(fSpecProj_K->GetMaximumBin()) < cutoffFreq_K_lo && fSpecProj_K->GetBinCenter(fSpecProj_K->GetMaximumBin()) > 0))
    {
      coldtowers[i]++;
    }
  }
  // go through and look for hot/cold "regions", e.g. sectors and interface boards
  hot_regions = 0;
  cold_regions = 0;

  // Interface board QA, create projection

  // hotIB
  for (int j = 0; j < nIB; j++)
  {
    TH1F *fSpecIBProj = (TH1F *) Fspeci_IB->ProjectionY("dummyIB", j + 1, j + 1);

    if (fSpecIBProj->GetBinCenter(fSpecIBProj->GetMaximumBin()) > cutoffFreq_IB)
    {
      hot_regions = 1;
      std::cout << "IB " << j << " is hot with ADC rate" << fSpecIBProj->GetBinCenter(fSpecIBProj->GetMaximumBin()) << " > " << cutoffFreq_IB << std::endl;

      hotIB[j]++;
      for (int j1 = 0; j1 < nTowersIB; j1++)
      {
        hottowers[j * nTowersIB + j1]++;
        for (int biny = 0; biny < Fspeci->GetNbinsY(); biny++)
        {
          Fspeci->SetBinContent(j * nTowersIB + j1 + 1, biny + 1, -10);
        }
        // towerF[j*nTowersIB+j1]=0;//mask tower's contribution so it's no longer used to calculate cutoffs
      }
      for (int biny = 0; biny < Fspeci_IB->GetNbinsY(); biny++)
      {
        Fspeci_IB->SetBinContent(j + 1, biny + 1, -10);
      }
      // ibF[j]=-1;
    }
  }

  // cold IB
  for (int j = 0; j < nIB; j++)
  {
    TH1F *fSpecIBProj = (TH1F *) Fspeci_IB->ProjectionY("dummyIB", j + 1, j + 1);

    if (fSpecIBProj->GetBinCenter(fSpecIBProj->GetMaximumBin()) < cutoffFreq_IB_lo)
    {
      cold_regions = 1;
      std::cout << "IB " << j << " is cold with ADC rate" << fSpecIBProj->GetBinCenter(fSpecIBProj->GetMaximumBin()) << " < " << cutoffFreq_IB_lo << std::endl;
      hotIB[j]++;

      for (int j1 = 0; j1 < nTowersIB; j1++)
      {
        hottowers[j * nTowersIB + j1]++;
        for (int biny = 0; biny < Fspeci->GetNbinsY(); biny++)
        {
          Fspeci->SetBinContent(j * nTowersIB + j1 + 1, biny + 1, -10);
        }
        // towerF[j*nTowersIB+j1]=0;
      }
      // ibF[j] = -1;
      for (int biny = 0; biny < Fspeci_IB->GetNbinsY(); biny++)
      {
        Fspeci_IB->SetBinContent(j + 1, biny + 1, -10);
      }
    }
  }

  /*Removing sector level qa
  Essentially the issue is that removing IB's biases the determination of hot and cold sectors.
  Like removing one IB's worth of towers is just going to make a sector cold, or having 1 hot IB can make a sector hot.
  One could, in principle, re-scale the response by the fraction of missing IB's, but it isn't clear that that's necessary.
  Towers are necessarily effected by the behavior of their IB, but it isn't clear that a single sector can be driven
  bad by the behavior of a single IB. - AH*/

  // QA the sectors
  // hot sector
  //  for(int k = 0; k < nSectors; k++){
  //    TH1F *fSpecSecProj = (TH1F*)Fspeci_sector->ProjectionY("dummySec",k+1,k+1);

  //   if(fSpecSecProj->GetBinCenter(fSpecSecProj->GetMaximumBin()) > cutoffFreq_sector){
  //     hot_regions = 1;
  //     std::cout << "sector " << k << "is hot with ADC rate " << fSpecSecProj->GetBinCenter(fSpecSecProj->GetMaximumBin()) << " > " << cutoffFreq_sector  << std::endl;
  //     hotsectors[k]++;

  //     for(int k1 = 0; k1 < nTowersSec; k1++){
  // 	hottowers[k*nTowersSec+k1]++;
  // 	towerF[k*nTowersSec+k1] = 0;
  // 	//for(int biny = 0; biny < Fspeci->GetNbinsY(); biny++) Fspeci->SetBinContent(k*nTowersSec+k1, biny+1, 0);
  //     }
  //     sectorF[k] = -1;
  //     //fSpecSecProj->SetBinContent(k,0);
  //     //for(int biny = 0; biny < Fspeci_sector->GetNbinsY(); biny++) Fspeci_sector->SetBinContent(k+1, biny+1, 0);
  //   }
  // }

  // //cold sector
  // for(int k = 0; k < nSectors; k++){
  //   TH1F *fSpecSecProj = (TH1F*)Fspeci_sector->ProjectionY("dummySec",k+1,k+1);

  //   if(fSpecSecProj->GetBinCenter(fSpecSecProj->GetMaximumBin()) < cutoffFreq_sector_lo){
  //     cold_regions = 1;
  //     std::cout << "sector " << k << "is cold  with ADC rate" << fSpecSecProj->GetBinCenter(fSpecSecProj->GetMaximumBin()) << " < " << cutoffFreq_sector_lo << std::endl;
  //     coldsectors[k]++;

  //     for(int k1 = 0; k1 < nIB; k1++){
  // 	coldtowers[k*nIB+k1]++;
  // 	//for(int biny = 0; biny < Fspeci->GetNbinsY(); biny++) Fspeci->SetBinContent(k*nTowersSec+k1, biny+1, 0);
  // 	towerF[k*nTowersSec+k1] = 0;
  //     }
  //     //for(int biny = 0; biny < Fspeci_sector->GetNbinsY(); biny++) Fspeci_sector->SetBinContent(k+1, biny+1, 0);
  //     sectorF[k]=-1;
  //   }
  // }

  std::cout << "removing hot/cold regions" << std::endl;

  // go through again and calculate cut offs with bad regions masked until the calculation is made without hot regions
  while (hot_regions == 1 || cold_regions == 1)
  {
    std::cout << "hot/cold IB or sector detected. Running another pass for hot towers" << std::endl;
    // Fspec->Reset();
    // Fspec_SG->Reset();
    // Fspec_K->Reset();
    // Fspec_IB->Reset();
    // Fspec_sector->Reset();

    // not used for now
    //  Espec->Reset();
    //  Espec_SG->Reset();
    //  Espec_K->Reset();
    //  Espec_IB->Reset();
    //  Espec_sector->Reset();
    // for(int i = 0; i < nTowers; i++){

    // Fspec->Fill(i,towerF[i]/goodevents);
    //  Espec->Fill(i,towerE[i]/goodevents);
    // channels->GetEntry(i);
    // if(fiber_type == 0){
    // Fspec_SG->Fill(i,towerF[i]/goodevents);
    // Espec_SG->Fill(i,towerE[i]/goodevents);
    //}
    // else{
    // Fspec_K->Fill(i,towerF[i]/goodevents);
    // Espec_K->Fill(i,towerE[i]/goodevents);
    //}
    //}

    // for(int j = 0; j<nIB; j++){
    // Fspec_IB->Fill(j,1.0*ibF[j]/nTowersIB/goodevents);
    // Espec_IB->Fill(1,1.0*ibE[j]/64/goodevents);
    //}
    // for(int k = 0; k<nSectors; k++){
    // Fspec_sector->Fill(k,1.0*sectorF[k]/nTowersSec/goodevents);
    // Espec_sector->Fill(1,1.0*sectorE[k]/384/goodevents);
    //}

    // kill zeros

    // Fspec->SetBinContent(1,0);
    // Fspec_SG->SetBinContent(1,0);
    // Fspec_K->SetBinContent(1,0);
    // Fspec_IB->SetBinContent(1,0);
    // Fspec_sector->SetBinContent(1,0);

    dummyProj = (TH1F *) Fspec_SG->ProjectionY("dummy");
    dummyProj->SetBinContent(1, 0);

    cutoffFreq_SG = dummyProj->GetStdDev() * sigmas_hi + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());
    cutoffFreq_SG_lo = -1 * dummyProj->GetStdDev() * sigmas_lo + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());

    dummyProj = (TH1F *) Fspec_K->ProjectionY("dummy");
    dummyProj->SetBinContent(1, 0);

    cutoffFreq_K = dummyProj->GetStdDev() * sigmas_hi + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());
    cutoffFreq_K_lo = -1 * Fspec_K->GetStdDev() * sigmas_lo + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());

    dummyProj = (TH1F *) Fspec->ProjectionY("dummy");
    dummyProj->SetBinContent(1, 0);
    cutoffFreq = dummyProj->GetStdDev() * sigmas_hi + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());

    dummyProj = (TH1F *) Fspec_IB->ProjectionY("dummy");
    dummyProj->SetBinContent(1, 0);
    cutoffFreq_IB = dummyProj->GetStdDev() * sigmas_hi + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());
    cutoffFreq_IB_lo = -1 * Fspec_IB->GetStdDev() * sigmas_lo + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());

    // dummyProj = (TH1F*)Fspec_sector -> ProjectionY("dummy");
    // dummyProj -> SetBinContent(1,0);
    // cutoffFreq_sector = dummyProj->GetStdDev()*sigmas_hi + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());
    // cutoffFreq_sector_lo = -1*dummyProj->GetStdDev()*sigmas_lo + dummyProj->GetBinCenter(dummyProj->GetMaximumBin());

    for (int i = 0; i < nTowers; i++)
    {
      channels->GetEntry(i);
      TH1F *fSpecProj_SG = (TH1F *) Fspeci_SG->ProjectionY("dummySG", i + 1, i + 1);
      TH1F *fSpecProj_K = (TH1F *) Fspeci_K->ProjectionY("dummyK", i + 1, i + 1);
      TH1F *fSpecProj = (TH1F *) Fspeci->ProjectionY("dummySpec", i + 1, i + 1);

      if ((fiber_type == 0 && fSpecProj_SG->GetBinCenter(fSpecProj_SG->GetMaximumBin()) > cutoffFreq_SG) ||
          (fiber_type == 1 && fSpecProj_K->GetBinCenter(fSpecProj_K->GetMaximumBin()) > cutoffFreq_K))
      {
        hottowers[i]++;
      }
      else if (fSpecProj->GetBinCenter(fSpecProj->GetMaximumBin()) == 0)
      {
        deadtowers[i]++;
      }
      else if ((fiber_type == 0 && fSpecProj_SG->GetBinCenter(fSpecProj_SG->GetMaximumBin()) < cutoffFreq_SG_lo && fSpecProj_SG->GetBinCenter(fSpecProj_SG->GetMaximumBin()) > 0) ||
               (fiber_type == 1 && fSpecProj_K->GetBinCenter(fSpecProj_K->GetMaximumBin()) < cutoffFreq_K_lo && fSpecProj_K->GetBinCenter(fSpecProj_K->GetMaximumBin()) > 0))
      {
        coldtowers[i]++;
      }
    }

    hot_regions = 0;
    cold_regions = 0;

    // Re-QA IB's
    // hotIB
    for (int j = 0; j < nIB; j++)
    {
      if (ibF[j] == -1)
      {
        continue;  // IB has been marked hot already, don't bother with it
      }
      TH1F *fSpecIBProj = (TH1F *) Fspeci_IB->ProjectionY("dummyIB", j + 1, j + 1);

      if (fSpecIBProj->GetBinCenter(fSpecIBProj->GetMaximumBin()) > cutoffFreq_IB)
      {
        hot_regions = 1;
        std::cout << "IB " << j << "is hot with ADC rate " << fSpecIBProj->GetBinCenter(fSpecIBProj->GetMaximumBin()) << " > " << cutoffFreq_IB << std::endl;
        hotIB[j]++;

        for (int j1 = 0; j1 < nTowersIB; j1++)
        {
          hottowers[j * nTowersIB + j1]++;
          // for(int biny = 0; biny < Fspeci->GetNbinsY(); biny++) Fspeci -> SetBinContent(j*nTowersIB+j1+1,biny+1, 0);
          towerF[j * nTowersIB + j1] = 0;  // mask tower's contribution so it's no longer used to calculate cutoffs
        }
        // for(int biny = 0; biny < Fspeci_IB->GetNbinsY(); biny++)Fspeci_IB->SetBinContent(j+1,biny+1,0);
        ibF[j] = -1;
      }
    }

    // cold IB
    for (int j = 0; j < nIB; j++)
    {
      if (ibF[j] == -1)
      {
        continue;  // IB has been marked hot already, don't bother with it
      }

      TH1F *fSpecIBProj = (TH1F *) Fspeci_IB->ProjectionY("dummyIB", j + 1, j + 1);

      if (fSpecIBProj->GetBinCenter(fSpecIBProj->GetMaximumBin()) < cutoffFreq_IB_lo)
      {
        hot_regions = 1;
        std::cout << "IB " << j << " is cold with ADC rate " << fSpecIBProj->GetBinCenter(fSpecIBProj->GetMaximumBin()) << " < " << cutoffFreq_IB_lo << std::endl;
        hotIB[j]++;

        for (int j1 = 0; j1 < nTowersIB; j1++)
        {
          hottowers[j * nTowersIB + j1]++;
          // for(int biny = 0; biny < Fspeci->GetNbinsY(); biny++) Fspeci -> SetBinContent(j*nTowersIB+j1+1,biny+1,0);
          towerF[j * nTowersIB + j1] = 0;
        }
        ibF[j] = -1;
        // for(int biny = 0; biny < Fspeci_IB->GetNbinsY(); biny++)Fspeci_IB->SetBinContent(j+1, biny+1,0);
      }
    }

    // Re-QA sectors
    // hot sectors
    //  for(int k = 0; k<nSectors; k++){
    //    if(sectorF[k]==-1) continue;//sector has been marked hot already
    //    TH1F *fSpecSecProj = (TH1F*)Fspeci_sector->ProjectionY("dummySec",k+1,k+1);

    //   if(fSpecSecProj->GetBinCenter(fSpecSecProj->GetMaximumBin()) > cutoffFreq_sector){
    // 	hot_regions = 1;
    // 	std::cout << "sector " << k << " is hot with ADC rate " << fSpecSecProj->GetBinCenter(fSpecSecProj->GetMaximumBin()) << " > " <<cutoffFreq_sector  << std::endl;
    // 	hotsectors[k]++;

    // 	for(int k1 = 0; k1 < nTowersSec; k1++){
    // 	  hottowers[k*nTowersSec+k1]++;

    // 	  //for(int biny = 0; biny < Fspeci->GetNbinsY(); biny++) Fspeci->SetBinContent(k*nTowersSec+k1, biny+1, 0);
    // 	  towerF[k*nTowersSec+k1] = 0;
    // 	}
    // 	sectorF[k] = -1;
    // 	//fSpecSecProj->SetBinContent(k,0);
    // 	//for(int biny = 0; biny < Fspeci_sector->GetNbinsY(); biny++) Fspeci_sector->SetBinContent(k+1, biny+1, 0);
    //   }
    // }

    // cold sector
    //  for(int k = 0; k<nSectors; k++){
    //    if(sectorF[k]==-1) continue;
    //    TH1F *fSpecSecProj = (TH1F*)Fspeci_sector->ProjectionY("dummySec",k+1,k+1);

    //   if(fSpecSecProj->GetBinCenter(fSpecSecProj->GetMaximumBin()) <cutoffFreq_sector_lo){
    // 	cold_regions = 1;
    // 	std::cout << "sector " << k << " is cold  with ADC rate " << fSpecSecProj->GetBinCenter(fSpecSecProj->GetMaximumBin())  << " < " <<cutoffFreq_sector_lo << std::endl;
    // 	coldsectors[k]++;

    // 	for(int k1 = 0; k1 < nIB; k1++){
    // 	  coldtowers[k*nIB+k1]++;
    // 	  //for(int biny = 0; biny < Fspeci->GetNbinsY(); biny++) Fspeci->SetBinContent(k*nTowersSec+k1, biny+1, 0);
    // 	   towerF[k*nTowersSec+k1] = 0;
    // 	}
    // 	//for(int biny = 0; biny < Fspeci_sector->GetNbinsY(); biny++) Fspeci_sector->SetBinContent(k+1, biny+1, 0);
    // 	sectorF[k] = -1;
    //   }
    // }
  }

  std::cout << "emcNoisyTowerFinder::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  std::cout << "Saint Gobain Cutoff Frequency: " << cutoffFreq_SG << std::endl;
  std::cout << "Kurary Cutoff Frequency: " << cutoffFreq_K << std::endl;
  std::cout << "Overall Tower Cutoff Frequency: " << cutoffFreq << std::endl;
  std::cout << "IB Cutoff hit Frequency: " << cutoffFreq_IB << std::endl;
  // std::cout <<  "Sector Cutoff hit Frequency: " << cutoffFreq_sector << std::endl;
}

//____________________________________________________________________________..
void emcNoisyTowerFinder::WriteCDBTree(const int runnumber)
{
  // TFile *cdbOut = new TFile(Form("cdbMaps/%d/DST_CALOR-%08d_badTowerMapCDBTTree.root",runnumber,runnumber),"RECREATE");
  TFile *cdbOut = new TFile(boost::str(boost::format("cdbMaps/%d/CEMC_%08d_badTowerMapCDBTTree.root") % runnumber % runnumber).c_str(), "RECREATE");

  T = new TTree("T", "T");
  T->Branch("hot_channels", &m_hot_channels);

  CDBTTree *cdbttree = new CDBTTree(cdbtreename);
  std::string fieldname = "status";
  for (int i = 0; i < nTowers; i++)
  {
    unsigned int key = TowerInfoDefs::encode_emcal(i);

    if (hottowers[i] >= 0.5)
    {
      m_hot_channels = 2;
      T->Fill();
      cdbttree->SetIntValue(key, fieldname, 2);
    }
    else if (deadtowers[i] >= 0.5)
    {
      m_hot_channels = 1;
      T->Fill();
      cdbttree->SetIntValue(key, fieldname, 1);
    }
    else if (coldtowers[i] > 0.5)
    {
      m_hot_channels = 3;
      T->Fill();
      cdbttree->SetIntValue(key, fieldname, 3);
    }
    else
    {
      m_hot_channels = 0;
      T->Fill();
      cdbttree->SetIntValue(key, fieldname, 0);
    }
  }

  // fchannels->Close();
  // delete fchannels;
  // fchannels=NULL;

  cdbOut->cd();
  T->Write();
  delete T;
  cdbttree->Commit();
  // cdbttree->Print();
  cdbttree->WriteCDBTTree();
  delete cdbttree;
  cdbOut->Close();
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
