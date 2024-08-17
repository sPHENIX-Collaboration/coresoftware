#include "TFile.h"
#include "TTree.h"

float PiRange(float deltaPhi)
{
  if(deltaPhi > M_PI) deltaPhi -= 2*M_PI;
  if(deltaPhi < -M_PI) deltaPhi += 2*M_PI;

  return deltaPhi;
}

void saveroot(int runnumber)
{
  gStyle->SetOptStat(0);

  TString inputfile = Form("./Reconstructed/%d/final_%d_ana.root",runnumber,runnumber);

  TFile *file = new TFile(inputfile, "READ");
  TTree *tree = (TTree*)file->Get("tree");

  std::vector<float> *_track_phi_emc = 0;
  std::vector<float> *_track_z_emc = 0;
  std::vector<float> *_emcal_phi = 0;
  std::vector<float> *_emcal_z = 0;

  tree->SetBranchAddress("_track_phi_emc", &_track_phi_emc);
  tree->SetBranchAddress("_track_z_emc", &_track_z_emc);
  tree->SetBranchAddress("_emcal_phi", &_emcal_phi);
  tree->SetBranchAddress("_emcal_z", &_emcal_z);

  TFile* outputfile = new TFile(Form("root/%d.root",runnumber),"recreate");
  TTree* outtree = new TTree("tree","");

  float m_dphi, m_dz, m_calo_z, m_track_z;

  outtree->Branch("dphi",&m_dphi,"dphi/F");
  outtree->Branch("dz",&m_dz,"dz/F");
  outtree->Branch("calo_z",&m_calo_z,"calo_z/F");
  outtree->Branch("track_z",&m_track_z,"track_z/F");

  for(int i = 0; i < tree->GetEntries(); i++)
  {
    tree->GetEntry(i);

    for(unsigned int itrack = 0; itrack < _track_phi_emc->size(); itrack++)
    {
      if(isnan(_track_phi_emc->at(itrack)))
      {
        continue;
      }

      std::pair<float, float> TrackProjsEMCal;
      TrackProjsEMCal = std::make_pair(_track_phi_emc->at(itrack), _track_z_emc->at(itrack));

      for(unsigned int iem = 0; iem < _emcal_z->size(); iem++)
      {
        std::pair<float, float> EMCalPos;
        EMCalPos = std::make_pair(_emcal_phi->at(iem), _emcal_z->at(iem));

        m_dphi = PiRange(TrackProjsEMCal.first - EMCalPos.first);
        m_dz = TrackProjsEMCal.second - EMCalPos.second;
        m_calo_z = EMCalPos.second;
        m_track_z = TrackProjsEMCal.second;
        outtree->Fill();
      }
    }
  }
  outputfile->Write();
  outputfile->Close();
}
