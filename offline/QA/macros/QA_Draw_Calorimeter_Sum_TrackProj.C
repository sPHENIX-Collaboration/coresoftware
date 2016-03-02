// $Id: $

/*!
 * \file QA_Draw_Calorimeter_Sum_TrackProj.C
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include <cmath>
#include <vector>
#include <string>
#include <TFile.h>
#include <TString.h>
#include <TLine.h>
#include <TTree.h>
#include <cassert>

//some common style files
#include "SaveCanvas.C"
#include "SetOKStyle.C"
#include "QA_Draw_Utility.C"
using namespace std;

void
QA_Draw_Calorimeter_Sum_TrackProj(
    const char * qa_file_name_new =
        "/phenix/u/jinhuang/links/ePHENIX_work/sPHENIX_work/production_analysis_updates/spacal1d/fieldmap/G4Hits_sPHENIX_pi-_eta0.30_32GeV-0000.root_qa.root",
    const char * qa_file_name_ref =
        "/phenix/u/jinhuang/links/ePHENIX_work/sPHENIX_work/production_analysis_updates/spacal1d/fieldmap/G4Hits_sPHENIX_pi+_eta0.30_32GeV-0000.root_qa.root")
//QA_Draw_Calorimeter_Sum_TrackProj(const char * qa_file_name_new =
//    "data/G4sPHENIXCells_100e24GeV.root_qa.root",
//    const char * qa_file_name_ref =
//        "data/G4Hits_sPHENIX_e-_eta0_24GeV-0000.root_qa.root")
{

  SetOKStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  TVirtualFitter::SetDefaultFitter("Minuit2");

  TFile * qa_file_new = new TFile(qa_file_name_new);
  assert(qa_file_new->IsOpen());

  TFile * qa_file_ref = NULL;
  if (qa_file_name_ref)
    {
      qa_file_ref = new TFile(qa_file_name_ref);
      assert(qa_file_ref->IsOpen());
    }

  // obtain normalization
  double Nevent_new = 1;
  double Nevent_ref = 1;

  // obtain normalization
  double Ntrack_new = 0;
  double Ntrack_ref = 0;

  if (qa_file_new)
    {
      TH1D * h_norm = (TH1D *) qa_file_new->GetObjectChecked(
          TString("h_QAG4Sim_CalorimeterSum_Normalization"), "TH1D");
      assert(h_norm);

      Nevent_new = h_norm->GetBinContent(h_norm->GetXaxis()->FindBin("Event"));
      Ntrack_new = h_norm->GetBinContent(h_norm->GetXaxis()->FindBin("Track"));
    }
  if (qa_file_ref)
    {
      TH1D * h_norm = (TH1D *) qa_file_ref->GetObjectChecked(
          TString("h_QAG4Sim_CalorimeterSum_Normalization"), "TH1D");
      assert(h_norm);

      Nevent_ref = h_norm->GetBinContent(h_norm->GetXaxis()->FindBin("Event"));
      Ntrack_ref = h_norm->GetBinContent(h_norm->GetXaxis()->FindBin("Track"));
    }

  vector<string> subsystems;
  subsystems.push_back("CEMC");
  subsystems.push_back("HCALIN");
  subsystems.push_back("HCALOUT");

  TCanvas *c1 = new TCanvas("QA_Draw_Calorimeter_Sum_TrackProj",
      "QA_Draw_Calorimeter_Sum_TrackProj", 1100, 1100);
  c1->Divide(3, 3);
  int idx = 1;
  TPad * p;

  for (int i = 0; i < subsystems.size(); ++i)
    {
      const TString subsystem(subsystems[i].c_str());

      p = (TPad *) c1->cd(idx++);
      c1->Update();
      p->SetLogz();

      TH2F * h_QAG4Sim_CalorimeterSum_TrackProj =
          (TH2F *) qa_file_new->GetObjectChecked(
              "h_QAG4Sim_CalorimeterSum_" + subsystem + "_TrackProj", "TH2F");
      assert(h_QAG4Sim_CalorimeterSum_TrackProj);
      h_QAG4Sim_CalorimeterSum_TrackProj->GetYaxis()->SetTitleOffset(1.5);
      h_QAG4Sim_CalorimeterSum_TrackProj->Draw("COLZ");

      TH2F * h_QAG4Sim_CalorimeterSum_TrackProj_Ref = NULL;
      if (qa_file_ref)
        {
          h_QAG4Sim_CalorimeterSum_TrackProj_Ref =
              (TH2F *) qa_file_ref->GetObjectChecked(
                  "h_QAG4Sim_CalorimeterSum_" + subsystem + "_TrackProj",
                  "TH2F");
          assert(h_QAG4Sim_CalorimeterSum_TrackProj);
        }

      p = (TPad *) c1->cd(idx++);
      c1->Update();
      p->SetLogz();

      // x projections
        {
          TH1 * h_new = h_QAG4Sim_CalorimeterSum_TrackProj->ProjectionX(h_QAG4Sim_CalorimeterSum_TrackProj->GetName()+TString("_px"));

          h_new->Scale(1. / Ntrack_new);

          TH1 * h_ref = NULL;
          if (h_QAG4Sim_CalorimeterSum_TrackProj_Ref)
            {
              TH1 * h_ref = h_QAG4Sim_CalorimeterSum_TrackProj_Ref->ProjectionX(h_QAG4Sim_CalorimeterSum_TrackProj->GetName()+TString("_px_ref"));

              h_ref->Scale(1. / Ntrack_ref);
            }

          h_new->GetXaxis()->SetTitleOffset(1.05);
          h_new->GetYaxis()->SetTitle("Energy / track / bin");

          DrawReference(h_new, h_ref);
        }

        p = (TPad *) c1->cd(idx++);
        c1->Update();
        p->SetLogz();

        // y projections
          {
            TH1 * h_new = h_QAG4Sim_CalorimeterSum_TrackProj->ProjectionY(h_QAG4Sim_CalorimeterSum_TrackProj->GetName()+TString("_py"));

            h_new->Scale(1. / Ntrack_new);

            TH1 * h_ref = NULL;
            if (h_QAG4Sim_CalorimeterSum_TrackProj_Ref)
              {
                TH1 * h_ref = h_QAG4Sim_CalorimeterSum_TrackProj_Ref->ProjectionY(h_QAG4Sim_CalorimeterSum_TrackProj->GetName()+TString("_py_ref"));

                h_ref->Scale(1. / Ntrack_ref);
              }

            h_new->GetXaxis()->SetTitleOffset(1.05);
            h_new->GetYaxis()->SetTitle("Energy / track / bin");

            DrawReference(h_new, h_ref);
          }

    }

  PutInputFileName(c1, .03, qa_file_name_new, qa_file_name_ref);

  SaveCanvas(c1, TString(qa_file_name_new) + TString(c1->GetName()), true);
}

