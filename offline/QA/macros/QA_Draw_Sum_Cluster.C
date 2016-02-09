// $Id: $

/*!
 * \file QA_Draw_Sum_Cluster.C
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include <cmath>
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
QA_Draw_Sum_Cluster(const char * qa_file_name_new =
    "data/G4sPHENIXCells_1000pi24GeV.root_qa.root",
    const char * qa_file_name_ref =
        "data/G4Hits_sPHENIX_pi-_eta0_24GeV-0000.root_qa.root")
//QA_Draw_Sum_Cluster(const char * qa_file_name_new =
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

  if (qa_file_new)
    {
      TH1D * h_norm = (TH1D *) qa_file_new->GetObjectChecked(
          TString("h_QAG4Sim_CalorimeterSum_Normalization"), "TH1D");
      assert(h_norm);

      Nevent_new = h_norm->GetBinContent(h_norm->GetXaxis()->FindBin("Event"));
    }
  if (qa_file_ref)
    {
      TH1D * h_norm = (TH1D *) qa_file_ref->GetObjectChecked(
          TString("h_QAG4Sim_CalorimeterSum_Normalization"), "TH1D");
      assert(h_norm);

      Nevent_ref = h_norm->GetBinContent(h_norm->GetXaxis()->FindBin("Event"));
    }

  TCanvas *c1 = new TCanvas("QA_Draw_Sum_Cluster", "QA_Draw_Sum_Cluster", 1800,
      900);
  c1->Divide(3, 2);
  int idx = 1;
  TPad * p;

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  p->SetLogz();

    {

      TH2F * h_QAG4Sim_CalorimeterSum_Cluster_CEMC_HCALIN =
          (TH2F *) qa_file_new->GetObjectChecked(
              "h_QAG4Sim_CalorimeterSum_Cluster_CEMC_HCALIN", "TH2F");
      assert(h_QAG4Sim_CalorimeterSum_Cluster_CEMC_HCALIN);
      h_QAG4Sim_CalorimeterSum_Cluster_CEMC_HCALIN->GetYaxis()->SetTitleOffset(
          1.5);
      h_QAG4Sim_CalorimeterSum_Cluster_CEMC_HCALIN->Draw("COLZ");

    }

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  p->SetLogz();

    {

      TH2F * h_QAG4Sim_CalorimeterSum_Cluster_CEMC_HCALIN_HCALOUT =
          (TH2F *) qa_file_new->GetObjectChecked(
              "h_QAG4Sim_CalorimeterSum_Cluster_CEMC_HCALIN_HCALOUT", "TH2F");
      assert(h_QAG4Sim_CalorimeterSum_Cluster_CEMC_HCALIN_HCALOUT);
      h_QAG4Sim_CalorimeterSum_Cluster_CEMC_HCALIN_HCALOUT->GetYaxis()->SetTitleOffset(
          1.5);
      h_QAG4Sim_CalorimeterSum_Cluster_CEMC_HCALIN_HCALOUT->Draw("COLZ");

    }

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  // empty pannel

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  p->SetLogy();

    {
      TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
          "h_QAG4Sim_CalorimeterSum_Cluster_Ratio_CEMC_HCALIN", "TH1F");
      assert(h_new);

      h_new->Sumw2();
      h_new->Scale(1. / Nevent_new);

      TH1F * h_ref = NULL;
      if (qa_file_ref)
        {
          TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
              "h_QAG4Sim_CalorimeterSum_Cluster_Ratio_CEMC_HCALIN", "TH1F");
          assert(h_ref);

          h_ref->Scale(1. / Nevent_ref);
        }

      h_new->GetYaxis()->SetTitleOffset(1.5);
      h_new->GetYaxis()->SetTitle("Count / event / bin");
      //      h_new->GetXaxis()->SetRangeUser(-0, .1);

      DrawReference(h_new, h_ref);
    }

    p = (TPad *) c1->cd(idx++);
    c1->Update();
    p->SetLogy();

      {
        TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
            "h_QAG4Sim_CalorimeterSum_Cluster_Ratio_CEMC_HCALIN_HCALOUT", "TH1F");
        assert(h_new);

        h_new->Sumw2();
        h_new->Scale(1. / Nevent_new);

        TH1F * h_ref = NULL;
        if (qa_file_ref)
          {
            TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
                "h_QAG4Sim_CalorimeterSum_Cluster_Ratio_CEMC_HCALIN_HCALOUT", "TH1F");
            assert(h_ref);

            h_ref->Scale(1. / Nevent_ref);
          }

        h_new->GetYaxis()->SetTitleOffset(1.5);
        h_new->GetYaxis()->SetTitle("Count / event / bin");
        //      h_new->GetXaxis()->SetRangeUser(-0, .1);

        DrawReference(h_new, h_ref);
      }

      p = (TPad *) c1->cd(idx++);
      c1->Update();
      p->SetLogy();

        {
          TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
              "h_QAG4Sim_CalorimeterSum_Cluster_EP", "TH1F");
          assert(h_new);

          h_new->Sumw2();
          h_new->Scale(1. / Nevent_new);

          TH1F * h_ref = NULL;
          if (qa_file_ref)
            {
              TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
                  "h_QAG4Sim_CalorimeterSum_Cluster_EP", "TH1F");
              assert(h_ref);

              h_ref->Scale(1. / Nevent_ref);
            }

          h_new->GetYaxis()->SetTitleOffset(1.5);
          h_new->GetYaxis()->SetTitle("Count / event / bin");
          //      h_new->GetXaxis()->SetRangeUser(-0, .1);

          DrawReference(h_new, h_ref);
        }
  SaveCanvas(c1, TString(qa_file_name_new) + TString(c1->GetName()), true);
}

