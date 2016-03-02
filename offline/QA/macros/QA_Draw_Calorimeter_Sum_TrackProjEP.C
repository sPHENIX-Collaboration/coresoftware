// $Id: $

/*!
 * \file QA_Draw_Calorimeter_Sum_TrackProjEP.C
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
QA_Draw_Calorimeter_Sum_TrackProjEP(
    const char * qa_file_name_new =
        "/phenix/u/jinhuang/links/ePHENIX_work/sPHENIX_work/production_analysis_updates/spacal1d/fieldmap/G4Hits_sPHENIX_pi-_eta0.30_32GeV-0000.root_qa.root",
    const char * qa_file_name_ref =
        "/phenix/u/jinhuang/links/ePHENIX_work/sPHENIX_work/production_analysis_updates/spacal1d/fieldmap/G4Hits_sPHENIX_pi+_eta0.30_32GeV-0000.root_qa.root")
//QA_Draw_Calorimeter_Sum_TrackProjEP(const char * qa_file_name_new =
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

  TCanvas *c1 = new TCanvas("QA_Draw_Calorimeter_Sum_TrackProjEP",
      "QA_Draw_Calorimeter_Sum_TrackProjEP", 1800, 600);
  c1->Divide(3, 1);
  int idx = 1;
  TPad * p;

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  p->SetLogy();

  if (Ntrack_new>0)
    {
      TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
          "h_QAG4Sim_CalorimeterSum_TrackProj_3x3Tower_EP", "TH1F");
      assert(h_new);

      h_new->Sumw2();
      h_new->Scale(1. / Ntrack_new);

      TH1F * h_ref = NULL;
      if (qa_file_ref)
        {
          TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
              "h_QAG4Sim_CalorimeterSum_TrackProj_3x3Tower_EP", "TH1F");
          assert(h_ref);

          h_ref->Scale(1. / Ntrack_ref);
        }

      h_new->GetYaxis()->SetTitleOffset(1.5);
      h_new->GetYaxis()->SetTitle("Count / track / bin");
      //      h_new->GetXaxis()->SetRangeUser(-0, .1);

      DrawReference(h_new, h_ref);
    }


  p = (TPad *) c1->cd(idx++);
  c1->Update();
  p->SetLogy();

  if (Ntrack_new>0)
    {
      TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
          "h_QAG4Sim_CalorimeterSum_TrackProj_5x5Tower_EP", "TH1F");
      assert(h_new);

      h_new->Sumw2();
      h_new->Scale(1. / Ntrack_new);

      TH1F * h_ref = NULL;
      if (qa_file_ref)
        {
          TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
              "h_QAG4Sim_CalorimeterSum_TrackProj_5x5Tower_EP", "TH1F");
          assert(h_ref);

          h_ref->Scale(1. / Ntrack_ref);
        }

      h_new->GetYaxis()->SetTitleOffset(1.5);
      h_new->GetYaxis()->SetTitle("Count / track / bin");
      //      h_new->GetXaxis()->SetRangeUser(-0, .1);

      DrawReference(h_new, h_ref);
    }

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  p->SetLogy();

  if (Nevent_new>0)
    {
      TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
          "h_QAG4Sim_CalorimeterSum_Cluster_EP", "TH1F");
      if (h_new)
        {

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
    }


  PutInputFileName(c1,0.07, qa_file_name_new, qa_file_name_ref);
  SaveCanvas(c1, TString(qa_file_name_new) + TString(c1->GetName()), true);
}

