// $Id: $

/*!
 * \file QA_Draw_HCALOUT_TowerCluster.C
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
QA_Draw_Jet_TruthMatching(const char * jet =
    "h_QAG4SimJet_AntiKt_Truth_r07_AntiKt_Tower_r07",
    const char * qa_file_name_new = "data/G4sPHENIXCells_2000jets25GeV.root_qa.root",
    const char * qa_file_name_ref = "data/G4sPHENIXCells_250jets25GeV.root_qa.root")
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
  const double Nevent_new = 2000; // TODO: need to use normalization histos
  const double Nevent_ref = 250; // TODO: need to use normalization histos

  TCanvas *c1 = new TCanvas(
      TString("QA_Draw_Jet_TruthMatching_") + TString(jet),
      TString("QA_Draw_Jet_TruthMatching_") + TString(jet), 1800, 900);
  c1->Divide(3, 2);
  int idx = 1;
  TPad * p;

  p = (TPad *) c1->cd(idx++);
  c1->Update();
//  p->SetLogz();

    {

      TH2F * proj_new = (TH2F *) qa_file_new->GetObjectChecked(
          TString(jet) + "_Matching_dPhi", "TH2F");

      assert(proj_new);

      proj_new->Rebin2D(1, 5);

      TGraphErrors * ge = FitProfile(proj_new);

      proj_new->Draw("COLZ");
      ge->Draw("p");

    }
  TLine * l = new TLine(0, 0, 100, 00);
  l->Draw();

  p = (TPad *) c1->cd(idx++);
  c1->Update();
//  p->SetLogz();

    {

      TH2F * proj_new = (TH2F *) qa_file_new->GetObjectChecked(
          TString(jet) + "_Matching_dEta", "TH2F");

      assert(proj_new);

      proj_new->Rebin2D(1, 5);
      TGraphErrors * ge = FitProfile(proj_new);

      proj_new->Draw("COLZ");
      ge->Draw("p");
    }
  TLine * l = new TLine(0, 0, 100, 00);
  l->Draw();

  p = (TPad *) c1->cd(idx++);
  c1->Update();
//  p->SetLogz();

    {

      TH2F * proj_new = (TH2F *) qa_file_new->GetObjectChecked(
          TString(jet) + "_Matching_dE", "TH2F");

      assert(proj_new);

//    proj_new->Rebin2D(1,5);

      TGraphErrors * ge = FitProfile(proj_new);
      proj_new->Draw("COLZ");
      ge->Draw("p");
    }
  TLine * l = new TLine(0, 1, 100, 1);
  l->Draw();

  p = (TPad *) c1->cd(idx++);
  c1->Update();
//  p->SetLogz();

    {

      TH2F * proj_new = (TH2F *) qa_file_new->GetObjectChecked(
          TString(jet) + "_Matching_dEt", "TH2F");

      assert(proj_new);

//    proj_new->Rebin2D(1,5);
      TGraphErrors * ge = FitProfile(proj_new);
      proj_new->Draw("COLZ");
      ge->Draw("p");
    }
  TLine * l = new TLine(0, 1, 100, 1);
  l->Draw();

  p = (TPad *) c1->cd(idx++);
  c1->Update();
//  p->SetLogz();

    {

      TH2F * h2 = (TH2F *) qa_file_new->GetObjectChecked(
          TString(jet) + "_Matching_Count_Truth_Et", "TH2F");

      assert(h2);

      TH1 * h_norm = h2->ProjectionX(
          TString(jet) + "_Matching_Count_Truth_Et" + "_All", 1, 1);
      TH1 * h_pass = h2->ProjectionX(
          TString(jet) + "_Matching_Count_Truth_Et" + "_Matched", 2, 2);

      h_pass->Sumw2();

      h_pass->Divide(h_norm);

      h_pass->GetYaxis()->SetTitle("Reco efficiency");

      DrawReference(h_pass, NULL);
    }

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  //  p->SetLogz();

    {

      TH2F * h2 = (TH2F *) qa_file_new->GetObjectChecked(
          TString(jet) + "_Matching_Count_Reco_Et", "TH2F");

      assert(h2);

      TH1 * h_norm = h2->ProjectionX(
          TString(jet) + "_Matching_Count_Reco_Et" + "_All", 1, 1);
      TH1 * h_pass = h2->ProjectionX(
          TString(jet) + "_Matching_Count_Reco_Et" + "_Matched", 2, 2);

      h_pass->Sumw2();

      h_pass->Divide(h_norm);

      h_pass->GetYaxis()->SetTitle("Reconstruction Purity");

      DrawReference(h_pass, NULL);
    }

  SaveCanvas(c1, TString(qa_file_name_new) + TString(c1->GetName()), true);
}

