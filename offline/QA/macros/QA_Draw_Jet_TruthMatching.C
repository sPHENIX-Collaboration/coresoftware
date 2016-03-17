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
#include <TGraphErrors.h>
#include <TLine.h>
#include <TTree.h>
#include <cassert>
#include <vector>

//some common style files
#include "SaveCanvas.C"
#include "SetOKStyle.C"
#include "QA_Draw_Utility.C"
using namespace std;

vector<TGraphErrors *>
QA_Draw_Jet_TruthMatching(const char * jet =
    "h_QAG4SimJet_AntiKt_Truth_r07_AntiKt_Tower_r07",
    const char * qa_file_name_new =
        "data/G4sPHENIXCells_250jets25GeV.root_qa.root",
    const char * qa_file_name_ref =
        "data/G4sPHENIXCells_2000jets25GeV.root_qa.root")
{
  //! drawing energy range
  const double min_Et = 10;
  const double max_Et = 80;

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

  vector<TGraphErrors *> resolution_collections;

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

      proj_new->GetXaxis()->SetRangeUser(min_Et, max_Et);
      proj_new->GetYaxis()->SetTitleOffset(1.5);
      proj_new->Draw("COLZ");

      TGraphErrors * ge_ref = NULL;
      if (qa_file_ref)
        {
          TH2F * proj_ref = (TH2F *) qa_file_ref->GetObjectChecked(
              TString(jet) + "_Matching_dPhi", "TH2F");
          assert(proj_ref);
          proj_ref->Rebin2D(1, 5);
          TGraphErrors * ge_ref = FitProfile(proj_ref);
        }
      DrawReference(ge, ge_ref);

      resolution_collections.push_back(ge);

    }
  TLine * l = new TLine(min_Et, 0, max_Et, 00);
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

      proj_new->GetXaxis()->SetRangeUser(min_Et, max_Et);
      proj_new->GetYaxis()->SetTitleOffset(1.5);
      proj_new->Draw("COLZ");

      TGraphErrors * ge_ref = NULL;
      if (qa_file_ref)
        {
          TH2F * proj_ref = (TH2F *) qa_file_ref->GetObjectChecked(
              TString(jet) + "_Matching_dEta", "TH2F");
          assert(proj_ref);
          proj_ref->Rebin2D(1, 5);
          TGraphErrors * ge_ref = FitProfile(proj_ref);
        }
      DrawReference(ge, ge_ref);

      resolution_collections.push_back(ge);
    }
  TLine * l = new TLine(min_Et, 0, max_Et, 00);
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

      proj_new->GetXaxis()->SetRangeUser(min_Et, max_Et);
      proj_new->GetYaxis()->SetTitleOffset(1.5);
      proj_new->Draw("COLZ");

      TGraphErrors * ge_ref = NULL;
      if (qa_file_ref)
        {
          TH2F * proj_ref = (TH2F *) qa_file_ref->GetObjectChecked(
              TString(jet) + "_Matching_dE", "TH2F");
          assert(proj_ref);
          proj_ref->Rebin2D(1, 5);
          TGraphErrors * ge_ref = FitProfile(proj_ref);
        }
      DrawReference(ge, ge_ref);

      resolution_collections.push_back(ge);
    }
  TLine * l = new TLine(min_Et, 1, max_Et, 1);
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

      proj_new->GetXaxis()->SetRangeUser(min_Et, max_Et);
      proj_new->GetYaxis()->SetTitleOffset(1.5);
      proj_new->Draw("COLZ");

      TGraphErrors * ge_ref = NULL;
      if (qa_file_ref)
        {
          TH2F * proj_ref = (TH2F *) qa_file_ref->GetObjectChecked(
              TString(jet) + "_Matching_dEt", "TH2F");
          assert(proj_ref);
          proj_ref->Rebin2D(1, 5);
          TGraphErrors * ge_ref = FitProfile(proj_ref);
        }
      DrawReference(ge, ge_ref);

      resolution_collections.push_back(ge);
    }
  TLine * l = new TLine(min_Et, 1, max_Et, 1);
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
      //          TH1 * h_pass = h2->ProjectionX(
      //              TString(jet) + "_Matching_Count_Truth_Et" + "_Matched", 2, 2);// inclusive match
      TH1 * h_pass = h2->ProjectionX(
          TString(jet) + "_Matching_Count_Truth_Et" + "_Matched", 3, 3); // unique match
      assert(h_norm);
      assert(h_pass);
      TH1 * h_ratio = GetBinominalRatio(h_pass, h_norm);

      h_ratio->GetXaxis()->SetRangeUser(min_Et, max_Et);
      h_ratio->GetYaxis()->SetTitle("Reco efficiency");
      h_ratio->GetYaxis()->SetRangeUser(-0, 1.2);

      TH1 * h_ratio_ref = NULL;
      if (qa_file_ref)
        {
          TH2F * h2 = (TH2F *) qa_file_ref->GetObjectChecked(
              TString(jet) + "_Matching_Count_Truth_Et", "TH2F");
          assert(h2);
          TH1 * h_norm = h2->ProjectionX(
              TString(jet) + "_Matching_Count_Truth_Et" + "_All", 1, 1);
//          TH1 * h_pass = h2->ProjectionX(
//              TString(jet) + "_Matching_Count_Truth_Et" + "_Matched", 2, 2);// inclusive match
          TH1 * h_pass = h2->ProjectionX(
              TString(jet) + "_Matching_Count_Truth_Et" + "_Matched", 3, 3); // unique match
          assert(h_norm);
          assert(h_pass);
          h_ratio_ref = GetBinominalRatio(h_pass, h_norm, true);
        }

      DrawReference(h_ratio, h_ratio_ref, true);

      resolution_collections.push_back(new TGraphErrors(h_ratio));
    }

  TLine * l = new TLine(min_Et, 1, max_Et, 1);
  l->Draw();

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  //  p->SetLogz();

    {

      TH2F * h2 = (TH2F *) qa_file_new->GetObjectChecked(
          TString(jet) + "_Matching_Count_Reco_Et", "TH2F");
      assert(h2);

      TH1 * h_norm = h2->ProjectionX(
          TString(jet) + "_Matching_Count_Reco_Et" + "_All", 1, 1);
//      TH1 * h_pass = h2->ProjectionX(
//          TString(jet) + "_Matching_Count_Reco_Et" + "_Matched", 2, 2); // inclusive match
      TH1 * h_pass = h2->ProjectionX(
          TString(jet) + "_Matching_Count_Reco_Et" + "_Matched", 3, 3); // unique match
      assert(h_norm);
      assert(h_pass);
      TH1 * h_ratio = GetBinominalRatio(h_pass, h_norm);

      h_ratio->GetXaxis()->SetRangeUser(min_Et, max_Et);
      h_ratio->GetYaxis()->SetTitle("Reconstruction Purity");
      h_ratio->GetYaxis()->SetRangeUser(-0, 1.2);

      TH1 * h_ratio_ref = NULL;
      if (qa_file_ref)
        {
          TH2F * h2 = (TH2F *) qa_file_ref->GetObjectChecked(
              TString(jet) + "_Matching_Count_Reco_Et", "TH2F");
          assert(h2);

          TH1 * h_norm = h2->ProjectionX(
              TString(jet) + "_Matching_Count_Reco_Et" + "_All", 1, 1);
          //      TH1 * h_pass = h2->ProjectionX(
          //          TString(jet) + "_Matching_Count_Reco_Et" + "_Matched", 2, 2); // inclusive match
          TH1 * h_pass = h2->ProjectionX(
              TString(jet) + "_Matching_Count_Reco_Et" + "_Matched", 3, 3); // unique match
          assert(h_norm);
          assert(h_pass);
          h_ratio_ref = GetBinominalRatio(h_pass, h_norm, true);
        }

      DrawReference(h_ratio, h_ratio_ref, true);

      resolution_collections.push_back(new TGraphErrors(h_ratio));
    }

  TLine * l = new TLine(min_Et, 1, max_Et, 1);
  l->Draw();

  PutInputFileName(c1, .04, qa_file_name_new, qa_file_name_ref);
  SaveCanvas(c1, TString(qa_file_name_new) + TString(c1->GetName()), true);

  return resolution_collections;
}

