// $Id: $

/*!
 * \file QA_Draw_CEMC_G4Hit.C
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
QA_Draw_CEMC_G4Hit(
    const char * qa_file_name_new =
        "/phenix/u/jinhuang/links/ePHENIX_work/sPHENIX_work/production_analysis_updates/spacal1d/fieldmap/G4Hits_sPHENIX_pi-_eta0.30_32GeV-0000.root_qa.root",
    const char * qa_file_name_ref =
        "/phenix/u/jinhuang/links/ePHENIX_work/sPHENIX_work/production_analysis_updates/spacal1d/fieldmap/G4Hits_sPHENIX_pi+_eta0.30_32GeV-0000.root_qa.root")

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

  TCanvas *c1 = new TCanvas("QA_Draw_CEMC_G4Hit", "QA_Draw_CEMC_G4Hit", 1800,
      900);
  c1->Divide(4, 2);
  int idx = 1;
  TPad * p;

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  p->SetLogz();

  TH2F * h_QAG4Sim_CEMC_G4Hit_XY = (TH2F *) qa_file_new->GetObjectChecked(
      "h_QAG4Sim_CEMC_G4Hit_XY", "TH2F");
  assert(h_QAG4Sim_CEMC_G4Hit_XY);
  h_QAG4Sim_CEMC_G4Hit_XY->GetYaxis()->SetTitleOffset(1.5);
  h_QAG4Sim_CEMC_G4Hit_XY->Draw("COLZ");

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  p->SetLogz();

  TH2F * h_QAG4Sim_CEMC_G4Hit_RZ = (TH2F *) qa_file_new->GetObjectChecked(
      "h_QAG4Sim_CEMC_G4Hit_RZ", "TH2F");
  assert(h_QAG4Sim_CEMC_G4Hit_RZ);
  h_QAG4Sim_CEMC_G4Hit_RZ->GetYaxis()->SetTitleOffset(1.5);
  h_QAG4Sim_CEMC_G4Hit_RZ->Draw("COLZ");

  p = (TPad *) c1->cd(idx++);
  c1->Update();
//  p->SetLogz();

    {

      TH2F * h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection =
          (TH2F *) qa_file_new->GetObjectChecked(
              "h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection", "TH2F");
      assert(h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection);

      TH1D * proj_new =
          h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection->ProjectionX(
              "qa_file_new_h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection_px");

      proj_new->Scale(1. / proj_new->GetSum());

      TH1D * proj_ref = NULL;
      if (qa_file_ref)
        {
          TH2F * h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection =
              (TH2F *) qa_file_ref->GetObjectChecked(
                  "h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection", "TH2F");
          assert(h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection);

          proj_ref = h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection->ProjectionX(
              "qa_file_ref_h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection_px");
          proj_ref->Scale(1. / proj_ref->GetSum());

        }

      proj_new->GetYaxis()->SetTitle("Normalized energy distribution");
      proj_new->GetXaxis()->SetRangeUser(-10, 10);

      DrawReference(proj_new, proj_ref);
    }

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  //  p->SetLogz();

    {

      TH2F * h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection =
          (TH2F *) qa_file_new->GetObjectChecked(
              "h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection", "TH2F");
      assert(h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection);

      TH1D * proj_new =
          h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection->ProjectionY(
              "qa_file_new_h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection_py");

      proj_new->Scale(1. / proj_new->GetSum());

      TH1D * proj_ref = NULL;
      if (qa_file_ref)
        {
          TH2F * h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection =
              (TH2F *) qa_file_ref->GetObjectChecked(
                  "h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection", "TH2F");
          assert(h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection);

          proj_ref = h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection->ProjectionY(
              "qa_file_ref_h_QAG4Sim_CEMC_G4Hit_LateralTruthProjection_py");
          proj_ref->Scale(1. / proj_ref->GetSum());

        }

      proj_new->GetYaxis()->SetTitle("Normalized energy distribution");
      proj_new->GetXaxis()->SetRangeUser(-10, 10);

      DrawReference(proj_new, proj_ref);
    }

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  p->SetLogx();
  p->SetLogy();

    {

      TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
          "h_QAG4Sim_CEMC_G4Hit_HitTime", "TH1F");
      assert(h_new);

      h_new->Scale(1. / h_new->GetSum());

      TH1F * h_ref = NULL;
      if (qa_file_ref)
        {
          TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
              "h_QAG4Sim_CEMC_G4Hit_HitTime", "TH1F");
          assert(h_ref);

          h_ref->Scale(1. / h_ref->GetSum());
        }

      h_new->GetYaxis()->SetTitleOffset(1.5);
      h_new->GetYaxis()->SetTitle("Normalized energy per bin");
//      h_new->GetXaxis()->SetRangeUser(-0, .1);

      DrawReference(h_new, h_ref);
    }

  p = (TPad *) c1->cd(idx++);
  c1->Update();
//    p->SetLogx();
  p->SetLogy();

    {

      TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
          "h_QAG4Sim_CEMC_G4Hit_FractionTruthEnergy", "TH1F");
      assert(h_new);

      h_new->Rebin(20);
      h_new->Sumw2();
      h_new->Scale(1. / h_new->GetSum());

      TH1F * h_ref = NULL;
      if (qa_file_ref)
        {
          TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
              "h_QAG4Sim_CEMC_G4Hit_FractionTruthEnergy", "TH1F");
          assert(h_ref);

          h_ref->Rebin(20);
          h_ref->Scale(1. / h_ref->GetSum());
        }

      h_new->GetYaxis()->SetTitleOffset(1.5);
      h_new->GetYaxis()->SetTitle("Probability per bin");
      //      h_new->GetXaxis()->SetRangeUser(-0, .1);

      DrawReference(h_new, h_ref);
    }

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  //  p->SetLogz();

    {

      TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
          "h_QAG4Sim_CEMC_G4Hit_VSF", "TH1F");
      assert(h_new);

      h_new->Rebin(2);
      h_new->Sumw2();
      h_new->Scale(1. / h_new->GetSum());

      TH1F * h_ref = NULL;
      if (qa_file_ref)
        {
          TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
              "h_QAG4Sim_CEMC_G4Hit_VSF", "TH1F");
          assert(h_ref);

          h_ref->Rebin(2);
          h_ref->Sumw2();
          h_ref->Scale(1. / h_ref->GetSum());
        }

      h_new->GetYaxis()->SetTitleOffset(1.5);
      h_new->GetYaxis()->SetTitle("Probability per bin");
      h_new->GetXaxis()->SetRangeUser(-0, .1);

      DrawReference(h_new, h_ref);
    }

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  //  p->SetLogz();

    {

      TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
          "h_QAG4Sim_CEMC_G4Hit_FractionEMVisibleEnergy", "TH1F");
      assert(h_new);

      h_new->Rebin(4);
      h_new->Sumw2();
      h_new->Scale(1. / h_new->GetSum());

      TH1F * h_ref = NULL;
      if (qa_file_ref)
        {
          TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
              "h_QAG4Sim_CEMC_G4Hit_FractionEMVisibleEnergy", "TH1F");
          assert(h_ref);

          h_ref->Rebin(4);
          h_ref->Scale(1. / h_ref->GetSum());
        }

      h_new->GetYaxis()->SetTitleOffset(1.5);
      h_new->GetYaxis()->SetTitle("Probability per bin");
//          h_new->GetXaxis()->SetRangeUser(-0, .1);

      DrawReference(h_new, h_ref);
    }

    PutInputFileName(c1, .04, qa_file_name_new, qa_file_name_ref);
  SaveCanvas(c1, TString(qa_file_name_new) + TString(c1->GetName()), true);
}

