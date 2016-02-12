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
QA_Draw_Jet_Spectrum(//
    const char * jet = "h_QAG4SimJet_AntiKt_Tower_r07",
//    const char * jet = "h_QAG4SimJet_AntiKt_Truth_r07",
    const char * qa_file_name_new =
        "data/G4sPHENIXCells_250jets25GeV.root_qa.root",
    const char * qa_file_name_ref =
        "data/G4sPHENIXCells_2000jets25GeV.root_qa.root")
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
          TString(jet) + TString("_Normalization"), "TH1D");
      assert(h_norm);

      Nevent_new = h_norm->GetBinContent(h_norm->GetXaxis()->FindBin("Event"));
    }
  if (qa_file_ref)
    {
      TH1D * h_norm = (TH1D *) qa_file_ref->GetObjectChecked(
          TString(jet) + TString("_Normalization"), "TH1D");
      assert(h_norm);

      Nevent_ref = h_norm->GetBinContent(h_norm->GetXaxis()->FindBin("Event"));
    }

  TCanvas *c1 = new TCanvas(TString("QA_Draw_Jet_Spectrum_") + TString(jet),
      TString("QA_Draw_Jet_Spectrum_") + TString(jet), 1800, 1000);
  c1->Divide(4, 2);
  int idx = 1;
  TPad * p;

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  p->SetLogy();

    {
      TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
          TString(jet) + TString("_Leading_eta"), "TH1F");
      assert(h_new);

      h_new -> Rebin(2);
      h_new->Sumw2();
      h_new->Scale(1. / Nevent_new);

      TH1F * h_ref = NULL;
      if (qa_file_ref)
        {
          TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
              TString(jet) + TString("_Leading_eta"), "TH1F");
          assert(h_ref);

          h_ref -> Rebin(2);
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
          TString(jet) + TString("_Leading_phi"), "TH1F");
      assert(h_new);

      h_new -> Rebin(2);
      h_new->Sumw2();
      h_new->Scale(1. / Nevent_new);

      TH1F * h_ref = NULL;
      if (qa_file_ref)
        {
          TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
              TString(jet) + TString("_Leading_phi"), "TH1F");
          assert(h_ref);

          h_ref -> Rebin(2);
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
          TString(jet) + TString("_Leading_Et"), "TH1F");
      assert(h_new);

      h_new->Sumw2();
      h_new->Scale(1. / Nevent_new);

      TH1F * h_ref = NULL;
      if (qa_file_ref)
        {
          TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
              TString(jet) + TString("_Leading_Et"), "TH1F");
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
          TString(jet) + TString("_Leading_Mass"), "TH1F");
      assert(h_new);

      h_new->Rebin(2);
      h_new->Sumw2();
      h_new->Scale(1. / Nevent_new);

      TH1F * h_ref = NULL;
      if (qa_file_ref)
        {
          TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
              TString(jet) + TString("_Leading_Mass"), "TH1F");
          assert(h_ref);

          h_ref->Rebin(2);
          h_ref->Scale(1. / Nevent_ref);
        }

      h_new->GetYaxis()->SetTitleOffset(1.5);
      h_new->GetYaxis()->SetTitle("Count / event / bin");
      //      h_new->GetXaxis()->SetRangeUser(-0, .1);

      DrawReference(h_new, h_ref);
    }

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  p->SetLogx();

    {
      TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
          TString(jet) + TString("_Leading_CompSize"), "TH1F");
      assert(h_new);

      h_new->Rebin(4);
      h_new->Sumw2();
      h_new->Scale(1. / Nevent_new);

      TH1F * h_ref = NULL;
      if (qa_file_ref)
        {
          TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
              TString(jet) + TString("_Leading_CompSize"), "TH1F");
          assert(h_ref);

          h_ref->Rebin(4);
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
          TString(jet) + TString("_Leading_CEMC_Ratio"), "TH1F");
      assert(h_new);

      h_new->Rebin(2);
      h_new->Sumw2();
      h_new->Scale(1. / Nevent_new);

      TH1F * h_ref = NULL;
      if (qa_file_ref)
        {
          TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
              TString(jet) + TString("_Leading_CEMC_Ratio"), "TH1F");
          assert(h_ref);

          h_ref->Rebin(2);
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
          TString(jet) + TString("_Leading_CEMC_HCalIN_Ratio"), "TH1F");
      assert(h_new);

      h_new->Rebin(2);
      h_new->Sumw2();
      h_new->Scale(1. / Nevent_new);

      TH1F * h_ref = NULL;
      if (qa_file_ref)
        {
          TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
              TString(jet) + TString("_Leading_CEMC_HCalIN_Ratio"), "TH1F");
          assert(h_ref);

          h_ref->Rebin(2);
          h_ref->Scale(1. / Nevent_ref);
        }

      h_new->GetYaxis()->SetTitleOffset(1.5);
      h_new->GetYaxis()->SetTitle("Count / event / bin");
      //      h_new->GetXaxis()->SetRangeUser(-0, .1);

      DrawReference(h_new, h_ref);
    }

    if (TString(jet).Contains("Truth"))
      {
        // truth jets

        p = (TPad *) c1->cd(idx++);
        c1->Update();
        p->SetLogy();

          {
            TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
                TString(jet) + TString("_Leading_Leakage_Ratio"), "TH1F");
            assert(h_new);

            h_new->Rebin(2);
            h_new->Sumw2();
            h_new->Scale(1. / Nevent_new);

            TH1F * h_ref = NULL;
            if (qa_file_ref)
              {
                TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
                    TString(jet) + TString("_Leading_Leakage_Ratio"), "TH1F");
                assert(h_ref);

                h_ref->Rebin(2);
                h_ref->Scale(1. / Nevent_ref);
              }

            h_new->GetYaxis()->SetTitleOffset(1.5);
            h_new->GetYaxis()->SetTitle("Count / event / bin");
            //      h_new->GetXaxis()->SetRangeUser(-0, .1);

            DrawReference(h_new, h_ref);
          }


      }

// inclusive jet stuff. Not very interesting.

//    p = (TPad *) c1->cd(idx++);
//    c1->Update();
//    p->SetLogy();
//
//      {
//        TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
//            TString(jet) + TString("_Inclusive_eta"), "TH1F");
//        assert(h_new);
//
////        h_new->Sumw2();
//        h_new->Scale(1. / Nevent_new);
//
//        TH1F * h_ref = NULL;
//        if (qa_file_ref)
//          {
//            TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
//                TString(jet) + TString("_Inclusive_eta"), "TH1F");
//            assert(h_ref);
//
//            h_ref->Scale(1. / Nevent_ref);
//          }
//
//        h_new->GetYaxis()->SetTitleOffset(1.5);
//        h_new->GetYaxis()->SetTitle("Energy (GeV) / event / bin");
//        //      h_new->GetXaxis()->SetRangeUser(-0, .1);
//
//        DrawReference(h_new, h_ref);
//      }
//
//    p = (TPad *) c1->cd(idx++);
//    c1->Update();
//    p->SetLogy();
//
//      {
//        TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
//            TString(jet) + TString("_Inclusive_phi"), "TH1F");
//        assert(h_new);
//
////        h_new->Sumw2();
//        h_new->Scale(1. / Nevent_new);
//
//        TH1F * h_ref = NULL;
//        if (qa_file_ref)
//          {
//            TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
//                TString(jet) + TString("_Inclusive_phi"), "TH1F");
//            assert(h_ref);
//
//            h_ref->Scale(1. / Nevent_ref);
//          }
//
//        h_new->GetYaxis()->SetTitleOffset(1.5);
//        h_new->GetYaxis()->SetTitle("Energy (GeV) / event / bin");
//        //      h_new->GetXaxis()->SetRangeUser(-0, .1);
//
//        DrawReference(h_new, h_ref);
//      }
//
//    p = (TPad *) c1->cd(idx++);
//    c1->Update();
//    p->SetLogy();
//    p->SetLogx();
//
//      {
//        TH1F * h_new = (TH1F *) qa_file_new->GetObjectChecked(
//            TString(jet) + TString("_Inclusive_E"), "TH1F");
//        assert(h_new);
//
//        h_new->Sumw2();
//        h_new->Scale(1. / Nevent_new);
//
//        TH1F * h_ref = NULL;
//        if (qa_file_ref)
//          {
//            TH1F * h_ref = (TH1F *) qa_file_ref->GetObjectChecked(
//                TString(jet) + TString("_Inclusive_E"), "TH1F");
//            assert(h_ref);
//
//            h_ref->Scale(1. / Nevent_ref);
//          }
//
//        h_new->GetYaxis()->SetTitleOffset(1.5);
//        h_new->GetYaxis()->SetTitle("Count / event / bin");
//        //      h_new->GetXaxis()->SetRangeUser(-0, .1);
//
//        DrawReference(h_new, h_ref);
//      }

  SaveCanvas(c1, TString(qa_file_name_new) + TString(c1->GetName()), true);
}

