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
using namespace std;

void
QA_Draw_Jet_TruthMatching(const char * jet =
    "h_QAG4SimJet_AntiKt_Truth_r07_AntiKt_Tower_r07",
    const char * qa_file_name_new = "G4sPHENIXCells_250jets25GeV.root_qa.root",
    const char * qa_file_name_ref = "G4sPHENIXCells_250jets25GeV.root_qa.root")
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

  const double Nevent_new = 250; // TODO: need to use normalization histos
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
//    TGraphErrors * ge = FitResolution(proj_new, false);

      proj_new->Draw("COLZ");
//    ge->Draw("p");

//    DrawReference(proj_new, proj_ref);
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
//    TGraphErrors * ge = FitResolution(proj_new, false);

      proj_new->Draw("COLZ");
//    ge->Draw("p");

//    DrawReference(proj_new, proj_ref);
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
//    TGraphErrors * ge = FitResolution(proj_new, false);

      proj_new->Draw("COLZ");
//    ge->Draw("p");

//    DrawReference(proj_new, proj_ref);
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
//    TGraphErrors * ge = FitResolution(proj_new, false);

      proj_new->Draw("COLZ");
//    ge->Draw("p");

//    DrawReference(proj_new, proj_ref);
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

void
DrawReference(TH1 * hnew, TH1 * href)
{

  hnew->SetLineColor(kBlue + 3);
  hnew->SetMarkerColor(kBlue + 3);
  hnew->SetLineWidth(2);
  hnew->SetMarkerStyle(kFullCircle);
  hnew->SetMarkerSize(1);

  if (href)
    {
      href->SetLineColor(kGreen + 1);
      href->SetFillColor(kGreen + 1);
      href->SetLineStyle(0);
      href->SetMarkerColor(kGreen + 1);
      href->SetLineWidth(0);
      href->SetMarkerStyle(kDot);
      href->SetMarkerSize(0);
    }

  hnew->Draw(); // set scale
  if (href)
    {
      href->Draw("HIST same");
      hnew->Draw("same"); // over lay data points
    }
}

TGraphErrors *
FitResolution(const TH2F * h2, const bool normalize_mean = true)
{

  TProfile * p2 = h2->ProfileX();

  int n = 0;
  double x[1000];
  double ex[1000];
  double y[1000];
  double ey[1000];

  for (int i = 1; i <= h2->GetNbinsX(); i++)
    {
      TH1D * h1 = h2->ProjectionY(Form("htmp_%d", rand()), i, i);

      if (h1->GetSum() < 10)
        continue;

      TF1 fgaus("fgaus", "gaus", -p2->GetBinError(i) * 4,
          p2->GetBinError(i) * 4);

      TF1 f2(Form("dgaus"), "gaus + [3]*exp(-0.5*((x-[1])/[4])**2) + [5]",
          -p2->GetBinError(i) * 4, p2->GetBinError(i) * 4);

      fgaus.SetParameter(1, p2->GetBinContent(i));
      fgaus.SetParameter(2, p2->GetBinError(i));

      h1->Fit(&fgaus, "MQ0");

//      f2.SetParameters(fgaus.GetParameter(0) / 2, fgaus.GetParameter(1),
//          fgaus.GetParameter(2), fgaus.GetParameter(0) / 2,
//          fgaus.GetParameter(2) / 4, 0);
//
//      h1->Fit(&f2, "MQ0");

      new TCanvas(Form("htmp_%d", rand()), Form("htmp_%d", rand()));
      h1->Draw();
      fgaus.Draw("same");
      break;

      x[n] = p2->GetBinCenter(i);
      ex[n] = (p2->GetBinCenter(2) - p2->GetBinCenter(1)) / 2;

      const double norm = normalize_mean ? fgaus.GetParameter(1) : 1;

      y[n] = fgaus.GetParameter(2) / norm;
      ey[n] = fgaus.GetParError(2) / norm;

      n++;
      delete h1;
    }

  TGraphErrors * ge = new TGraphErrors(n, x, y, 0, ey);
  ge->SetName(TString(h2->GetName()) + "_FitResolution");

  ge->SetLineColor(kBlue + 3);
  ge->SetMarkerColor(kBlue + 3);
  ge->SetLineWidth(2);
  ge->SetMarkerStyle(kFullCircle);
  ge->SetMarkerSize(1);
  return ge;
}
