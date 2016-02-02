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
QA_Draw_Jet_Spectrum(const char * jet =
    "h_QAG4SimJet_AntiKt_Truth_r07_AntiKt_Tower_r07",
    const char * qa_file_name_new = "G4sPHENIXCells_250jets25GeV.root_qa.root",
    const char * qa_file_name_ref = NULL)
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

  TCanvas *c1 = new TCanvas(TString("QA_Draw_Jet_Spectrum_") + TString(jet),
      TString("QA_Draw_Jet_Spectrum_") + TString(jet), 1800, 900);
  c1->Divide(4, 2);
  int idx = 1;
  TPad * p;

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
