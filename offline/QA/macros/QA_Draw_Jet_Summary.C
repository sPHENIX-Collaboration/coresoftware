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
#include <map>

//some common style files
#include "SaveCanvas.C"
#include "SetOKStyle.C"
#include "QA_Draw_Utility.C"
#include "QA_Draw_Jet_TruthMatching.C"
using namespace std;

void
QA_Draw_Jet_Summary(const char * jet_family = "AntiKt_Tower",
    const char * qa_file_name_new =
        "data/G4sPHENIXCells_2000jets25GeV.root_qa.root",
    const char * qa_file_name_ref =
        "data/G4sPHENIXCells_2000jets25GeV.root_qa.root")
{
  //! drawing energy range
  const double min_Et = 10;
  const double max_Et = 80;

  // style sets
  SetOKStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  TVirtualFitter::SetDefaultFitter("Minuit2");

  // file IO
  TFile * qa_file_new = new TFile(qa_file_name_new);
  assert(qa_file_new->IsOpen());

  // buffer for results
  vector<float> vec_radius;
  vector<TGraphErrors *> vec_phi_res;
  vector<TGraphErrors *> vec_eta_res;
  vector<TGraphErrors *> vec_e_res;
  vector<TGraphErrors *> vec_et_res;
  vector<TGraphErrors *> vec_reco_eff;
  vector<TGraphErrors *> vec_purity;

  // list and process all jets
  TList* hist_key_list = qa_file_new->GetListOfKeys();
  for (int i = 0; i < hist_key_list->GetSize(); ++i)
    {
      TString key_name = hist_key_list->At(i)->GetName();

      TString s_re_fullname = Form(
          "h_QAG4SimJet_.*_r[0-9]*_%s_r[0-9]*_Matching_Count_Truth_Et",
          jet_family); // regular expression for search
      TRegexp re_fullname(s_re_fullname, false);
      if (key_name.Index(re_fullname) == kNPOS)
        continue;

//      cout << " key_name = " << key_name << endl;
      TString jet_pair_name = key_name(0,
          key_name.Length() - TString("_Matching_Count_Truth_Et").Length()); // remove suffix

//      cout << " jet_pair_name = " << jet_pair_name << endl;

      //get jet radius
      TRegexp re_jetradius("_r[0-9]*", false);
      Ssiz_t index_radius = key_name.Index(re_jetradius); // first radius
      index_radius = key_name.Index(re_jetradius, index_radius + 1); // second radius
      assert(index_radius != kNPOS);
      float radius = 0;
      sscanf(key_name(index_radius, 100).Data(), "_r%f", &radius);
//      cout << " index_radius = " << index_radius << endl;
      assert(radius != 0);
      radius /= 10; // jet radius convention in DST names

      cout << "QA_Draw_Jet_Summary - process jet pair " << jet_pair_name
          << " with radius = " << radius << endl;

      vector<TGraphErrors *> resolution_efficiency_summary(
          QA_Draw_Jet_TruthMatching(jet_pair_name, qa_file_name_new,
              qa_file_name_ref));

      //save results
      vec_radius.push_back(radius);
      vec_phi_res.push_back(resolution_efficiency_summary[0]);
      vec_eta_res.push_back(resolution_efficiency_summary[1]);
      vec_e_res.push_back(resolution_efficiency_summary[2]);
      vec_et_res.push_back(resolution_efficiency_summary[3]);
      vec_reco_eff.push_back(resolution_efficiency_summary[4]);
      vec_purity.push_back(resolution_efficiency_summary[5]);

//      break;
    }

  // plot
  TCanvas *c1 = new TCanvas(
      TString("QA_Draw_Jet_Summary_") + TString(jet_family),
      TString("QA_Draw_Jet_Summary_") + TString(jet_family), 1800, 900);
  c1->Divide(3, 2);
  int idx = 1;
  TPad * p;

  // ------------------------------------
  p = (TPad *) c1->cd(idx++);
  c1->Update();
//  p->SetLogz();

  TH1 * h_frame =
      p->DrawFrame(min_Et, -.1, max_Et, .1,
          TString(jet_family)
              + " #phi Reconstruction;E_{T, Truth} (GeV);#phi_{Reco} - #phi_{Truth} (rad)");
  h_frame->GetYaxis()->SetTitleOffset(1.01);
  TLine * l = new TLine(min_Et, 0, max_Et, 0);
  l->Draw();
  p->SetGridx(0);
  p->SetGridy(0);
  TLegend * legend = new TLegend(0.7, 0.2, .95, 0.5);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(1001);
  legend->SetLineWidth(2);
  legend->SetLineColor(kBlack);
  legend->SetLineStyle(kSolid);
  for (int i = 0; i < vec_radius.size(); ++i)
    {
      const float radius = vec_radius[i];

      TGraphErrors* ge = vec_phi_res[i];
      assert(ge);
      ge = new TGraphErrors(*ge); // make a copy

      ge->SetLineColor(i + 2); // automatic color scheme from ROOT
      ge->SetMarkerColor(i + 2); // automatic color scheme from ROOT
      for (int idata = 0; idata < ge->GetN(); ++idata)
        {
          (ge->GetX())[idata] += i * 0.5; // shift x a little bit
          (ge->GetEX())[idata] = 0; // no x error bar
        }
      ge->Draw("p E l");
      legend->AddEntry(ge, Form("r = %.1f", radius), "elp");
    }
  legend->Draw();

  // ------------------------------------
  p = (TPad *) c1->cd(idx++);
  c1->Update();
//  p->SetLogz();

  TH1 * h_frame =
      p->DrawFrame(min_Et, -.1, max_Et, .1,
          TString(jet_family)
              + " #eta Reconstruction;E_{T, Truth} (GeV);#eta_{Reco} - #eta_{Truth}");
  h_frame->GetYaxis()->SetTitleOffset(1.01);
  TLine * l = new TLine(min_Et, 0, max_Et, 0);
  l->Draw();
  p->SetGridx(0);
  p->SetGridy(0);
  TLegend * legend = new TLegend(0.7, 0.2, .95, 0.5);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(1001);
  legend->SetLineWidth(2);
  legend->SetLineColor(kBlack);
  legend->SetLineStyle(kSolid);
  for (int i = 0; i < vec_radius.size(); ++i)
    {
      const float radius = vec_radius[i];

      TGraphErrors* ge = vec_eta_res[i];
      assert(ge);
      ge = new TGraphErrors(*ge); // make a copy

      ge->SetLineColor(i + 2); // automatic color scheme from ROOT
      ge->SetMarkerColor(i + 2); // automatic color scheme from ROOT
      for (int idata = 0; idata < ge->GetN(); ++idata)
        {
          (ge->GetX())[idata] += i * 0.5; // shift x a little bit
          (ge->GetEX())[idata] = 0; // no x error bar
        }
      ge->Draw("p E l");
      legend->AddEntry(ge, Form("r = %.1f", radius), "elp");
    }
  legend->Draw();

  // ------------------------------------
  p = (TPad *) c1->cd(idx++);
  c1->Update();
//  p->SetLogz();

  TH1 * h_frame = p->DrawFrame(min_Et, 0, max_Et, 2,
      TString(jet_family)
          + " Jet Energy Reconstruction;E_{Truth} (GeV);E_{Reco} / E_{Truth}");
  h_frame->GetYaxis()->SetTitleOffset(1.01);
  TLine * l = new TLine(min_Et, 1, max_Et, 1);
  l->Draw();
  p->SetGridx(0);
  p->SetGridy(0);
  TLegend * legend = new TLegend(0.7, 0.2, .95, 0.5);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(1001);
  legend->SetLineWidth(2);
  legend->SetLineColor(kBlack);
  legend->SetLineStyle(kSolid);
  for (int i = 0; i < vec_radius.size(); ++i)
    {
      const float radius = vec_radius[i];

      TGraphErrors* ge = vec_e_res[i];
      assert(ge);
      ge = new TGraphErrors(*ge); // make a copy

      ge->SetLineColor(i + 2); // automatic color scheme from ROOT
      ge->SetMarkerColor(i + 2); // automatic color scheme from ROOT
      for (int idata = 0; idata < ge->GetN(); ++idata)
        {
          (ge->GetX())[idata] += i * 0.5; // shift x a little bit
          (ge->GetEX())[idata] = 0; // no x error bar
        }
      ge->Draw("p E l");
      legend->AddEntry(ge, Form("r = %.1f", radius), "elp");
    }
  legend->Draw();

  // ------------------------------------
  p = (TPad *) c1->cd(idx++);
  c1->Update();
//  p->SetLogz();

  TH1 * h_frame =
      p->DrawFrame(min_Et, 0, max_Et, 2,
          TString(jet_family)
              + " Jet E_{T} Reconstruction;E_{T, Truth} (GeV);E_{T, Reco} / E_{T, Truth}");
  h_frame->GetYaxis()->SetTitleOffset(1.01);
  TLine * l = new TLine(min_Et, 1, max_Et, 1);
  l->Draw();
  p->SetGridx(0);
  p->SetGridy(0);
  TLegend * legend = new TLegend(0.7, 0.2, .95, 0.5);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(1001);
  legend->SetLineWidth(2);
  legend->SetLineColor(kBlack);
  legend->SetLineStyle(kSolid);
  for (int i = 0; i < vec_radius.size(); ++i)
    {
      const float radius = vec_radius[i];

      TGraphErrors* ge = vec_et_res[i];
      assert(ge);
      ge = new TGraphErrors(*ge); // make a copy

      ge->SetLineColor(i + 2); // automatic color scheme from ROOT
      ge->SetMarkerColor(i + 2); // automatic color scheme from ROOT
      for (int idata = 0; idata < ge->GetN(); ++idata)
        {
          (ge->GetX())[idata] += i * 0.5; // shift x a little bit
          (ge->GetEX())[idata] = 0; // no x error bar
        }
      ge->Draw("p E l");
      legend->AddEntry(ge, Form("r = %.1f", radius), "elp");
    }
  legend->Draw();

  // ------------------------------------
  p = (TPad *) c1->cd(idx++);
  c1->Update();
//  p->SetLogz();

  TH1 * h_frame = p->DrawFrame(min_Et, 0, max_Et, 1.2,
      TString(jet_family)
          + " Reco Efficiency;E_{T, Truth} (GeV);Reco efficiency");
  h_frame->GetYaxis()->SetTitleOffset(1.01);
  TLine * l = new TLine(min_Et, 1, max_Et, 1);
  l->Draw();
  p->SetGridx(0);
  p->SetGridy(0);
  TLegend * legend = new TLegend(0.7, 0.2, .95, 0.5);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(1001);
  legend->SetLineWidth(2);
  legend->SetLineColor(kBlack);
  legend->SetLineStyle(kSolid);
  for (int i = 0; i < vec_radius.size(); ++i)
    {
      const float radius = vec_radius[i];

      TGraphErrors* ge = vec_reco_eff[i];
      assert(ge);
      ge = new TGraphErrors(*ge); // make a copy

      ge->SetLineColor(i + 2); // automatic color scheme from ROOT
      ge->SetMarkerColor(i + 2); // automatic color scheme from ROOT
      for (int idata = 0; idata < ge->GetN(); ++idata)
        {
          (ge->GetX())[idata] += i * 0.5; // shift x a little bit
          (ge->GetEX())[idata] = 0; // no x error bar
        }
      ge->Draw("p E l");
      legend->AddEntry(ge, Form("r = %.1f", radius), "elp");
    }
  legend->Draw();

  // ------------------------------------
  p = (TPad *) c1->cd(idx++);
  c1->Update();
//  p->SetLogz();

  TH1 * h_frame = p->DrawFrame(min_Et, 0, max_Et, 1.2,
      TString(jet_family)
          + " Reconstruction Purity;E_{T, Reco} (GeV);Reconstruction Purity");
  h_frame->GetYaxis()->SetTitleOffset(1.01);
  TLine * l = new TLine(min_Et, 1, max_Et, 1);
  l->Draw();
  p->SetGridx(0);
  p->SetGridy(0);
  TLegend * legend = new TLegend(0.7, 0.2, .95, 0.5);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(1001);
  legend->SetLineWidth(2);
  legend->SetLineColor(kBlack);
  legend->SetLineStyle(kSolid);
  for (int i = 0; i < vec_radius.size(); ++i)
    {
      const float radius = vec_radius[i];

      TGraphErrors* ge = vec_purity[i];
      assert(ge);
      ge = new TGraphErrors(*ge); // make a copy

      ge->SetLineColor(i + 2); // automatic color scheme from ROOT
      ge->SetMarkerColor(i + 2); // automatic color scheme from ROOT
      for (int idata = 0; idata < ge->GetN(); ++idata)
        {
          (ge->GetX())[idata] += i * 0.5; // shift x a little bit
          (ge->GetEX())[idata] = 0; // no x error bar
        }
      ge->Draw("p E l");
      legend->AddEntry(ge, Form("r = %.1f", radius), "elp");
    }
  legend->Draw();

  SaveCanvas(c1, TString(qa_file_name_new) + TString(c1->GetName()), true);

}

