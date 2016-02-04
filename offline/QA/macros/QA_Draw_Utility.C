// $Id: $

/*!
 * \file QA_Draw_Utility.C
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

//! Draw 1D histogram along with its reference as shade
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

//! Fit for resolution of TH2F
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

//! Fit for profile along the Y direction of TH2F
TGraphErrors *
FitProfile(const TH2F * h2)
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

//      new TCanvas;
//      h1->Draw();
//      fgaus.Draw("same");
//      break;

      x[n] = p2->GetBinCenter(i);
      ex[n] = (p2->GetBinCenter(2) - p2->GetBinCenter(1)) / 2;
      y[n] = fgaus.GetParameter(1);
      ey[n] = fgaus.GetParameter(2);

//      p2->SetBinContent(i, fgaus.GetParameter(1));
//      p2->SetBinError(i, fgaus.GetParameter(2));

      n++;
      delete h1;
    }

  TGraphErrors * ge = new TGraphErrors(n, x, y, ex, ey);
  ge->SetName(TString(h2->GetName()) + "_FitProfile");
  ge->SetLineColor(kBlue + 3);
  ge->SetMarkerColor(kBlue + 3);
  ge->SetLineWidth(2);
  ge->SetMarkerStyle(kFullCircle);
  ge->SetMarkerSize(1);
  return ge;
}

//!ratio between two histograms with binominal error based on Wilson score interval. Assuming each histogram is count.
TH1 *
GetBinominalRatio(TH1 * h_pass, TH1 * h_n_trial)
{
  assert(h_pass);
  assert(h_n_trial);

  assert(h_pass->GetNbinsX() == h_n_trial->GetNbinsX());
  assert(h_pass->GetNbinsY() == h_n_trial->GetNbinsY());
  assert(h_pass->GetNbinsZ() == h_n_trial->GetNbinsZ());

  TH1 * h_ratio = (TH1 *) h_pass->Clone(TString(h_pass->GetName()) + "_Ratio");
  assert(h_ratio);
  h_ratio->Divide(h_n_trial); // a rough estimation first, also taking care of the overflow bins and zero bins

  for (int x = 1; x <= h_n_trial->GetNbinsX(); ++x)
    for (int y = 1; y <= h_n_trial->GetNbinsY(); ++y)
      for (int z = 1; z <= h_n_trial->GetNbinsZ(); ++z)
        {
          const double n_trial = h_n_trial->GetBinContent(x, y, z);

          if (n_trial > 0)
            {
              const double p = h_pass->GetBinContent(x, y, z) / n_trial;

              // Wilson score interval
              h_ratio->SetBinContent(x, y, z, //
                  (p + 1 / (2 * n_trial)) / (1 + 1 / n_trial));
              h_ratio->SetBinError(x, y,
                  z, //
                  TMath::Sqrt(
                      1. / n_trial * p * (1 - p) + 1. / (4 * n_trial * n_trial))
                      / (1 + 1 / n_trial));
            }
        }

  return h_ratio;
}
