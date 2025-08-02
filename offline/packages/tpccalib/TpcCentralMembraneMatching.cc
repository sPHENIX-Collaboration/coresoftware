
/**
 * \file TpcCentralMembraneMatching.cc
 * \brief match reconstructed CM clusters to CM pads, calculate differences, store on the node tree and compute distortion reconstruction maps
 * \author Tony Frawley <frawley@fsunuc.physics.fsu.edu>, Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcCentralMembraneMatching.h"

#include <trackbase/CMFlashDifferenceContainerv1.h>
#include <trackbase/CMFlashDifferencev1.h>
#include <trackbase/LaserClusterContainerv1.h>
#include <trackbase/LaserCluster.h>
#include <trackbase/TpcDefs.h>

#include <ffaobjects/EventHeader.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVector3.h>

#include <boost/format.hpp>

#include <cmath>
#include <iomanip>
#include <set>
#include <string>

int layerMins[3] = {16,23,39};
int layerMaxes[3] = {22, 38, 54};

namespace
{
  template <class T>
  inline constexpr T delta_phi(const T& phi)
  {
    if (phi > M_PI)
    {
      return phi - 2. * M_PI;
    }
    else if (phi <= -M_PI)
    {
      return phi + 2. * M_PI;
    }
    else
    {
      return phi;
    }
  }

  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }

  template <class T>
  inline T get_r(const T& x, const T& y)
  {
    return std::sqrt(square(x) + square(y));
  }

  // stream acts vector3
  [[maybe_unused]] std::ostream& operator<<(std::ostream& out, const Acts::Vector3& v)
  {
    out << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return out;
  }

  /// normalize distortions based on the number of entries in each cell, as recorded in the m_hentries histogram
  [[maybe_unused]] void normalize_distortions(TpcDistortionCorrectionContainer* dcc)
  {
    // loop over side
    for (unsigned int i = 0; i < 2; ++i)
    {
      // loop over bins in entries
      for (int ip = 0; ip < dcc->m_hentries[i]->GetNbinsX(); ++ip)
      {
        for (int ir = 0; ir < dcc->m_hentries[i]->GetNbinsY(); ++ir)
        {
          // count number of times a given cell was filled
          const auto entries = dcc->m_hentries[i]->GetBinContent(ip + 1, ir + 1);
          if (entries <= 1)
          {
            continue;
          }

          // normalize histograms
          for (const auto& h : {dcc->m_hDRint[i], dcc->m_hDPint[i], dcc->m_hDZint[i]})
          {
            h->SetBinContent(ip + 1, ir + 1, h->GetBinContent(ip + 1, ir + 1) / entries);
            h->SetBinError(ip + 1, ir + 1, h->GetBinError(ip + 1, ir + 1) / entries);
          }
        }
      }
    }
  }

  /// fill distortion correction histograms' guarding bins, to allow ::Interpolate to work over the full acceptance
  [[maybe_unused]] void fill_guarding_bins(TpcDistortionCorrectionContainer* dcc)
  {
    // loop over side
    for (unsigned int i = 0; i < 2; ++i)
    {
      for (const auto& h : {dcc->m_hDRint[i], dcc->m_hDPint[i], dcc->m_hDZint[i], dcc->m_hentries[i]})
      {
        // fill guarding phi bins
        /*
         * we use 2pi periodicity to do that:
         * - last valid bin is copied to first guarding bin;
         * - first valid bin is copied to last guarding bin
         */
        const auto phibins = h->GetNbinsX();
        const auto rbins = h->GetNbinsY();
        for (int ir = 0; ir < rbins; ++ir)
        {
          // copy last valid bin to first guarding bin
          h->SetBinContent(1, ir + 1, h->GetBinContent(phibins - 1, ir + 1));
          h->SetBinError(1, ir + 1, h->GetBinError(phibins - 1, ir + 1));

          // copy first valid bin to last guarding bin
          h->SetBinContent(phibins, ir + 1, h->GetBinContent(2, ir + 1));
          h->SetBinError(phibins, ir + 1, h->GetBinError(2, ir + 1));
        }

        // fill guarding r bins
        for (int iphi = 0; iphi < phibins; ++iphi)
        {
          // copy first valid bin to first guarding bin
          h->SetBinContent(iphi + 1, 1, h->GetBinContent(iphi + 1, 2));
          h->SetBinError(iphi + 1, 1, h->GetBinError(iphi + 1, 2));

          // copy last valid bin to last guarding bin
          h->SetBinContent(iphi + 1, rbins, h->GetBinContent(iphi + 1, rbins - 1));
          h->SetBinError(iphi + 1, rbins, h->GetBinError(iphi + 1, rbins - 1));
        }
      }
    }
  }

}  // namespace

//____________________________________________________________________________..
TpcCentralMembraneMatching::TpcCentralMembraneMatching(const std::string& name)
: SubsysReco(name)
{
  // calculate stripes center positions
  CalculateCenters(nPads_R1, R1_e, nGoodStripes_R1_e, keepUntil_R1_e, nStripesIn_R1_e, nStripesBefore_R1_e, cx1_e, cy1_e);
  CalculateCenters(nPads_R1, R1, nGoodStripes_R1, keepUntil_R1, nStripesIn_R1, nStripesBefore_R1, cx1, cy1);
  CalculateCenters(nPads_R2, R2, nGoodStripes_R2, keepUntil_R2, nStripesIn_R2, nStripesBefore_R2, cx2, cy2);
  CalculateCenters(nPads_R3, R3, nGoodStripes_R3, keepUntil_R3, nStripesIn_R3, nStripesBefore_R3, cx3, cy3);
}

//___________________________________________________________
void TpcCentralMembraneMatching::set_grid_dimensions(int phibins, int rbins)
{
  m_phibins = phibins;
  m_rbins = rbins;
}

// get the average phi rotation using smoothed histograms
double TpcCentralMembraneMatching::getPhiRotation_smoothed(TH1* hitHist, TH1* clustHist, bool side)
{
  // smooth the truth and cluster histograms
  hitHist->Smooth();
  clustHist->Smooth();

  // make a TF1 with a lambda function to make a function out of the truth histogram and shift it by a constant
  TF1* f1 = new TF1(
    "f1", [&](double* x, double* p)
    { return p[0] * hitHist->GetBinContent(hitHist->FindBin((x[0] - p[1]) > M_PI ? x[0] - p[1] - 2 * M_PI : x[0] - p[1])); },
    -M_PI, M_PI, 2);
  f1->SetParNames("A", "shift");
  f1->SetParameters(1.0, 0.0);
  f1->SetParLimits(1, -M_PI / 18, M_PI / 18);

  if (side && m_recoRotation[1][1] != -999)
  {
    f1->SetParameter(1, m_recoRotation[1][1]);
  }
  if (!side && m_recoRotation[0][1] != -999)
  {
    f1->SetParameter(1, m_recoRotation[0][1]);
  }

  f1->SetNpx(500);

  clustHist->Fit("f1", "IQ");

  //  clustHist->Draw();
  // f1->Draw("same");

  return f1->GetParameter(1);
}

void TpcCentralMembraneMatching::getRegionPhiRotation(bool side)
{
  TH1* oddHist[3];
  TH1* evenHist[3];

  TH1* oddTruthHist[3];
  TH1* evenTruthHist[3];

  for (int i = 0; i < (int) m_reco_RPeaks[side].size(); i++)
  {
    int peak = reco_r_phi[side]->GetYaxis()->FindBin(m_reco_RPeaks[side][i]);

    int region = 2;
    if (m_reco_RPeaks[side][i] < 41)
    {
      region = 0;
    }
    else if (m_reco_RPeaks[side][i] < 58)
    {
      region = 1;
    }

    if (m_reco_RMatches[side][i] % 2)
    {
      if (!oddHist[region])
      {
        oddHist[region] = (TH1*) reco_r_phi[side]->ProjectionX((boost::format("oddR%d_%d") %region %side).str().c_str(), peak - 2, peak + 2);
      }
      else
      {
        oddHist[region]->Add(reco_r_phi[side]->ProjectionX((boost::format("oddR%d_%d") %region %side).str().c_str(), peak - 2, peak + 2));
      }
    }
    else
    {
      if (!evenHist[region])
      {
        evenHist[region] = (TH1*) reco_r_phi[side]->ProjectionX((boost::format("evenR%d_%d") %region %side).str().c_str(), peak - 2, peak + 2);
      }
      else
      {
        evenHist[region]->Add(reco_r_phi[side]->ProjectionX((boost::format("evenR%d_%d") %region %side).str().c_str(), peak - 2, peak + 2));
      }
    }
  }

  for (int i = 0; i < (int) m_truth_RPeaks.size(); i++)
  {
    int peak = truth_r_phi[side]->GetYaxis()->FindBin(m_truth_RPeaks[i]);

    int region = 2;
    if (m_truth_RPeaks[i] < 41)
    {
      region = 0;
    }
    else if (m_truth_RPeaks[i] < 58)
    {
      region = 1;
    }

    if (i % 2)
    {
      if (!oddTruthHist[region])
      {
        oddTruthHist[region] = (TH1*) truth_r_phi[side]->ProjectionX((boost::format("oddTR%d_%d") %region %side).str().c_str(), peak - 1, peak + 1);
      }
      else
      {
        oddTruthHist[region]->Add(truth_r_phi[side]->ProjectionX((boost::format("oddTR%d_%d") %region %side).str().c_str(), peak - 1, peak + 1));
      }
    }
    else
    {
      if (!evenTruthHist[region])
      {
        evenTruthHist[region] = (TH1*) truth_r_phi[side]->ProjectionX((boost::format("evenTR%d_%d") %region %side).str().c_str(), peak - 1, peak + 1);
      }
      else
      {
        evenTruthHist[region]->Add(truth_r_phi[side]->ProjectionX((boost::format("evenTR%d_%d") %region %side).str().c_str(), peak - 1, peak + 1));
      }
    }
  }

  int regions[3]{1, 0, 2};

  double rotationEven[3]{-999, -999, -999};
  double rotationOdd[3]{-999, -999, -999};

  for (int region : regions)
  {
    TF1* fEven = new TF1(
      "fEven", [&](double* x, double* p)
      { return p[0] * evenTruthHist[region]->GetBinContent(evenTruthHist[region]->FindBin((x[0] - p[1]) > M_PI ? x[0] - p[1] - 2 * M_PI : x[0] - p[1])); },
      -M_PI, M_PI, 2);
    fEven->SetParNames("A", "shift");
    fEven->SetParameters(1.0, 0.0);
    fEven->SetParLimits(1, -M_PI / 18, M_PI / 18);
    fEven->SetNpx(1000);

    TF1* fOdd = new TF1(
      "fOdd", [&](double* x, double* p)
      { return p[0] * oddTruthHist[region]->GetBinContent(oddTruthHist[region]->FindBin((x[0] - p[1]) > M_PI ? x[0] - p[1] - 2 * M_PI : x[0] - p[1])); },
      -M_PI, M_PI, 2);
    fOdd->SetParNames("A", "shift");
    fOdd->SetParameters(1.0, 0.0);
    fOdd->SetParLimits(1, -M_PI / 18, M_PI / 18);
    fOdd->SetNpx(1000);

    if (region != 1)
    {
      if (rotationEven[1] != -999)
      {
        fEven->SetParameter(1, rotationEven[1]);
      }

      if (rotationOdd[1] != -999)
      {
        fOdd->SetParameter(1, rotationOdd[1]);
      }
    }

    evenHist[region]->Fit("fEven", "IQ");
    oddHist[region]->Fit("fOdd", "IQ");

    rotationEven[region] = fEven->GetParameter(1);
    rotationOdd[region] = fOdd->GetParameter(1);

    m_recoRotation[side][region] = (rotationOdd[region] + rotationEven[region]) / 2.0;
  }

  /*
    int regions[3] {1,0,2};

    for(int region : regions)
    {

      double evenMaxTruth = evenTruthHist[region]->GetMaximum();
      double oddMaxTruth = oddTruthHist[region]->GetMaximum();

      std::vector<double> evenTruthPeaks;
      std::vector<double> oddTruthPeaks;

      for(int i=2; i<oddTruthHist[region]->GetNbinsX(); i++)
      {

        if(evenTruthHist[region]->GetBinContent(i) > 0.15*evenMaxTruth)
        {
          evenTruthPeaks.push_back(evenTruthHist[region]->GetBinCenter(i));
        }

        if(oddTruthHist[region]->GetBinContent(i) > 0.15*oddMaxTruth)
        {
          oddTruthPeaks.push_back(oddTruthHist[region]->GetBinCenter(i));
        }

      }

      double evenMax = evenHist[region]->GetMaximum();
      double oddMax = oddHist[region]->GetMaximum();

      std::vector<double> evenRecoPeaks;
      std::vector<double> oddRecoPeaks;


    }

    double threshold = 0.02;
    std::vector<std::vector<int>> groupPhi;

    std::vector<double> finalEvenTruthPeaks;

    for (int i = 0; i < (int) evenTruthPeaks.size(); i++)
    {
      std::vector<int> tmpPhi;
      tmpPhi.push_back(i);

      bool closePeak = false;
      int currentPeak = -1;

      for (int j = 0; j < (int) finalEvenTruthPeaks.size(); j++)
      {
        for (int k = 0; k < (int) groupPhi[j].size(); k++)
        {
          if (fabs(evenTruthPeaks[i] - evenTruthPeaks[groupPhi[j][k]]) <= threshold || fabs(evenTruthPeaks[i] - finalEvenTruthPeaks[j]) <= threshold)
          {
            closePeak = true;
            currentPeak = j;
            break;
          }
        }
        if (closePeak)
        {
          break;
        }
      }

      if (!closePeak)
      {
        finalEvenTruthPeaks.push_back(evenTruthPeaks[i]);
        groupPhi.push_back(tmpPhi);
        tmpPhi.clear();
        continue;
      }

      groupPhi[currentPeak].push_back(i);
      double num = 0.0;
      double den = 0.0;
      for (int j : groupPhi[currentPeak])
      {
        double rHeight = evenTruthHist[region]->GetBinContent(evenTruthHist[region]->FindBin(evenTruthPeaks[j]));
        num += evenTruthPeaks[j] * rHeight;
        den += rHeight;
      }

      finalEvenTruthPeaks[currentPeak] = num / den;
    }

    std::vector<double> finalOddTruthPeaks;
    groupPhi.clear();

    for (int i = 0; i < (int) oddTruthPeaks.size(); i++)
    {
      std::vector<int> tmpPhi;
      tmpPhi.push_back(i);

      bool closePeak = false;
      int currentPeak = -1;

      for (int j = 0; j < (int) finalOddTruthPeaks.size(); j++)
      {
        for (int k = 0; k < (int) groupPhi[j].size(); k++)
        {
          if (fabs(oddTruthPeaks[i] - oddTruthPeaks[groupPhi[j][k]]) <= threshold || fabs(oddTruthPeaks[i] - finalOddTruthPeaks[j]) <= threshold)
          {
            closePeak = true;
            currentPeak = j;
            break;
          }
        }
        if (closePeak)
        {
          break;
        }
      }

      if (!closePeak)
      {
        finalOddTruthPeaks.push_back(oddTruthPeaks[i]);
        groupPhi.push_back(tmpPhi);
        tmpPhi.clear();
        continue;
      }

      groupPhi[currentPeak].push_back(i);
      double num = 0.0;
      double den = 0.0;
      for (int j : groupPhi[currentPeak])
      {
        double rHeight = oddTruthHist[region]->GetBinContent(oddTruthHist[region]->FindBin(oddTruthPeaks[j]));
        num += oddTruthPeaks[j] * rHeight;
        den += rHeight;
      }

      finalOddTruthPeaks[currentPeak] = num / den;
    }

    //reco
    std::vector<double> finalEvenRecoPeaks;
    groupPhi.clear();

    for (int i = 0; i < (int) evenRecoPeaks.size(); i++)
    {
      std::vector<int> tmpPhi;
      tmpPhi.push_back(i);

      bool closePeak = false;
      int currentPeak = -1;

      for (int j = 0; j < (int) finalEvenRecoPeaks.size(); j++)
      {
        for (int k = 0; k < (int) groupPhi[j].size(); k++)
        {
          if (fabs(evenRecoPeaks[i] - evenRecoPeaks[groupPhi[j][k]]) <= threshold || fabs(evenRecoPeaks[i] - finalEvenRecoPeaks[j]) <= threshold)
          {
            closePeak = true;
            currentPeak = j;
            break;
          }
        }
        if (closePeak)
        {
          break;
        }
      }

      if (!closePeak)
      {
        finalEvenRecoPeaks.push_back(evenRecoPeaks[i]);
        groupPhi.push_back(tmpPhi);
        tmpPhi.clear();
        continue;
      }

      groupPhi[currentPeak].push_back(i);
      double num = 0.0;
      double den = 0.0;
      for (int j : groupPhi[currentPeak])
      {
        double rHeight = evenHist[region]->GetBinContent(evenHist[region]->FindBin(evenRecoPeaks[j]));
        num += evenRecoPeaks[j] * rHeight;
        den += rHeight;
      }

      finalEvenRecoPeaks[currentPeak] = num / den;
    }

    std::vector<double> finalOddRecoPeaks;
    groupPhi.clear();

    for (int i = 0; i < (int) oddRecoPeaks.size(); i++)
    {
      std::vector<int> tmpPhi;
      tmpPhi.push_back(i);

      bool closePeak = false;
      int currentPeak = -1;

      for (int j = 0; j < (int) finalOddRecoPeaks.size(); j++)
      {
        for (int k = 0; k < (int) groupPhi[j].size(); k++)
        {
          if (fabs(oddRecoPeaks[i] - oddRecoPeaks[groupPhi[j][k]]) <= threshold || fabs(oddRecoPeaks[i] - finalOddRecoPeaks[j]) <= threshold)
          {
            closePeak = true;
            currentPeak = j;
            break;
          }
        }
        if (closePeak)
        {
          break;
        }
      }

      if (!closePeak)
      {
        finalOddRecoPeaks.push_back(oddRecoPeaks[i]);
        groupPhi.push_back(tmpPhi);
        tmpPhi.clear();
        continue;
      }

      groupPhi[currentPeak].push_back(i);
      double num = 0.0;
      double den = 0.0;
      for (int j : groupPhi[currentPeak])
      {
        double rHeight = oddHist[region]->GetBinContent(oddHist[region]->FindBin(oddRecoPeaks[j]));
        num += oddRecoPeaks[j] * rHeight;
        den += rHeight;
      }

      finalOddRecoPeaks[currentPeak] = num / den;
    }
    */
}

std::vector<int> TpcCentralMembraneMatching::doGlobalRMatching(TH2* r_phi, bool side)
{
  TH1D* proj = r_phi->ProjectionY("R_proj");

  if (side)
  {
    m_m[1] = 0.0;
    m_b[1] = 0.0;
  }
  else
  {
    m_m[0] = 0.0;
    m_b[0] = 0.0;
  }

  std::vector<double> rPeaks;
  std::vector<double> rHeights;
  std::vector<double> finalRPeaks;
  std::vector<double> finalRHeights;
  std::vector<std::vector<int>> groupR;

  proj->GetXaxis()->SetRangeUser(0, 41);
  double maxR1 = proj->GetMaximum();

  proj->GetXaxis()->SetRangeUser(41, 58);
  double maxR2 = proj->GetMaximum();

  proj->GetXaxis()->SetRangeUser(58, 100);
  double maxR3 = proj->GetMaximum();

  proj->GetXaxis()->SetRange(0, 0);

  std::vector<int> hitMatches;

  for (int i = 2; i < proj->GetNbinsX(); i++)
  {
    // peak is when content is higher than 0.15* maximum value and content is greater than or equal to both adjacent bins
    if ((proj->GetBinCenter(i) < 41.0 && proj->GetBinContent(i) > 0.15 * maxR1) || (proj->GetBinCenter(i) >= 41.0 && proj->GetBinCenter(i) < 58.0 && proj->GetBinContent(i) > 0.15 * maxR2) || (proj->GetBinCenter(i) >= 58.0 && proj->GetBinContent(i) > 0.15 * maxR3))
    {
      rPeaks.push_back(proj->GetBinCenter(i));
      rHeights.push_back(proj->GetBinContent(i));
    }
  }

  if (rPeaks.size() < 5)
  {
    if (side)
    {
      m_reco_RPeaks[1].clear();
    }
    else
    {
      m_reco_RPeaks[0].clear();
    }

    if (Verbosity())
    {
      std::cout << (side ? "side 1" : "side 0") << " has fewer than 5 radial peaks (only has " << rPeaks.size() << "). Returning empty vectors" << std::endl;
    }
    return hitMatches;
  }

  double threshold = 0.75;

  for (int i = 0; i < (int) rPeaks.size(); i++)
  {
    std::vector<int> tmpR;
    tmpR.push_back(i);

    bool closePeak = false;
    int currentPeak = -1;

    if (rPeaks[i] > 41.0)
    {
      threshold = 1.0;
    }

    for (int j = 0; j < (int) finalRPeaks.size(); j++)
    {
      for (int k = 0; k < (int) groupR[j].size(); k++)
      {
        if (fabs(rPeaks[i] - rPeaks[groupR[j][k]]) <= threshold || fabs(rPeaks[i] - finalRPeaks[j]) <= threshold)
        {
          closePeak = true;
          currentPeak = j;
          break;
        }
      }
      if (closePeak)
      {
        break;
      }
    }

    if (!closePeak)
    {
      finalRPeaks.push_back(rPeaks[i]);
      finalRHeights.push_back(rHeights[i]);
      groupR.push_back(tmpR);
      tmpR.clear();
      continue;
    }

    groupR[currentPeak].push_back(i);
    double num = 0.0;
    double den = 0.0;
    for (int j : groupR[currentPeak])
    {
      double rHeight = proj->GetBinContent(proj->FindBin(rPeaks[j]));
      num += rPeaks[j] * rHeight;
      den += rHeight;
    }

    finalRPeaks[currentPeak] = num / den;
    finalRHeights[currentPeak] = den;
  }

  if (side)
  {
    m_reco_RPeaks[1].clear();

    for (double& finalRPeak : finalRPeaks)
    {
      m_reco_RPeaks[1].push_back(finalRPeak);
    }
  }
  else
  {
    m_reco_RPeaks[0].clear();

    for (double& finalRPeak : finalRPeaks)
    {
      m_reco_RPeaks[0].push_back(finalRPeak);
    }
  }

  if (Verbosity())
  {
    std::cout << "finalRPeaks: {";
    for (int i = 0; i < (int) finalRPeaks.size() - 1; i++)
    {
      std::cout << finalRPeaks[i] << ", ";
    }
    std::cout << finalRPeaks[finalRPeaks.size() - 1] << "}" << std::endl;

    if (side)
    {
      std::cout << "m_reco_RPeaks[1]: {";
      for (int i = 0; i < (int) m_reco_RPeaks[1].size() - 1; i++)
      {
        std::cout << m_reco_RPeaks[1][i] << ", ";
      }
      std::cout << m_reco_RPeaks[1][m_reco_RPeaks[1].size() - 1] << "}" << std::endl;
    }
    else
    {
      std::cout << "m_reco_RPeaks[0]: {";
      for (int i = 0; i < (int) m_reco_RPeaks[0].size() - 1; i++)
      {
        std::cout << m_reco_RPeaks[0][i] << ", ";
      }
      std::cout << m_reco_RPeaks[0][m_reco_RPeaks[0].size() - 1] << "}" << std::endl;
    }

    std::cout << "rHeights: {";
    for (int i = 0; i < (int) finalRHeights.size() - 1; i++)
    {
      std::cout << finalRHeights[i] << ", ";
    }
    std::cout << finalRHeights[finalRHeights.size() - 1] << "}" << std::endl;
  }

  /*
  int middle_peak = -1;
  double closestPeak = 100000000.0;
  for (int i = 0; i < (int) finalRPeaks.size(); i++)
  {
    if (fabs(m_truth_RPeaks[21] - finalRPeaks[i]) < closestPeak)
    {
      closestPeak = fabs(m_truth_RPeaks[21] - finalRPeaks[i]);
      middle_peak = i;
    }
  }

  double middle_NN = 100000000.0;
  int middle_match = -1;
  for (int i = 0; i < (int) m_truth_RPeaks.size(); i++)
  {
    if (fabs(finalRPeaks[middle_peak] - m_truth_RPeaks[i]) < middle_NN)
    {
      middle_NN = fabs(finalRPeaks[middle_peak] - m_truth_RPeaks[i]);
      middle_match = i;
    }
  }

  */

  double first_NN = 100000000.0;
  int first_match = -1;
  for (int i = 0; i < (int) m_truth_RPeaks.size(); i++)
  {
    if (fabs(finalRPeaks[0] - m_truth_RPeaks[i]) < first_NN)
    {
      first_NN = fabs(finalRPeaks[0] - m_truth_RPeaks[i]);
      first_match = i;
    }
  }

  double last_NN = 100000000.0;
  int last_match = -1;
  for (int i = 0; i < (int) m_truth_RPeaks.size(); i++)
  {
    if (fabs(finalRPeaks[finalRPeaks.size() - 1] - m_truth_RPeaks[i]) < last_NN)
    {
      last_NN = fabs(finalRPeaks[finalRPeaks.size() - 1] - m_truth_RPeaks[i]);
      last_match = i;
    }
  }

  double bestSum = 100000000.0;
  int match_i = 0;
  int match_j = 0;

  int minI = -3;
  if (first_match + minI < 0)
  {
    minI = -first_match;
  }

  int maxJ = 3;
  if (last_match + maxJ >= (int) m_truth_RPeaks.size())
  {
    maxJ = ((int) m_truth_RPeaks.size()) - 1 - last_match;
  }

  std::vector<std::vector<double>> matches(3 - minI + 1, std::vector<double>(maxJ + 3 + 1, -999.0));
  std::vector<std::vector<std::vector<int>>> NNMatches(3 - minI + 1, std::vector<std::vector<int>>(maxJ + 3 + 1, std::vector<int>(finalRPeaks.size(), -1)));

  for (int i = minI; i <= 3; i++)
  {
    for (int j = -3; j <= maxJ; j++)
    {
      bool goodPair = true;

      if (m_fixShifts)
      {
        if (m_fieldOn)
        {
          if (i != -1 || j != -1)
          {
            goodPair = false;
          }
        }

        if (!m_fieldOn)
        {
          if (i != 0 || j != 0)
          {
            goodPair = false;
          }
        }
      }

      if (!goodPair)
      {
        continue;
      }

      double sum = 0.0;
      double m = (m_truth_RPeaks[last_match + j] - m_truth_RPeaks[first_match + i]) / (finalRPeaks[finalRPeaks.size() - 1] - finalRPeaks[0]);
      double b = ((m_truth_RPeaks[first_match + i] * finalRPeaks[finalRPeaks.size() - 1]) - (m_truth_RPeaks[last_match + j] * finalRPeaks[0])) / (finalRPeaks[finalRPeaks.size() - 1] - finalRPeaks[0]);
      // std::vector<int> tmpMatches;
      int recoIndex = 0;
      for (double finalRPeak : finalRPeaks)
      {
        int minMatch = 0;
        double minResidual = 1000000000.0;
        for (int k = 0; k < (int) m_truth_RPeaks.size(); k++)
        {
          double residual = fabs((m * finalRPeak + b) - m_truth_RPeaks[k]);
          if (residual < minResidual)
          {
            minResidual = residual;
            minMatch = k;
            if (minResidual < 0.5)
            {
              break;
            }
          }
        }
        sum += fabs((m * finalRPeak + b) - m_truth_RPeaks[minMatch]);
        // tmpMatches.push_back(minMatch);
        NNMatches[i - minI][j + 3][recoIndex] = minMatch;
        recoIndex++;
      }
      // NNMatches.push_back(tmpMatches);
      // matches.push_back(sum);
      matches[i - minI][j + 3] = sum;
      if (m_savehistograms)
      {
        if (side)
        {
          m_matchResiduals[1]->Fill(i, j, sum);
        }
        else
        {
          m_matchResiduals[0]->Fill(i, j, sum);
        }
      }
      if (sum < bestSum)
      {
        bestSum = sum;
        match_i = i - minI;
        match_j = j + 3;
      }
    }
  }

  if (Verbosity())
  {
    std::cout << "best total residual = " << bestSum << "   at i " << match_i + minI << "   j " << match_j - 3 << std::endl;
    for (int i = 0; i < (int) matches.size(); i++)
    {
      for (int j = 0; j < (int) matches[i].size(); j++)
      {
        if (matches[i][j] == -999.0)
        {
          continue;
        }
        std::cout << "total Residual match i=" << i + minI << "   j=" << j - 3 << "   : " << matches[i][j] << std::endl;
      }
    }
  }

  for (int& i : NNMatches[match_i][match_j])
  {
    hitMatches.push_back(i);
  }

  if (Verbosity())
  {
    for (int i = 0; i < (int) NNMatches.size(); i++)
    {
      for (int j = 0; j < (int) NNMatches[i].size(); j++)
      {
        if (matches[i][j] == -999.0)
        {
          continue;
        }
        double m = (m_truth_RPeaks[last_match + j - 3] - m_truth_RPeaks[first_match + i + minI]) / (finalRPeaks[finalRPeaks.size() - 1] - finalRPeaks[0]);
        double b = ((m_truth_RPeaks[first_match + i + minI] * finalRPeaks[finalRPeaks.size() - 1]) - (m_truth_RPeaks[last_match + j - 3] * finalRPeaks[0])) / (finalRPeaks[finalRPeaks.size() - 1] - finalRPeaks[0]);
        for (int k = 0; k < (int) NNMatches[i][j].size(); k++)
        {
          std::cout << "i " << i + minI << "   j " << j - 3 << "   Reco index " << k << "   recoR=" << finalRPeaks[k] << "   shifted R=" << (m * finalRPeaks[k] + b) << "   matchIndex=" << NNMatches[i][j][k] << "   truth R=" << m_truth_RPeaks[NNMatches[i][j][k]] << "   residual=" << (m * finalRPeaks[k] + b) - m_truth_RPeaks[NNMatches[i][j][k]] << std::endl;
        }
      }
    }
  }

  double final_m = (m_truth_RPeaks[last_match + match_j - 3] - m_truth_RPeaks[first_match + match_i + minI]) / (finalRPeaks[finalRPeaks.size() - 1] - finalRPeaks[0]);
  double final_b = ((m_truth_RPeaks[first_match + match_i + minI] * finalRPeaks[finalRPeaks.size() - 1]) - (m_truth_RPeaks[last_match + match_j - 3] * finalRPeaks[0])) / (finalRPeaks[finalRPeaks.size() - 1] - finalRPeaks[0]);

  if (side)
  {
    m_m[1] = final_m;
    m_b[1] = final_b;
    m_matchLow[1] = match_i + minI;
    m_matchHigh[1] = match_j - 3;
  }
  else
  {
    m_m[0] = final_m;
    m_b[0] = final_b;
    m_matchLow[0] = match_i + minI;
    m_matchHigh[0] = match_j - 3;
  }

  return hitMatches;
}

int TpcCentralMembraneMatching::getClusterRMatch(double clusterR, int side)
{
  double closestDist = 100.;
  int closestPeak = -1;

  // find cluster peak closest to position of passed cluster
  for (int j = 0; j < (int) m_reco_RPeaks[side].size(); j++)
  {
    if (std::abs(clusterR - m_reco_RPeaks[side][j]) < closestDist)
    {
      closestDist = std::abs(clusterR - m_reco_RPeaks[side][j]);
      closestPeak = j;
    }
  }

  // return hit match to cluster peak or -1 if closest peak failed (shouldn't be possible)
  if (closestPeak != -1)
  {
    return m_reco_RMatches[side][closestPeak];
  }
  else
  {
    return -1;
  }
}

//____________________________________________________________________________..
int TpcCentralMembraneMatching::InitRun(PHCompositeNode* topNode)
{
  if(!m_fieldOn)
  {
    m_useHeader = false;
  }

  if (m_savehistograms)
  {
    static constexpr float max_dr = 5.0;
    static constexpr float max_dphi = 0.05;

    fout = new TFile(m_histogramfilename.c_str(), "RECREATE");
    //fout.reset(new TFile(m_histogramfilename.c_str(), "RECREATE"));
    hxy_reco = new TH2F("hxy_reco", "reco cluster x:y", 800, -100, 100, 800, -80, 80);
    hxy_truth = new TH2F("hxy_truth", "truth cluster x:y", 800, -100, 100, 800, -80, 80);

    hdrdphi = new TH2F("hdrdphi", "dr vs dphi", 800, -max_dr, max_dr, 800, -max_dphi, max_dphi);
    hdrdphi->GetXaxis()->SetTitle("dr");
    hdrdphi->GetYaxis()->SetTitle("dphi");

    hrdr = new TH2F("hrdr", "dr vs r", 800, 0.0, 80.0, 800, -max_dr, max_dr);
    hrdr->GetXaxis()->SetTitle("r");
    hrdr->GetYaxis()->SetTitle("dr");

    hrdphi = new TH2F("hrdphi", "dphi vs r", 800, 0.0, 80.0, 800, -max_dphi, max_dphi);
    hrdphi->GetXaxis()->SetTitle("r");
    hrdphi->GetYaxis()->SetTitle("dphi");

    hdphi = new TH1F("hdphi", "dph", 800, -max_dphi, max_dphi);
    hdphi->GetXaxis()->SetTitle("dphi");

    hdr1_single = new TH1F("hdr1_single", "innner dr single", 200, -max_dr, max_dr);
    hdr2_single = new TH1F("hdr2_single", "mid dr single", 200, -max_dr, max_dr);
    hdr3_single = new TH1F("hdr3_single", "outer dr single", 200, -max_dr, max_dr);
    hdr1_double = new TH1F("hdr1_double", "innner dr double", 200, -max_dr, max_dr);
    hdr2_double = new TH1F("hdr2_double", "mid dr double", 200, -max_dr, max_dr);
    hdr3_double = new TH1F("hdr3_double", "outer dr double", 200, -max_dr, max_dr);
    hdrphi = new TH1F("hdrphi", "r * dphi", 200, -0.05, 0.05);
    hnclus = new TH1F("hnclus", " nclusters ", 3, 0., 3.);

    m_debugfile = new TFile(m_debugfilename.c_str(), "RECREATE");
    //m_debugfile.reset(new TFile(m_debugfilename.c_str(), "RECREATE"));
    match_tree = new TTree("match_tree", "Match TTree");

    match_tree->Branch("event", &m_event_index);
    match_tree->Branch("matched", &m_matched);
    match_tree->Branch("truthIndex", &m_truthIndex);
    match_tree->Branch("truthR", &m_truthR);
    match_tree->Branch("truthPhi", &m_truthPhi);
    match_tree->Branch("recoR", &m_recoR);
    match_tree->Branch("recoPhi", &m_recoPhi);
    match_tree->Branch("recoZ", &m_recoZ);
    match_tree->Branch("rawR", &m_rawR);
    match_tree->Branch("rawPhi", &m_rawPhi);
    match_tree->Branch("staticR", &m_staticR);
    match_tree->Branch("staticPhi", &m_staticPhi);
    match_tree->Branch("fitMode", &m_fitMode);
    match_tree->Branch("side", &m_side);
    match_tree->Branch("adc", &m_adc);
    match_tree->Branch("nhits", &m_nhits);
    match_tree->Branch("nLayers", &m_nLayers);
    match_tree->Branch("nIPhi", &m_nIPhi);
    match_tree->Branch("nIT", &m_nIT);
    match_tree->Branch("layerSD", &m_layersSD);
    match_tree->Branch("IPhiSD", &m_IPhiSD);
    match_tree->Branch("ITSD", &m_ITSD);
    match_tree->Branch("layerWeightedSD", &m_layersWeightedSD);
    match_tree->Branch("IPhiWeightedSD", &m_IPhiWeightedSD);
    match_tree->Branch("ITWeightedSD", &m_ITWeightedSD);
    match_tree->Branch("lowShift", &m_lowShift);
    match_tree->Branch("highShift", &m_highShift);
    match_tree->Branch("phiRotation", &m_phiRotation);
    match_tree->Branch("distanceToTruth", &m_distanceToTruth);
    match_tree->Branch("NNDistance", &m_NNDistance);
    match_tree->Branch("NNR", &m_NNR);
    match_tree->Branch("NNPhi", &m_NNPhi);
    match_tree->Branch("NNIndex", &m_NNIndex);

    // match_ntup = new TNtuple("match_ntup", "Match NTuple", "event:truthR:truthPhi:truthIndex:recoR:recoPhi:recoZ:side:adc:nhits:nLayers:nIPhi:nIT");
    m_matchResiduals[0] = new TH2F("matchResiduals_0", "Matching Residuals TPC South;Shift of smallest R from NN match;Shift of largest R from NN match", 7, -3.5, 3.5, 7, -3.5, 3.5);
    m_matchResiduals[1] = new TH2F("matchResiduals_1", "Matching Residuals TPC North;Shift of smallest R from NN match;Shift of largest R from NN match", 7, -3.5, 3.5, 7, -3.5, 3.5);
  }

  truth_r_phi[0] = new TH2F("truth_r_phi_0", "truth r vs #phi side 0;#phi (rad); r (cm)", 360, -M_PI, M_PI, 500, 0, 100);
  truth_r_phi[1] = new TH2F("truth_r_phi_1", "truth r vs #phi side 1;#phi (rad); r (cm)", 360, -M_PI, M_PI, 500, 0, 100);

  reco_r_phi[0] = new TH2F("reco_r_phi_0", "reco R vs #phi side 0;#phi (rad); r (cm)", 360, -M_PI, M_PI, 350, 20, 90);
  reco_r_phi[1] = new TH2F("reco_r_phi_1", "reco R vs #phi side 1;#phi (rad); r (cm)", 360, -M_PI, M_PI, 350, 20, 90);

  // Get truth cluster positions
  //=====================

  const double phi_petal = M_PI / 9.0;  // angle span of one petal

  /*
   * utility function to
   * - duplicate generated truth position to cover both sides of the central membrane
   * - assign proper z,
   * - insert in container
   */
  auto save_truth_position = [&](TVector3 source)
  {
    source.SetZ(-1);
    source.RotateZ(-M_PI / 18);
    m_truth_pos.push_back(source);
    truth_r_phi[0]->Fill(source.Phi(), source.Perp());

    source.SetZ(+1);
    source.RotateZ(M_PI / 18);
    source.SetX(-1*source.X());
    m_truth_pos.push_back(source);
    truth_r_phi[1]->Fill(source.Phi(), source.Perp());

  };

  // inner region extended is the 8 layers inside 30 cm
  for (int j = 0; j < nRadii; ++j)
  {
    for (int i = 0; i < nGoodStripes_R1_e[j]; ++i)
    {
      for (int k = 0; k < 18; ++k)
      {
        TVector3 dummyPos(cx1_e[i][j], cy1_e[i][j], 0.0);
        dummyPos.RotateZ(k * phi_petal);
        save_truth_position(dummyPos);
        //int truth_index_1 = k*10000 + j*100 + i;
        //int truth_index_1 = k*10000 + j*100 + (nGoodStripes_R1_e[j]-i-1);
        int truth_index_1;
        int truth_index_0;
        if(k<=8)
        {
         truth_index_0 = (26-k)*10000 + j*100 + (nGoodStripes_R1_e[j]-i-1);
	 truth_index_1 = (8-k)*10000 + j*100 + (nGoodStripes_R1_e[j]-i-1);
       }
	else
	{
	  truth_index_0 = (44-k)*10000 + j*100 + (nGoodStripes_R1_e[j]-i-1);
	  truth_index_1 = (26-k)*10000 + j*100 + (nGoodStripes_R1_e[j]-i-1);
	}
       m_truth_index.push_back(truth_index_0);
       m_truth_index.push_back(truth_index_1);

       if (Verbosity() > 2)
       {
        std::cout << " i " << i << " j " << j << " k " << k << " x1 " << dummyPos.X() << " y1 " << dummyPos.Y()
        << " theta " << std::atan2(dummyPos.Y(), dummyPos.X())
        << " radius " << get_r(dummyPos.X(), dummyPos.y()) << std::endl;
      }
      if (m_savehistograms)
      {
        hxy_truth->Fill(dummyPos.X(), dummyPos.Y());
      }
    }
  }
}

  // inner region is the 8 layers outside 30 cm
for (int j = 0; j < nRadii; ++j)
{
  for (int i = 0; i < nGoodStripes_R1[j]; ++i)
  {
    for (int k = 0; k < 18; ++k)
    {
      TVector3 dummyPos(cx1[i][j], cy1[i][j], 0.0);
      dummyPos.RotateZ(k * phi_petal);
      save_truth_position(dummyPos);
      //int truth_index_1 = k*10000 + (j+8)*100 + i;
      //int truth_index_1 = k*10000 + (j+8)*100 + (nGoodStripes_R1[j]-i-1);
      int truth_index_1;
      int truth_index_0;
      if(k<=8)
      {
        truth_index_0 = (26-k)*10000 + (j+8)*100 + (nGoodStripes_R1[j]-i-1);
        truth_index_1 = (8-k)*10000 + (j+8)*100 + (nGoodStripes_R1[j]-i-1);
      }
      else
      {
	truth_index_0 = (44-k)*10000 + (j+8)*100 + (nGoodStripes_R1[j]-i-1);
	truth_index_1 = (26-k)*10000 + (j+8)*100 + (nGoodStripes_R1[j]-i-1);
      }
      m_truth_index.push_back(truth_index_0);
      m_truth_index.push_back(truth_index_1);

      if (Verbosity() > 2)
      {
        std::cout << " i " << i << " j " << j << " k " << k << " x1 " << dummyPos.X() << " y1 " << dummyPos.Y()
        << " theta " << std::atan2(dummyPos.Y(), dummyPos.X())
        << " radius " << get_r(dummyPos.X(), dummyPos.y()) << std::endl;
      }
      if (m_savehistograms)
      {
        hxy_truth->Fill(dummyPos.X(), dummyPos.Y());
      }
    }
  }
}

for (int j = 0; j < nRadii; ++j)
{
  for (int i = 0; i < nGoodStripes_R2[j]; ++i)
  {
    for (int k = 0; k < 18; ++k)
    {
      TVector3 dummyPos(cx2[i][j], cy2[i][j], 0.0);
      dummyPos.RotateZ(k * phi_petal);
      save_truth_position(dummyPos);
      //int truth_index_1 = k*10000 + (j+16)*100 + i;
      //int truth_index_1 = k*10000 + (j+16)*100 + (nGoodStripes_R2[j]-i-1);
      int truth_index_1;
      int truth_index_0;
      if(k<=8)
      {
        truth_index_0 = (26-k)*10000 + (j+16)*100 + (nGoodStripes_R2[j]-i-1);
        truth_index_1 = (8-k)*10000 + (j+16)*100 + (nGoodStripes_R2[j]-i-1);
      }
      else
      {
	truth_index_0 = (44-k)*10000 + (j+16)*100 + (nGoodStripes_R2[j]-i-1);
	truth_index_1 = (26-k)*10000 + (j+16)*100 + (nGoodStripes_R2[j]-i-1);
      }
      m_truth_index.push_back(truth_index_0);
      m_truth_index.push_back(truth_index_1);

      if (Verbosity() > 2)
      {
        std::cout << " i " << i << " j " << j << " k " << k << " x1 " << dummyPos.X() << " y1 " << dummyPos.Y()
        << " theta " << std::atan2(dummyPos.Y(), dummyPos.X())
        << " radius " << get_r(dummyPos.X(), dummyPos.y()) << std::endl;
      }
      if (m_savehistograms)
      {
        hxy_truth->Fill(dummyPos.X(), dummyPos.Y());
      }
    }
  }
}

for (int j = 0; j < nRadii; ++j)
{
  for (int i = 0; i < nGoodStripes_R3[j]; ++i)
  {
    for (int k = 0; k < 18; ++k)
    {
      TVector3 dummyPos(cx3[i][j], cy3[i][j], 0.0);
      dummyPos.RotateZ(k * phi_petal);
      save_truth_position(dummyPos);
      //int truth_index_1 = k*10000 + (j+24)*100 + i;
      //int truth_index_1 = k*10000 + (j+24)*100 + (nGoodStripes_R3[j]-i-1);
      int truth_index_1;
      int truth_index_0;
      if(k<=8)
      {
        truth_index_0 = (26-k)*10000 + (j+24)*100 + (nGoodStripes_R3[j]-i-1);
        truth_index_1 = (8-k)*10000 + (j+24)*100 + (nGoodStripes_R3[j]-i-1);
      }
      else
      {
	truth_index_0 = (44-k)*10000 + (j+24)*100 + (nGoodStripes_R3[j]-i-1);
	truth_index_1 = (26-k)*10000 + (j+24)*100 + (nGoodStripes_R3[j]-i-1);
      }
      m_truth_index.push_back(truth_index_0);
      m_truth_index.push_back(truth_index_1);

      if (Verbosity() > 2)
      {
        std::cout << " i " << i << " j " << j << " k " << k << " x1 " << dummyPos.X() << " y1 " << dummyPos.Y()
        << " theta " << std::atan2(dummyPos.Y(), dummyPos.X())
        << " radius " << get_r(dummyPos.X(), dummyPos.y()) << std::endl;
      }
      if (m_savehistograms)
      {
        hxy_truth->Fill(dummyPos.X(), dummyPos.Y());
      }
    }
  }
}

int ret = GetNodes(topNode);
return ret;
}

//____________________________________________________________________________..
int TpcCentralMembraneMatching::process_event(PHCompositeNode* topNode)
{
  std::vector<TVector3> reco_pos;
  std::vector<TVector3> static_pos;
  std::vector<TVector3> raw_pos;
  std::vector<bool> reco_side;
  std::vector<bool> fitMode;
  std::vector<unsigned int> reco_nhits;
  std::vector<unsigned int> reco_adc;
  std::vector<unsigned int> reco_nLayers;
  std::vector<unsigned int> reco_nIPhi;
  std::vector<unsigned int> reco_nIT;
  std::vector<float> reco_SDLayer;
  std::vector<float> reco_SDIPhi;
  std::vector<float> reco_SDIT;
  std::vector<float> reco_SDWeightedLayer;
  std::vector<float> reco_SDWeightedIPhi;
  std::vector<float> reco_SDWeightedIT;

  eventHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if(!eventHeader)
  {
    std::cout << PHWHERE << " EventHeader Node missing, abort" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if(m_useHeader && eventHeader->get_EvtSequence() == 0)
  {
    m_useHeader = false;
  }


  if(m_useHeader)
  {
    m_event_index = eventHeader->get_EvtSequence();
  }

  if (!m_corrected_CMcluster_map || m_corrected_CMcluster_map->size() < 1000)
  {
    if(!m_useHeader)
    {
      m_event_index++;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }
  
  
  if (Verbosity())
  {
    std::cout << PHWHERE << "   working on event " << m_event_index << " which has " << m_corrected_CMcluster_map->size() << " clusters" << std::endl;
  }

  // reset output distortion correction container histograms
  for (const auto& harray : {m_dcc_out->m_hDRint, m_dcc_out->m_hDPint, m_dcc_out->m_hDZint, m_dcc_out->m_hentries})
  {
    for (const auto& h : harray)
    {
      h->Reset();
    }
  }

  reco_r_phi[0]->Reset();
  reco_r_phi[1]->Reset();

  int nClus_gtMin = 0;
  int clusterIndex = 0;

  // read the reconstructed CM clusters
  auto clusrange = m_corrected_CMcluster_map->getClusters();
  for (auto cmitr = clusrange.first;
   cmitr != clusrange.second;
   ++cmitr)
  {
    const auto& [cmkey, cmclus_orig] = *cmitr;
    // CMFlashCluster *cmclus = cmclus_orig;
    LaserCluster* cmclus = cmclus_orig;
    const unsigned int nhits = cmclus->getNhits();
    // const unsigned int nclus = cmclus->getNclusters();
    const unsigned int adc = cmclus->getAdc();
    bool side = (bool) TpcDefs::getSide(cmkey);


    // std::cout << "cluster " << clusterIndex << " side=" << side << "   sideKey=" << sideKey << "   are they the same? " << (side == sideKey ? "true" : "false") << std::endl;
    //  if(m_useOnly_nClus2 && nclus != 2) continue;


    nClus_gtMin++;

    // Do the static + average distortion corrections if the container was found
    //since incorrect z values are in cluster do to wrong t0 of laser flash, fixing based on the side for now
    //Acts::Vector3 pos(cmclus->getX(), cmclus->getY(), cmclus->getZ());
    Acts::Vector3 pos(cmclus->getX(), cmclus->getY(), (side ? 1.0 : -1.0));
    TVector3 tmp_raw(pos[0], pos[1], pos[2]);
    if (m_dcc_in_module_edge)
    {
      pos = m_distortionCorrection.get_corrected_position(pos, m_dcc_in_module_edge);
    }
    if (m_dcc_in_static)
    {
      pos = m_distortionCorrection.get_corrected_position(pos, m_dcc_in_static);
    }
    TVector3 tmp_static(pos[0], pos[1], pos[2]);

    if (m_dcc_in_average)
    {
      pos = m_distortionCorrection.get_corrected_position(pos, m_dcc_in_average);
    }
    TVector3 tmp_pos(pos[0], pos[1], pos[2]);


    // if(nclus == 1 && isRGap) continue;

    reco_pos.push_back(tmp_pos);
    static_pos.push_back(tmp_static);
    raw_pos.push_back(tmp_raw);
    reco_side.push_back(side);
    fitMode.push_back(cmclus->getFitMode());
    reco_nhits.push_back(nhits);
    reco_adc.push_back(adc);
    reco_nLayers.push_back(cmclus->getNLayers());
    reco_nIPhi.push_back(cmclus->getNIPhi());
    reco_nIT.push_back(cmclus->getNIT());
    reco_SDLayer.push_back(cmclus->getSDLayer());
    reco_SDIPhi.push_back(cmclus->getSDIPhi());
    reco_SDIT.push_back(cmclus->getSDIT());
    reco_SDWeightedLayer.push_back(cmclus->getSDWeightedLayer());
    reco_SDWeightedIPhi.push_back(cmclus->getSDWeightedIPhi());
    reco_SDWeightedIT.push_back(cmclus->getSDWeightedIT());

    if (side == 0)
    {
      reco_r_phi[0]->Fill(tmp_pos.Phi(), tmp_pos.Perp());
    }
    else
    {
      reco_r_phi[1]->Fill(tmp_pos.Phi(), tmp_pos.Perp());
    }

    if (Verbosity())
    {
      double raw_rad = sqrt(cmclus->getX() * cmclus->getX() + cmclus->getY() * cmclus->getY());
      double static_rad = sqrt(tmp_static.X() * tmp_static.X() + tmp_static.Y() * tmp_static.Y());
      double corr_rad = sqrt(tmp_pos.X() * tmp_pos.X() + tmp_pos.Y() * tmp_pos.Y());
      std::cout << "cluster " << clusterIndex << std::endl;
      clusterIndex++;
      std::cout << "found raw cluster " << cmkey << " side " << side << " with x " << cmclus->getX() << " y " << cmclus->getY() << " z " << cmclus->getZ() << " radius " << raw_rad << std::endl;
      std::cout << "        --- static corrected positions: " << tmp_static.X() << "  " << tmp_static.Y() << "  " << tmp_static.Z() << " radius " << static_rad << std::endl;
      std::cout << "                --- corrected positions: " << tmp_pos.X() << "  " << tmp_pos.Y() << "  " << tmp_pos.Z() << " radius " << corr_rad << std::endl;
    }

    if (m_savehistograms)
    {
      hxy_reco->Fill(tmp_pos.X(), tmp_pos.Y());
    }
  }

  if (nClus_gtMin < 25)
  {
    if(!m_useHeader)
    {
      m_event_index++;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  int truth_index = 0;
  int nMatched = 0;
  std::vector<bool> truth_matched(m_truth_pos.size(), false);
  std::vector<bool> reco_matched(reco_pos.size(), false);
  std::vector<int> truth_matchedRecoIndex(m_truth_pos.size(), -1);
  std::vector<int> reco_matchedTruthIndex(reco_pos.size(), -1);

  std::vector<float> reco_distToNN(reco_pos.size(), 10000.0);
  std::vector<float> NNDist(reco_pos.size(), 100000.0);
  std::vector<float> NNR(reco_pos.size(), 100000.0);
  std::vector<float> NNPhi(reco_pos.size(), 100000.0);
  std::vector<int> NNIndex(reco_pos.size(), -1);

  std::vector<std::vector<int>> truth_NNRecoIndex(m_truth_pos.size(), std::vector<int>(0));


  if(m_doFancy){

    // get global phi rotation for each module
    m_recoRotation[0][1] = getPhiRotation_smoothed(truth_r_phi[0]->ProjectionX("hR2", 206, 290), reco_r_phi[0]->ProjectionX("cR2_0", 206, 290), false);
    m_recoRotation[0][0] = getPhiRotation_smoothed(truth_r_phi[0]->ProjectionX("hR1", 151, 206), reco_r_phi[0]->ProjectionX("cR1_0", 151, 206), false);
    m_recoRotation[0][2] = getPhiRotation_smoothed(truth_r_phi[0]->ProjectionX("hR3", 290, 499), reco_r_phi[0]->ProjectionX("cR3_0", 290, 499), false);

    m_recoRotation[1][1] = getPhiRotation_smoothed(truth_r_phi[1]->ProjectionX("hR2", 206, 290), reco_r_phi[1]->ProjectionX("cR2_1", 206, 290), true);
    m_recoRotation[1][0] = getPhiRotation_smoothed(truth_r_phi[1]->ProjectionX("hR1", 151, 206), reco_r_phi[1]->ProjectionX("cR1_1", 151, 206), true);
    m_recoRotation[1][2] = getPhiRotation_smoothed(truth_r_phi[1]->ProjectionX("hR3", 290, 499), reco_r_phi[1]->ProjectionX("cR3_1", 290, 499), true);

    m_reco_RMatches[0] = doGlobalRMatching(reco_r_phi[0], false);
    m_reco_RMatches[1] = doGlobalRMatching(reco_r_phi[1], true);

    if (Verbosity())
    {
      for (int i = 0; i < (int) m_reco_RMatches[0].size(); i++)
      {
        std::cout << "side 0 cluster index " << i << "   hit match " << m_reco_RMatches[0][i] << "   recoPeak=" << m_reco_RPeaks[0][i] << "   shifted recoPeak=" << (m_m[0] * m_reco_RPeaks[0][i] + m_b[0]) << "   truthPeak=" << m_truth_RPeaks[m_reco_RMatches[0][i]] << "   residual=" << m_truth_RPeaks[m_reco_RMatches[0][i]] - (m_m[0] * m_reco_RPeaks[0][i] + m_b[0]) << std::endl;
      }

      for (int i = 0; i < (int) m_reco_RMatches[1].size(); i++)
      {
        std::cout << "side 1 cluster index " << i << "   hit match " << m_reco_RMatches[1][i] << "   recoPeak=" << m_reco_RPeaks[1][i] << "   shifted recoPeak=" << (m_m[1] * m_reco_RPeaks[1][i] + m_b[1]) << "   truthPeak=" << m_truth_RPeaks[m_reco_RMatches[1][i]] << "   residual=" << m_truth_RPeaks[m_reco_RMatches[1][i]] - (m_m[1] * m_reco_RPeaks[1][i] + m_b[1]) << std::endl;
      }
    }




    for (const auto& truth : m_truth_pos)
    {
      double tR = get_r(truth.X(), truth.Y());
      double tPhi = truth.Phi();
      double tZ = truth.Z();

      int truthRIndex = -1;

      // get which hit radial index this it
      for (int k = 0; k < (int) m_truth_RPeaks.size(); k++)
      {
        if (std::abs(tR - m_truth_RPeaks[k]) < 0.5)
        {
          truthRIndex = k;
          break;
        }
      }

      if (truthRIndex == -1)
      {
        truth_index++;
        continue;
      }

      double prev_dphi = 10000.0;

      int recoMatchIndex = -1;
      unsigned int reco_index = 0;
      for (const auto& reco : reco_pos)
      {
        if (reco_matched[reco_index] || reco_nhits[reco_index] < m_nHitsInCuster_minimum)
        {
          reco_index++;
          continue;
        }

        double rR = get_r(reco.X(), reco.Y());
        double rPhi = reco.Phi();
        bool side = reco_side[reco_index];

        int region = -1;

        if (rR < 41)
        {
          region = 0;
        }
        else if (rR >= 41 && rR < 58)
        {
          region = 1;
        }
        else if (rR >= 58)
        {
          region = 2;
        }

        if (region != -1)
        {
          if (side)
          {
            rPhi -= m_recoRotation[1][region];
          }
          else
          {
            rPhi -= m_recoRotation[0][region];
          }
        }

        int clustRMatchIndex = -1;
        clustRMatchIndex = getClusterRMatch(rR, (int) (side ? 1 : 0));

        if (clustRMatchIndex == -1 || truthRIndex != clustRMatchIndex)
        {
          reco_index++;
          continue;
        }

        if ((!side && tZ > 0) || (side && tZ < 0))
        {
          reco_index++;
          continue;
        }

        auto dphi = delta_phi(tPhi - rPhi);
        if (fabs(dphi) > m_phi_cut)
        {
          reco_index++;
          continue;
        }

        if (fabs(dphi) < fabs(prev_dphi))
        {
          prev_dphi = dphi;
          recoMatchIndex = reco_index;
          truth_matched[truth_index] = true;
        }

        reco_index++;
      }  // end loop over reco

      if (recoMatchIndex != -1)
      {
        truth_matchedRecoIndex[truth_index] = recoMatchIndex;
        reco_matched[recoMatchIndex] = true;
        reco_matchedTruthIndex[recoMatchIndex] = truth_index;
        nMatched++;

        if (Verbosity() > 2)
        {
          std::cout << "truth " << truth_index << " matched to reco " << recoMatchIndex << " and there are now " << nMatched << " matches" << std::endl;
          std::cout << "tR=" << std::setw(10) << tR << " tPhi=" << std::setw(10) << tPhi << " tZ=" << std::setw(10) << tZ << std::endl;
          std::cout << "rR=" << std::setw(10) << get_r(reco_pos[recoMatchIndex].X(), reco_pos[recoMatchIndex].Y()) << " rPhi=" << std::setw(10) << reco_pos[recoMatchIndex].Phi() << " rZ=" << std::setw(10) << reco_pos[recoMatchIndex].Z() << " rSide=" << reco_side[recoMatchIndex] << "   dPhi=" << prev_dphi << std::endl;
        }
      }

      truth_index++;
    }  // end loop over truth

    // loop again to find nearest neighbor for unmatched reco clusters

    int recoIndex = 0;
    for (const auto& reco : reco_pos)
    {
      double rR = get_r(reco.X(), reco.Y());
      double rPhi = reco.Phi();
      bool side = reco_side[recoIndex];

      int region = -1;

      if (rR < 41)
      {
        region = 0;
      }
      else if (rR >= 41 && rR < 58)
      {
        region = 1;
      }
      else if (rR >= 58)
      {
        region = 2;
      }

      if (region != -1)
      {
        if (side)
        {
          rPhi -= m_recoRotation[1][region];
        }
        else
        {
          rPhi -= m_recoRotation[0][region];
        }
      }

      int clustRMatchIndex = -1;
      clustRMatchIndex = getClusterRMatch(rR, (int) (side ? 1 : 0));


      int truthMatchIndex = -1;
      truth_index = 0;
      for (const auto& truth : m_truth_pos)
      {
        double tR = get_r(truth.X(), truth.Y());
        double tPhi = truth.Phi();
        double tZ = truth.Z();

	if( (!side && tZ > 0) || (side && tZ < 0) )
	{
	  truth_index++;
	  continue;
	}

        if(sqrt(pow(truth.X() - reco.X(),2) + pow(truth.Y() - reco.Y(),2)) < NNDist[recoIndex])
        {
         NNDist[recoIndex] = sqrt(pow(truth.X() - reco.X(),2) + pow(truth.Y() - reco.Y(),2));
         NNR[recoIndex] = tR;
         NNPhi[recoIndex] = tPhi;
         NNIndex[recoIndex] = m_truth_index[truth_index];
       }

       if (clustRMatchIndex == -1)
       {
         continue;
       }

       if (reco_matched[recoIndex])
       {
        if (reco_matchedTruthIndex[recoIndex] == truth_index)
        {
          reco_distToNN[recoIndex] = sqrt(pow(truth.X() - reco.X(), 2) + pow(truth.Y() - reco.Y(), 2));
          break;
        }
        else
        {
          truth_index++;
          continue;
        }
       }

      int truthRIndex = -1;

      // get which hit radial index this it
      for (int k = 0; k < (int) m_truth_RPeaks.size(); k++)
      {
        if (std::abs(tR - m_truth_RPeaks[k]) < 0.5)
        {
          truthRIndex = k;
          break;
        }
      }

      if (truthRIndex == -1 || truthRIndex != clustRMatchIndex)
      {
        truth_index++;
        continue;
      }

      if ((!side && tZ > 0) || (side && tZ < 0))
      {
        truth_index++;
        continue;
      }

      auto dphi = delta_phi(tPhi - rPhi);
      if (fabs(dphi) > m_phi_cut)
      {
        truth_index++;
        continue;
      }

      float dist = sqrt(pow(tR * sin(tPhi) - rR * sin(rPhi), 2) + pow(tR * cos(tPhi) - rR * cos(rPhi), 2));

      if (dist < reco_distToNN[recoIndex])
      {
        reco_distToNN[recoIndex] = dist;
        truthMatchIndex = truth_index;
      }

      truth_index++;

      }  // end of truth loop

      if (!reco_matched[recoIndex])
      {
        reco_matchedTruthIndex[recoIndex] = truthMatchIndex;
      }
      recoIndex++;
    }
  } //end fancy
  else
  {
    int reco_index = 0;
    for (const auto& reco : reco_pos)
    {
      double rR = get_r(reco.X(), reco.Y());
      double rPhi = reco.Phi();
      bool side = reco_side[reco_index];

      truth_index = 0;
      double minNNDist = 100000.0;
      int match_localTruth = -1;
      for (const auto& truth : m_truth_pos)
      {
	//std::cout << "reco_index " << reco_index << "   truth_index " << truth_index << std::endl;
	double tR = get_r(truth.X(), truth.Y());
        double tPhi = truth.Phi();
        double tZ = truth.Z();

        if ((!side && tZ > 0) || (side && tZ < 0))
        {
          truth_index++;
          continue;
        }

        auto dR = fabs(tR - rR);
        if (dR > 5.0)
        {
          truth_index++;
          continue;
        }

        auto dphi = delta_phi(tPhi - rPhi);
        if (fabs(dphi) > 0.05)
        {
          truth_index++;
          continue;
        }

        double dist = sqrt(pow(truth.X() - reco.X(),2) + pow(truth.Y() - reco.Y(),2));
        if (dist < minNNDist)
        {
          minNNDist = dist;
          match_localTruth = truth_index;
        }
	truth_index++;
      } // end truth loop

      if(match_localTruth == -1)
      {
	reco_index++;
	continue;
      }

      truth_NNRecoIndex[match_localTruth].push_back(reco_index);
      NNDist[reco_index] = sqrt(pow(m_truth_pos[match_localTruth].X() - reco.X(),2) + pow(m_truth_pos[match_localTruth].Y() - reco.Y(),2));
      NNR[reco_index] = get_r(m_truth_pos[match_localTruth].X(), m_truth_pos[match_localTruth].Y());
      NNPhi[reco_index] = m_truth_pos[match_localTruth].Phi();
      NNIndex[reco_index] = m_truth_index[match_localTruth];

      reco_index++;
    } // end reco loop

    truth_index = 0;
    for(const auto& truthIndex : truth_NNRecoIndex)
    {
      double maxADC = 0.0;
      int recoMatchIndex = -1;
      for(const auto& recoIndex : truthIndex)
      {

            if(reco_nLayers[recoIndex] < 2 || reco_nLayers[recoIndex] > 3) continue;

        int lowRMatch = 0;
        int highRMatch = (int)m_truth_RPeaks.size()-1;
        int upperBound = (int)m_truth_RPeaks.size();
        while(lowRMatch <= highRMatch)
        {
          int mid = lowRMatch + (highRMatch - lowRMatch)/2;
          if(m_truth_RPeaks[mid] >= reco_pos[recoIndex].Perp())
          {
            upperBound = mid;
            highRMatch = mid - 1;
          }
          else
          {
            lowRMatch = mid + 1;
          }
        }
        if(upperBound < 1) upperBound = 1;

	double dR = m_truth_pos[truth_index].Perp() - reco_pos[recoIndex].Perp();
	if(fabs(dR) > 0.5 * (m_truth_RPeaks[upperBound] - m_truth_RPeaks[upperBound-1]))
	{
	  continue;
	}

	double dphi = delta_phi(m_truth_pos[truth_index].Phi() - reco_pos[recoIndex].Phi());
	if(fabs(dphi) > 0.5*phiSpacing[upperBound-1])
	{
	  continue;
	}

        if(reco_adc[recoIndex] > maxADC)
        {
          maxADC = reco_adc[recoIndex];
          recoMatchIndex = recoIndex;
        }
      }

      if(recoMatchIndex >= 0)
      {
        truth_matchedRecoIndex[truth_index] = recoMatchIndex;
        reco_matched[recoMatchIndex] = true;
        reco_matchedTruthIndex[recoMatchIndex] = truth_index;
        truth_matched[truth_index] = true;
	nMatched++;
      }

      truth_index++;
    }
  } // end else for fancy

  // print some statistics:
  if (Verbosity())
  {
    const auto n_valid_truth = std::count_if(m_truth_pos.begin(), m_truth_pos.end(), [](const TVector3& pos)
     { return get_r(pos.x(), pos.y()) > 30; });
    std::cout << "TpcCentralMembraneMatching::process_event - m_truth_pos size: " << m_truth_pos.size() << std::endl;
    std::cout << "TpcCentralMembraneMatching::process_event - m_truth_pos size, r>30cm: " << n_valid_truth << std::endl;
    std::cout << "TpcCentralMembraneMatching::process_event - reco_pos size: " << reco_pos.size() << std::endl;
    std::cout << "TpcCentralMembraneMatching::process_event - matched_pair size: " << nMatched << std::endl;
  }

  if (m_savehistograms)
  {
    for (int i = 0; i < (int) reco_pos.size(); i++)
    {
      m_matched = reco_matched[i];
      m_truthR = m_truth_pos[reco_matchedTruthIndex[i]].Perp();
      m_truthPhi = m_truth_pos[reco_matchedTruthIndex[i]].Phi();
      m_distanceToTruth = sqrt(pow(m_truth_pos[reco_matchedTruthIndex[i]].Y() - reco_pos[i].Y(), 2) + pow(m_truth_pos[reco_matchedTruthIndex[i]].X() - reco_pos[i].X(), 2));
      m_truthIndex = m_truth_index[reco_matchedTruthIndex[i]];
      m_recoR = reco_pos[i].Perp();
      m_recoPhi = reco_pos[i].Phi();
      m_recoZ = reco_pos[i].Z();
      m_rawR = raw_pos[i].Perp();
      m_rawPhi = raw_pos[i].Phi();
      m_staticR = static_pos[i].Perp();
      m_staticPhi = static_pos[i].Phi();
      m_side = reco_side[i];
      m_fitMode = fitMode[i];
      m_adc = reco_adc[i];
      m_nhits = reco_nhits[i];
      m_nLayers = reco_nLayers[i];
      m_nIPhi = reco_nIPhi[i];
      m_nIT = reco_nIT[i];
      m_layersSD = reco_SDLayer[i];
      m_IPhiSD = reco_SDIPhi[i];
      m_ITSD = reco_SDIT[i];
      m_layersWeightedSD = reco_SDWeightedLayer[i];
      m_IPhiWeightedSD = reco_SDWeightedIPhi[i];
      m_ITWeightedSD = reco_SDWeightedIT[i];


      if (m_side)
      {
        m_lowShift = m_matchLow[1];
        m_highShift = m_matchHigh[1];
      }
      else
      {
        m_lowShift = m_matchLow[0];
        m_highShift = m_matchHigh[0];
      }

      int region = -1;

      if (m_recoR < 41)
      {
        region = 0;
      }
      else if (m_recoR >= 41 && m_recoR < 58)
      {
        region = 1;
      }
      else if (m_recoR >= 58)
      {
        region = 2;
      }

      if (region == -1)
      {
        m_phiRotation = -999.0;
      }
      else
      {
        m_phiRotation = m_recoRotation[m_side][region];
      }

      m_NNDistance = NNDist[i];
      m_NNR = NNR[i];
      m_NNPhi = NNPhi[i];
      m_NNIndex = NNIndex[i];

      match_tree->Fill();
    }
  }

  unsigned int ckey = 0;
  for (unsigned int i = 0; i < m_truth_pos.size(); i++)
  {
    if (!truth_matched[i])
    {
      continue;
    }

    //std::cout << "i=" << i << "   ckey=" << ckey << "   reco_index=" << truth_matchedRecoIndex[i] << std::endl;

    int reco_index = truth_matchedRecoIndex[i];

    if (reco_index == -1)
    {
      continue;
    }

    // std::cout << "good reco_index" << std::endl;

    auto cmdiff = new CMFlashDifferencev1();
    cmdiff->setTruthPhi(m_truth_pos[i].Phi());
    cmdiff->setTruthR(m_truth_pos[i].Perp());
    cmdiff->setTruthZ(m_truth_pos[i].Z());

    if (m_averageMode)
    {
      cmdiff->setRecoPhi(static_pos[reco_index].Phi());
      cmdiff->setRecoR(static_pos[reco_index].Perp());
      cmdiff->setRecoZ(static_pos[reco_index].Z());
      cmdiff->setNclusters(reco_nhits[reco_index]);
    }
    else
    {
      cmdiff->setRecoPhi(reco_pos[reco_index].Phi());
      cmdiff->setRecoR(reco_pos[reco_index].Perp());
      cmdiff->setRecoZ(reco_pos[reco_index].Z());
      cmdiff->setNclusters(reco_nhits[reco_index]);
    }

    if (Verbosity() > 1)
    {
      std::cout << "adding cmdiff to container with key " << ckey << " in event " << m_event_index << std::endl;
    }

    m_cm_flash_diffs->addDifferenceSpecifyKey(ckey, cmdiff);

    // store cluster position
    const double clus_r = reco_pos[reco_index].Perp();
    double clus_phi = reco_pos[reco_index].Phi();
    if (clus_phi < 0)
    {
      clus_phi += 2 * M_PI;
    }

    //const double clus_z = reco_pos[reco_index].z();
    const bool side = reco_side[reco_index];
    //const bool side = (clus_z < 0) ? false : true;
    // if(side != reco_side[reco_index]) std::cout << "sides do not match!" << std::endl;

    // calculate residuals (cluster - truth)
    double dr = reco_pos[reco_index].Perp() - m_truth_pos[i].Perp();
    double dphi = delta_phi(reco_pos[reco_index].Phi() - m_truth_pos[i].Phi());
    double rdphi = reco_pos[reco_index].Perp() * dphi;
    if (m_averageMode)
    {
      dr = static_pos[reco_index].Perp() - m_truth_pos[i].Perp();
      dphi = delta_phi(static_pos[reco_index].Phi() - m_truth_pos[i].Phi());
      rdphi = static_pos[reco_index].Perp() * dphi;
    }
    //currently, we cannot get any z distortion since we don't know when the laser actually flashed
    //so the distortion is set to 0 for now
    //const double dz = reco_pos[reco_index].z() - m_truth_pos[i].z();
    const double dz = 0.0;

    // fill distortion correction histograms
    /*
     * TODO:
     * - we might need to only fill the histograms for cm clusters that have 2 clusters only
     * - we might need a smoothing procedure to fill the bins that have no entries using neighbors
     */
    for (const auto& dcc : {m_dcc_out, m_dcc_out_aggregated.get()})
    {
      static_cast<TH2*>(dcc->m_hDRint[side])->Fill(clus_phi, clus_r, dr);
      static_cast<TH2*>(dcc->m_hDPint[side])->Fill(clus_phi, clus_r, rdphi);
      static_cast<TH2*>(dcc->m_hDZint[side])->Fill(clus_phi, clus_r, dz);
      static_cast<TH2*>(dcc->m_hentries[side])->Fill(clus_phi, clus_r);
    }

    ckey++;
  }

  if (Verbosity())
  {
    std::cout << "TpcCentralMembraneMatching::process_events - cmclusters: " << m_corrected_CMcluster_map->size() << std::endl;
    std::cout << "TpcCentralMembraneMatching::process_events - matched pairs: " << nMatched << std::endl;
    std::cout << "TpcCentralMembraneMatching::process_events - differences: " << m_cm_flash_diffs->size() << std::endl;
    std::cout << "TpcCentralMembraneMatching::process_events - entries: " << m_dcc_out->m_hentries[0]->GetEntries() << ", " << m_dcc_out->m_hentries[1]->GetEntries() << std::endl;
  }

  // normalize per-event distortion correction histograms and fill guarding bins
  normalize_distortions(m_dcc_out);
  fill_guarding_bins(m_dcc_out);

  if (Verbosity() > 2)
  {
    // read back differences from node tree as a check
    auto diffrange = m_cm_flash_diffs->getDifferences();
    for (auto cmitr = diffrange.first;
     cmitr != diffrange.second;
     ++cmitr)
    {
      auto key = cmitr->first;
      auto cmreco = cmitr->second;

      std::cout << " key " << key
      << " nclus " << cmreco->getNclusters()
      << " truth Phi " << cmreco->getTruthPhi() << " reco Phi " << cmreco->getRecoPhi()
      << " truth R " << cmreco->getTruthR() << " reco R " << cmreco->getRecoR()
      << " truth Z " << cmreco->getTruthZ() << " reco Z " << cmreco->getRecoZ()
      << std::endl;
    }
  }

  if(!m_useHeader){
    m_event_index++;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcCentralMembraneMatching::End(PHCompositeNode* /*topNode*/)
{
  // write distortion corrections
  if (m_dcc_out_aggregated)
  {
    // normalize aggregated distortion correction histograms and fill guarding bins
    if(m_doHadd)
    {
      normalize_distortions(m_dcc_out_aggregated.get());
    }
    fill_guarding_bins(m_dcc_out_aggregated.get());

    // create TFile and write all histograms
    std::unique_ptr<TFile> outputfile(TFile::Open(m_outputfile.c_str(), "RECREATE"));
    outputfile->cd();

    // loop over side
    for (unsigned int i = 0; i < 2; ++i)
    {
      for (const auto& h : {m_dcc_out_aggregated->m_hDRint[i], m_dcc_out_aggregated->m_hDPint[i], m_dcc_out_aggregated->m_hDZint[i], m_dcc_out_aggregated->m_hentries[i]})
      {
        if (h)
        {
          std::cout << "writing in " << m_outputfile << "   histogram " << h->GetName() << std::endl;
          h->Write();
        }
      }
    }
  }

  // write evaluation histograms
  if (m_savehistograms && fout)
  {
    fout->cd();

    hxy_reco->Write();
    hxy_truth->Write();
    hdrdphi->Write();
    hrdr->Write();
    hrdphi->Write();
    hdphi->Write();
    hdrphi->Write();
    hdr1_single->Write();
    hdr2_single->Write();
    hdr3_single->Write();
    hdr1_double->Write();
    hdr2_double->Write();
    hdr3_double->Write();
    hnclus->Write();

    fout->Close();
  }

  if (m_savehistograms && m_debugfile)
  {
    m_debugfile->cd();

    std::cout << "writing in " << m_debugfilename << "   match tree " << std::endl;

    /*
    match_tree->Write();
    // match_ntup->Write();
    std::cout << "writing in " << m_debugfilename << "   truth_r_phi " << std::endl;
    truth_r_phi[0]->Write();
    truth_r_phi[1]->Write();
    std::cout << "writing in " << m_debugfilename << "   reco_r_phi " << std::endl;
    reco_r_phi[0]->Write();
    reco_r_phi[1]->Write();
    std::cout << "writing in " << m_debugfilename << "   m_matchResiduals " << std::endl;
    m_matchResiduals[0]->Write();
    m_matchResiduals[1]->Write();
    */
    m_debugfile->Write();
    std::cout << "closing " << m_debugfilename << std::endl;
    m_debugfile->Close();
    std::cout << m_debugfilename << " has been closed" << std::endl;

  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..

int TpcCentralMembraneMatching::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  // m_corrected_CMcluster_map  = findNode::getClass<CMFlashClusterContainer>(topNode, "CORRECTED_CM_CLUSTER");
  m_corrected_CMcluster_map = findNode::getClass<LaserClusterContainer>(topNode, "LASER_CLUSTER");
  if (!m_corrected_CMcluster_map)
  {
    std::cout << PHWHERE << "CORRECTED_CM_CLUSTER Node missing, abort." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // input tpc distortion correction module edge
  m_dcc_in_module_edge = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerModuleEdge");
  if (m_dcc_in_module_edge)
  {
    std::cout << "TpcCentralMembraneMatching::GetNodes - found TPC distortion correction container module edge" << std::endl;
  }

  // input tpc distortion correction static
  m_dcc_in_static = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerStatic");
  if (m_dcc_in_static)
  {
    std::cout << "TpcCentralMembraneMatching::GetNodes - found TPC distortion correction container static" << std::endl;
  }

  // input tpc distortion correction average
  m_dcc_in_average = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerAverage");
  if (m_dcc_in_average)
  {
    std::cout << "TpcCentralMembraneMatching::GetNodes - found TPC distortion correction container average" << std::endl;
  }

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }



  // Looking for the RUN node
  // PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  // if (!runNode)
  //{
  // std::cout << PHWHERE << "RUN Node missing, doing nothing." << std::endl;
  // return Fun4AllReturnCodes::ABORTRUN;
  //}
  // PHNodeIterator runiter(runNode);
  // auto flashDiffContainer = findNode::getClass<CMFlashDifferenceContainerv1>(runNode, "CM_FLASH_DIFFERENCES");
  auto flashDiffContainer = findNode::getClass<CMFlashDifferenceContainerv1>(topNode, "CM_FLASH_DIFFERENCES");
  if (!flashDiffContainer)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode* DetNode =
    dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    flashDiffContainer = new CMFlashDifferenceContainerv1;
    PHIODataNode<PHObject>* CMFlashDifferenceNode =
    new PHIODataNode<PHObject>(flashDiffContainer, "CM_FLASH_DIFFERENCES", "PHObject");
    // runNode->addNode(CMFlashDifferenceNode);
    DetNode->addNode(CMFlashDifferenceNode);
  }

  // m_cm_flash_diffs = findNode::getClass<CMFlashDifferenceContainerv1>(runNode, "CM_FLASH_DIFFERENCES");
  m_cm_flash_diffs = findNode::getClass<CMFlashDifferenceContainerv1>(topNode, "CM_FLASH_DIFFERENCES");
  if (!m_cm_flash_diffs)
  {
    std::cout << PHWHERE << " ERROR: Can't find CM_FLASH_DIFFERENCES." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  /*
  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode =
    dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
  if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

  m_cm_flash_diffs = new CMFlashDifferenceContainerv1;
  PHIODataNode<PHObject> *CMFlashDifferenceNode =
    new PHIODataNode<PHObject>(m_cm_flash_diffs, "CM_FLASH_DIFFERENCES", "PHObject");
  DetNode->addNode(CMFlashDifferenceNode);
  */

  //// output tpc fluctuation distortion container
  //// this one is filled on the fly on a per-CM-event basis, and applied in the tracking chain
  const std::string dcc_out_node_name = "TpcDistortionCorrectionContainerFluctuation";
  m_dcc_out = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, dcc_out_node_name);
  if (!m_dcc_out)
  {
    /// Get the RUN node and check
    auto runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
    if (!runNode)
    {
      std::cout << "TpcCentralMembraneMatching::InitRun - RUN Node missing, quitting" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    std::cout << "TpcCentralMembraneMatching::GetNodes - creating TpcDistortionCorrectionContainer in node " << dcc_out_node_name << std::endl;
    m_dcc_out = new TpcDistortionCorrectionContainer;
    m_dcc_out->m_dimensions = 2;
    m_dcc_out->m_phi_hist_in_radians = false;
    m_dcc_out->m_interpolate_z = true;
    auto node = new PHDataNode<TpcDistortionCorrectionContainer>(m_dcc_out, dcc_out_node_name);
    runNode->addNode(node);
  }

  // create per event distortions. Do not put on the node tree
  // m_dcc_out = new TpcDistortionCorrectionContainer;

  // also prepare the local distortion container, used to aggregate multple events
  m_dcc_out_aggregated.reset(new TpcDistortionCorrectionContainer);

  // compute axis limits to include guarding bins, needed for TH2::Interpolate to work
  const float phiMin = m_phiMin - (m_phiMax - m_phiMin) / m_phibins;
  const float phiMax = m_phiMax + (m_phiMax - m_phiMin) / m_phibins;

  const float rMin = m_rMin - (m_rMax - m_rMin) / m_rbins;
  const float rMax = m_rMax + (m_rMax - m_rMin) / m_rbins;

  /*
  double r_bins_mm[69] = {217.83-2,217.83,
                       223.83, 229.83, 235.83, 241.83, 247.83, 253.83, 259.83, 265.83, 271.83, 277.83, 283.83, 289.83, 295.83, 301.83, 306.83,
                       311.05, 317.92, 323.31, 329.27, 334.63, 340.59, 345.95, 351.91, 357.27, 363.23, 368.59, 374.55, 379.91, 385.87, 391.23, 397.19, 402.49,
                       411.53, 421.70, 431.90, 442.11, 452.32, 462.52, 472.73, 482.94, 493.14, 503.35, 513.56, 523.76, 533.97, 544.18, 554.39, 564.59, 574.76,
                       583.67, 594.59, 605.57, 616.54, 627.51, 638.48, 649.45, 660.42, 671.39, 682.36, 693.33, 704.30, 715.27, 726.24, 737.21, 748.18, 759.11, 759.11+2};

  double r_bins[69];

  for(int i=0; i<69; i++){
    r_bins[i] = r_bins_mm[i]/10.0;
  }



  double phiBins[206] = { 0., 6.3083-2 * M_PI, 6.3401-2 * M_PI, 6.372-2 * M_PI, 6.4039-2 * M_PI, 6.4358-2 * M_PI, 6.4676-2 * M_PI, 6.4995-2 * M_PI, 6.5314-2 * M_PI,
                          0.2618, 0.2937, 0.3256, 0.3574, 0.3893, 0.4212, 0.453, 0.4849, 0.5168, 0.5487, 0.5805, 0.6124, 0.6443, 0.6762, 0.7081,
                          0.7399, 0.7718, 0.7854, 0.8173, 0.8491, 0.881, 0.9129, 0.9448, 0.9767, 1.0085, 1.0404, 1.0723, 1.1041, 1.136, 1.1679,
                          1.1998, 1.2317, 1.2635, 1.2954, 1.309, 1.3409, 1.3727, 1.4046, 1.4365, 1.4684, 1.5002, 1.5321, 1.564, 1.5959, 1.6277,
                          1.6596, 1.6915, 1.7234, 1.7552, 1.7871, 1.819, 1.8326, 1.8645, 1.8963, 1.9282, 1.9601, 1.992, 2.0238, 2.0557, 2.0876,
                          2.1195, 2.1513, 2.1832, 2.2151, 2.247, 2.2788, 2.3107, 2.3426, 2.3562, 2.3881, 2.42, 2.4518, 2.4837, 2.5156, 2.5474,
                          2.5793, 2.6112, 2.6431, 2.6749, 2.7068, 2.7387, 2.7706, 2.8024, 2.8343, 2.8662, 2.8798, 2.9117, 2.9436, 2.9754, 3.0073,
                          3.0392, 3.0711, 3.1029, 3.1348, 3.1667, 3.1986, 3.2304, 3.2623, 3.2942, 3.326, 3.3579, 3.3898, 3.4034, 3.4353, 3.4671,
                          3.499, 3.5309, 3.5628, 3.5946, 3.6265, 3.6584, 3.6903, 3.7221, 3.754, 3.7859, 3.8178, 3.8496, 3.8815, 3.9134, 3.927,
                          3.9589, 3.9907, 4.0226, 4.0545, 4.0864, 4.1182, 4.1501, 4.182, 4.2139, 4.2457, 4.2776, 4.3095, 4.3414, 4.3732, 4.4051,
                          4.437, 4.4506, 4.4825, 4.5143, 4.5462, 4.5781, 4.61, 4.6418, 4.6737, 4.7056, 4.7375, 4.7693, 4.8012, 4.8331, 4.865,
                          4.8968, 4.9287, 4.9606, 4.9742, 5.0061, 5.0379, 5.0698, 5.1017, 5.1336, 5.1654, 5.1973, 5.2292, 5.2611, 5.2929, 5.3248,
                          5.3567, 5.3886, 5.4204, 5.4523, 5.4842, 5.4978, 5.5297, 5.5615, 5.5934, 5.6253, 5.6572, 5.689, 5.7209, 5.7528, 5.7847,
                          5.8165, 5.8484, 5.8803, 5.9122, 5.944, 5.9759, 6.0078, 6.0214, 6.0533, 6.0851, 6.117, 6.1489, 6.1808, 6.2127, 6.2445,
                          6.2764, 2 * M_PI};
  */

  // reset all output distortion container so that they match the requested grid size
  const std::array<const std::string, 2> extension = {{"_negz", "_posz"}};
  for (const auto& dcc : {m_dcc_out, m_dcc_out_aggregated.get()})
  {
    // set dimensions to 2, since central membrane flashes only provide distortions at z = 0
    dcc->m_dimensions = 2;

    // create all histograms
    for (int i = 0; i < 2; ++i)
    {
      delete dcc->m_hDPint[i];
      dcc->m_hDPint[i] = new TH2F((boost::format("hIntDistortionP%s") % extension[i]).str().c_str(), (boost::format("hIntDistortionP%s") % extension[i]).str().c_str(), m_phibins + 2, phiMin, phiMax, m_rbins + 2, rMin, rMax);
      delete dcc->m_hDRint[i];
      dcc->m_hDRint[i] = new TH2F((boost::format("hIntDistortionR%s") % extension[i]).str().c_str(), (boost::format("hIntDistortionR%s") % extension[i]).str().c_str(), m_phibins + 2, phiMin, phiMax, m_rbins + 2, rMin, rMax);
      delete dcc->m_hDZint[i];
      dcc->m_hDZint[i] = new TH2F((boost::format("hIntDistortionZ%s") % extension[i]).str().c_str(), (boost::format("hIntDistortionZ%s") % extension[i]).str().c_str(), m_phibins + 2, phiMin, phiMax, m_rbins + 2, rMin, rMax);
      delete dcc->m_hentries[i];
      dcc->m_hentries[i] = new TH2I((boost::format("hEntries%s") % extension[i]).str().c_str(), (boost::format("hEntries%s") % extension[i]).str().c_str(), m_phibins + 2, phiMin, phiMax, m_rbins + 2, rMin, rMax);

      // delete dcc->m_hDPint[i]; dcc->m_hDPint[i] = new TH2F( (boost::format("hIntDistortionP%s") % extension[i]).str().c_str(), (boost::format("hIntDistortionP%s") % extension[i]).str().c_str(), 205, phiBins, 68, r_bins );
      // delete dcc->m_hDRint[i]; dcc->m_hDRint[i] = new TH2F( (boost::format("hIntDistortionR%s") % extension[i]).str().c_str(), (boost::format("hIntDistortionR%s") % extension[i]).str().c_str(), 205, phiBins, 68, r_bins );
      // delete dcc->m_hDZint[i]; dcc->m_hDZint[i] = new TH2F( (boost::format("hIntDistortionZ%s") % extension[i]).str().c_str(), (boost::format("hIntDistortionZ%s") % extension[i]).str().c_str(), 205, phiBins, 68, r_bins );
      // delete dcc->m_hentries[i]; dcc->m_hentries[i] = new TH2I( (boost::format("hEntries%s") % extension[i]).str().c_str(), (boost::format("hEntries%s") % extension[i]).str().c_str(), 205, phiBins, 68, r_bins);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________
void TpcCentralMembraneMatching::CalculateCenters(
  int nPads,
  const std::array<double, nRadii>& R,
  std::array<int, nRadii>& nGoodStripes,
  const std::array<int, nRadii>& keepUntil,
  std::array<int, nRadii>& nStripesIn,
  std::array<int, nRadii>& nStripesBefore,
  double cx[][nRadii], double cy[][nRadii])
{
  const double phi_module = M_PI / 6.0;  // angle span of a module
  const int pr_mult = 3;                 // multiples of intrinsic resolution of pads
  const int dw_mult = 8;                 // multiples of diffusion width
  const double diffwidth = 0.6 * mm;     // diffusion width
  const double adjust = 0.015;           // arbitrary angle to center the pattern in a petal

  double theta = 0.0;

  // center coords

  // calculate spacing first:
  std::array<double, nRadii> spacing{};
  for (int i = 0; i < nRadii; i++)
  {
    spacing[i] = 2.0 * ((dw_mult * diffwidth / R[i]) + (pr_mult * phi_module / nPads));
  }

  // center calculation
  for (int j = 0; j < nRadii; j++)
  {
    int i_out = 0;
    for (int i = keepThisAndAfter[j]; i < keepUntil[j]; i++)
    {
      if (j % 2 == 0)
      {
        theta = i * spacing[j] + (spacing[j] / 2) - adjust;
        cx[i_out][j] = R[j] * cos(theta) / cm;
        cy[i_out][j] = R[j] * sin(theta) / cm;
      }
      else
      {
        theta = (i + 1) * spacing[j] - adjust;
        cx[i_out][j] = R[j] * cos(theta) / cm;
        cy[i_out][j] = R[j] * sin(theta) / cm;
      }

      if (Verbosity() > 2)
      {
        std::cout << " j " << j << " i " << i << " i_out " << i_out << " theta " << theta << " cx " << cx[i_out][j] << " cy " << cy[i_out][j]
        << " radius " << sqrt(pow(cx[i_out][j], 2) + pow(cy[i_out][j], 2)) << std::endl;
      }

      i_out++;

      nStripesBefore_R1_e[0] = 0;

      nStripesIn[j] = keepUntil[j] - keepThisAndAfter[j];
      if (j == 0)
      {
        nStripesBefore[j] = 0;
      }
      else
      {
        nStripesBefore[j] = nStripesIn[j - 1] + nStripesBefore[j - 1];
      }
      nStripesBefore_R1_e[0] = 0;
    }
    nGoodStripes[j] = i_out;
  }
}
