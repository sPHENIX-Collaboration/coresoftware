/**
 * \file PHTpcCentralMembraneMatcher.cc
 * \brief match reconstructed CM clusters to CM pads, calculate differences, store on the node tree and compute distortion reconstruction maps
 * \author Tony Frawley <frawley@fsunuc.physics.fsu.edu>, Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHTpcCentralMembraneMatcher.h"

#include <trackbase/CMFlashClusterContainerv1.h>
#include <trackbase/CMFlashClusterv3.h>
#include <trackbase/CMFlashDifferenceContainerv1.h>
#include <trackbase/CMFlashDifferencev1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TString.h>
#include <TStyle.h>
#include <TVector3.h>

#include <boost/format.hpp>

#include <cmath>
#include <set>
#include <string>

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
PHTpcCentralMembraneMatcher::PHTpcCentralMembraneMatcher(const std::string& name)
  : SubsysReco(name)
{
  // calculate stripes center positions
  CalculateCenters(nPads_R1, R1_e, nGoodStripes_R1_e, keepUntil_R1_e, nStripesIn_R1_e, nStripesBefore_R1_e, cx1_e, cy1_e);
  CalculateCenters(nPads_R1, R1, nGoodStripes_R1, keepUntil_R1, nStripesIn_R1, nStripesBefore_R1, cx1, cy1);
  CalculateCenters(nPads_R2, R2, nGoodStripes_R2, keepUntil_R2, nStripesIn_R2, nStripesBefore_R2, cx2, cy2);
  CalculateCenters(nPads_R3, R3, nGoodStripes_R3, keepUntil_R3, nStripesIn_R3, nStripesBefore_R3, cx3, cy3);
}

//___________________________________________________________
void PHTpcCentralMembraneMatcher::set_grid_dimensions(int phibins, int rbins)
{
  m_phibins = phibins;
  m_rbins = rbins;
}

/*
std::vector<double> PHTpcCentralMembraneMatcher::getRGaps( TH2 *r_phi ){

  TH1D *proj = r_phi->ProjectionY("R_proj",1,360);

  std::vector<double> pass1;

  for(int i=2; i<proj->GetNbinsX(); i++){
    if(proj->GetBinContent(i) > 0.15*proj->GetMaximum() && proj->GetBinContent(i) >= proj->GetBinContent(i-1) && proj->GetBinContent(i) >= proj->GetBinContent(i+1)) pass1.push_back(proj->GetBinCenter(i));
  }

  for(int i=0; i<(int)pass1.size()-1; i++){
    if(pass1[i+1]-pass1[i] > 0.75) continue;

    if(proj->GetBinContent(proj->FindBin(pass1[i])) > proj->GetBinContent(proj->FindBin(pass1[i+1]))) pass1.erase(std::next(pass1.begin(), i+1));
    else pass1.erase(std::next(pass1.begin(), i));

    i--;
  }

  return pass1;

}
*/

// get the average phi rotation using smoothed histograms
double PHTpcCentralMembraneMatcher::getPhiRotation_smoothed(TH1* hitHist, TH1* clustHist)
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
  //  f1->SetParLimits(1,-M_PI/18,M_PI/18);

  f1->SetNpx(500);

  clustHist->Fit("f1", "IQ");

  //  clustHist->Draw();
  // f1->Draw("same");

  return f1->GetParameter(1);
}

// get vector with peak positions in Y (R) of histogram
std::vector<double> PHTpcCentralMembraneMatcher::getRPeaks(TH2* r_phi)
{
  TH1D* proj = r_phi->ProjectionY("R_proj");

  std::vector<double> rPeaks;

  for (int i = 2; i < proj->GetNbinsX(); i++)
  {
    // peak is when content is higher than 0.15* maximum value and content is greater than or equal to both adjacent bins
    if (proj->GetBinContent(i) > 0.15 * proj->GetMaximum() && proj->GetBinContent(i) >= proj->GetBinContent(i - 1) && proj->GetBinContent(i) >= proj->GetBinContent(i + 1))
    {
      rPeaks.push_back(proj->GetBinCenter(i));
    }
  }

  // if two peaks are within 0.75 cm of eachother, remove one with fewer counts
  for (int i = 0; i < (int) rPeaks.size() - 1; i++)
  {
    if (rPeaks[i + 1] - rPeaks[i] > 0.75)
    {
      continue;
    }
    if (proj->GetBinContent(proj->FindBin(rPeaks[i])) > proj->GetBinContent(proj->FindBin(rPeaks[i + 1])))
    {
      rPeaks.erase(rPeaks.begin() + i + 1);
    }
    else
    {
      rPeaks.erase(rPeaks.begin() + i);
    }
    i--;
  }
  return rPeaks;
}

int PHTpcCentralMembraneMatcher::getClusterRMatch(std::vector<int> hitMatches, std::vector<double> clusterPeaks, double clusterR)
{
  double closestDist = 100.;
  int closestPeak = -1;
  // find cluster peak closest to position of passed cluster
  for (int j = 0; j < (int) clusterPeaks.size(); j++)
  {
    if (std::abs(clusterR - clusterPeaks[j]) < closestDist)
    {
      closestDist = std::abs(clusterR - clusterPeaks[j]);
      closestPeak = j;
    }
  }

  // return hit match to cluster peak or -1 if closest peak failed (shouldn't be possible)
  if (closestPeak != -1)
  {
    return hitMatches[closestPeak];
  }
  else
  {
    return -1;
  }
}

//____________________________________________________________________________..
int PHTpcCentralMembraneMatcher::InitRun(PHCompositeNode* topNode)
{
  if (m_savehistograms)
  {
    static constexpr float max_dr = 5.0;
    static constexpr float max_dphi = 0.05;

    fout.reset(new TFile(m_histogramfilename.c_str(), "RECREATE"));
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

    fout2.reset(new TFile(m_histogramfilename2.c_str(), "RECREATE"));
    match_ntup = new TNtuple("match_ntup", "Match NTuple", "event:truthR:truthPhi:recoR:recoPhi:recoZ:nclus:r1:phi1:e1:layer1:r2:phi2:e2:layer2");
  }

  hit_r_phi = new TH2F("hit_r_phi", "hit r vs #phi;#phi (rad); r (cm)", 360, -M_PI, M_PI, 500, 0, 100);

  clust_r_phi_pos = new TH2F("clust_r_phi_pos", "clust R vs #phi Z>0;#phi (rad); r (cm)", 360, -M_PI, M_PI, 500, 0, 100);
  clust_r_phi_neg = new TH2F("clust_r_phi_neg", "clust R vs #phi Z<0;#phi (rad); r (cm)", 360, -M_PI, M_PI, 500, 0, 100);

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
    source.SetZ(+1);
    m_truth_pos.push_back(source);

    hit_r_phi->Fill(source.Phi(), source.Perp());

    source.SetZ(-1);
    m_truth_pos.push_back(source);

    hit_r_phi->Fill(source.Phi(), source.Perp());
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
int PHTpcCentralMembraneMatcher::process_event(PHCompositeNode* /*topNode*/)
{
  std::vector<TVector3> reco_pos;
  std::vector<TVector3> pos1;
  std::vector<TVector3> pos2;
  std::vector<unsigned int> reco_nclusters;
  std::vector<unsigned int> reco_adc;
  std::vector<unsigned int> adc1;
  std::vector<unsigned int> adc2;
  std::vector<unsigned int> layer1;
  std::vector<unsigned int> layer2;

  // reset output distortion correction container histograms
  for (const auto& harray : {m_dcc_out->m_hDRint, m_dcc_out->m_hDPint, m_dcc_out->m_hDZint, m_dcc_out->m_hentries})
  {
    clust_r_phi_pos->Reset();
    clust_r_phi_neg->Reset();

    if (!m_corrected_CMcluster_map || m_corrected_CMcluster_map->size() < 100)
    {
      m_event_index++;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    for (const auto& h : harray)
    {
      h->Reset();
    }
  }

  // read the reconstructed CM clusters
  auto clusrange = m_corrected_CMcluster_map->getClusters();
  for (auto cmitr = clusrange.first;
       cmitr != clusrange.second;
       ++cmitr)
  {
    const auto& [cmkey, cmclus_orig] = *cmitr;
    CMFlashCluster* cmclus = cmclus_orig;
    const unsigned int nclus = cmclus->getNclusters();
    const unsigned int adc = cmclus->getAdc();

    if (m_useOnly_nClus2 && nclus != 2)
    {
      continue;
    }

    const bool isRGap = cmclus->getIsRGap();

    // Do the static + average distortion corrections if the container was found
    Acts::Vector3 pos(cmclus->getX(), cmclus->getY(), cmclus->getZ());
    Acts::Vector3 apos1(cmclus->getX1(), cmclus->getY1(), cmclus->getZ1());
    Acts::Vector3 apos2(cmclus->getX2(), cmclus->getY2(), cmclus->getZ2());
    if (m_dcc_in_static)
    {
      pos = m_distortionCorrection.get_corrected_position(pos, m_dcc_in_static);
      apos1 = m_distortionCorrection.get_corrected_position(apos1, m_dcc_in_static);
      apos2 = m_distortionCorrection.get_corrected_position(apos2, m_dcc_in_static);
    }
    if (m_dcc_in_average)
    {
      pos = m_distortionCorrection.get_corrected_position(pos, m_dcc_in_average);
      apos1 = m_distortionCorrection.get_corrected_position(apos1, m_dcc_in_average);
      apos2 = m_distortionCorrection.get_corrected_position(apos2, m_dcc_in_average);
    }

    TVector3 tmp_pos(pos[0], pos[1], pos[2]);
    TVector3 tmp_pos1(apos1[0], apos1[1], apos1[2]);
    TVector3 tmp_pos2(apos2[0], apos2[1], apos2[2]);

    if (nclus == 1 && isRGap)
    {
      continue;
    }

    reco_pos.push_back(tmp_pos);
    pos1.push_back(tmp_pos1);
    pos2.push_back(tmp_pos2);
    reco_nclusters.push_back(nclus);
    reco_adc.push_back(adc);
    adc1.push_back(cmclus->getAdc1());
    adc2.push_back(cmclus->getAdc2());
    layer1.push_back(cmclus->getLayer1());
    layer2.push_back(cmclus->getLayer2());

    if (tmp_pos.Z() < 0)
    {
      clust_r_phi_neg->Fill(tmp_pos.Phi(), tmp_pos.Perp());
    }
    else
    {
      clust_r_phi_pos->Fill(tmp_pos.Phi(), tmp_pos.Perp());
    }

    if (Verbosity())
    {
      double raw_rad = sqrt(cmclus->getX() * cmclus->getX() + cmclus->getY() * cmclus->getY());
      double corr_rad = sqrt(tmp_pos.X() * tmp_pos.X() + tmp_pos.Y() * tmp_pos.Y());
      std::cout << "found raw cluster " << cmkey << " with x " << cmclus->getX() << " y " << cmclus->getY() << " z " << cmclus->getZ() << " radius " << raw_rad << std::endl;
      std::cout << "                --- corrected positions: " << tmp_pos.X() << "  " << tmp_pos.Y() << "  " << tmp_pos.Z() << " radius " << corr_rad << std::endl;
    }

    if (m_savehistograms)
    {
      hxy_reco->Fill(tmp_pos.X(), tmp_pos.Y());
    }
  }

  // get global phi rotation for each module
  m_clustRotation_pos[0] = getPhiRotation_smoothed(hit_r_phi->ProjectionX("hR1", 151, 206), clust_r_phi_pos->ProjectionX("cR1_pos", 151, 206));
  m_clustRotation_pos[1] = getPhiRotation_smoothed(hit_r_phi->ProjectionX("hR2", 206, 290), clust_r_phi_pos->ProjectionX("cR2_pos", 206, 290));
  m_clustRotation_pos[2] = getPhiRotation_smoothed(hit_r_phi->ProjectionX("hR3", 290, 499), clust_r_phi_pos->ProjectionX("cR3_pos", 290, 499));

  m_clustRotation_neg[0] = getPhiRotation_smoothed(hit_r_phi->ProjectionX("hR1", 151, 206), clust_r_phi_neg->ProjectionX("cR1_neg", 151, 206));
  m_clustRotation_neg[1] = getPhiRotation_smoothed(hit_r_phi->ProjectionX("hR2", 206, 290), clust_r_phi_neg->ProjectionX("cR2_neg", 206, 290));
  m_clustRotation_neg[2] = getPhiRotation_smoothed(hit_r_phi->ProjectionX("hR3", 290, 499), clust_r_phi_neg->ProjectionX("cR3_neg", 290, 499));

  // get hit and cluster peaks
  std::vector<double> hit_RPeaks = getRPeaks(hit_r_phi);
  std::vector<double> clust_RPeaks_pos = getRPeaks(clust_r_phi_pos);
  std::vector<double> clust_RPeaks_neg = getRPeaks(clust_r_phi_neg);

  // identify where gaps between module 1&2 and 2&3 are by finding largest distances between peaks
  std::vector<double> clust_RGaps_pos;
  int R12Gap_pos = -1;
  int R23Gap_pos = -1;
  clust_RGaps_pos.reserve((int) clust_RPeaks_pos.size() - 1);
  for (int i = 0; i < (int) clust_RPeaks_pos.size() - 1; i++)
  {
    clust_RGaps_pos.push_back(clust_RPeaks_pos[i + 1] - clust_RPeaks_pos[i]);
  }
  int tmpGap_pos = std::distance(clust_RGaps_pos.begin(), std::max_element(clust_RGaps_pos.begin(), clust_RGaps_pos.end()));
  if (tmpGap_pos > (int) clust_RGaps_pos.size() / 2)
  {
    R23Gap_pos = tmpGap_pos;
    R12Gap_pos = std::distance(clust_RGaps_pos.begin(), std::max_element(clust_RGaps_pos.begin(), clust_RGaps_pos.begin() + R23Gap_pos));
  }
  else
  {
    R12Gap_pos = tmpGap_pos;
    R23Gap_pos = std::distance(clust_RGaps_pos.begin(), std::max_element(clust_RGaps_pos.begin() + R12Gap_pos + 1, clust_RGaps_pos.end()));
  }

  std::vector<double> clust_RGaps_neg;
  int R12Gap_neg = -1;
  int R23Gap_neg = -1;
  clust_RGaps_neg.reserve((int) clust_RPeaks_neg.size() - 1);
  for (int i = 0; i < (int) clust_RPeaks_neg.size() - 1; i++)
  {
    clust_RGaps_neg.push_back(clust_RPeaks_neg[i + 1] - clust_RPeaks_neg[i]);
  }
  int tmpGap_neg = std::distance(clust_RGaps_neg.begin(), std::max_element(clust_RGaps_neg.begin(), clust_RGaps_neg.end()));
  if (tmpGap_neg > (int) clust_RGaps_neg.size() / 2)
  {
    R23Gap_neg = tmpGap_neg;
    R12Gap_neg = std::distance(clust_RGaps_neg.begin(), std::max_element(clust_RGaps_neg.begin(), clust_RGaps_neg.begin() + R23Gap_neg));
  }
  else
  {
    R12Gap_neg = tmpGap_neg;
    R23Gap_neg = std::distance(clust_RGaps_neg.begin(), std::max_element(clust_RGaps_neg.begin() + R12Gap_neg + 1, clust_RGaps_neg.end()));
  }

  std::vector<int> hitMatches_pos;
  // match cluster peaks to hit peaks using gap positions
  for (int i = 0; i < (int) clust_RPeaks_pos.size(); i++)
  {
    // Module 1
    if (i < (R12Gap_pos + 1))
    {
      // module 1-2 gap is between 15 & 16 in hitPeaks
      hitMatches_pos.push_back(15 + i - R12Gap_pos);
      // if multiple rows missing, offset for each one
      if (clust_RGaps_pos[R12Gap_pos] > 3.6)
      {
        hitMatches_pos[i] -= 1;
      }
      if (clust_RGaps_pos[R12Gap_pos] > 4.6)
      {
        hitMatches_pos[i] -= 1;
      }
      if (clust_RGaps_pos[R12Gap_pos] > 5.8)
      {
        hitMatches_pos[i] -= 1;
      }
    }
    // module 1-2 gap is between 15 & 16
    else if (i < (R23Gap_pos + 1) && i >= (R12Gap_pos + 1))
    {
      hitMatches_pos.push_back(15 + 1 + i - (R12Gap_pos + 1));
      // module 2-3 gap is between 22 & 23
    }
    else if (i >= R23Gap_pos + 1)
    {
      hitMatches_pos.push_back(23 + 1 + i - (R23Gap_pos + 1));
    }
  }

  std::vector<int> hitMatches_neg;
  for (int i = 0; i < (int) clust_RPeaks_neg.size(); i++)
  {
    if (i < (R12Gap_neg + 1))
    {
      hitMatches_neg.push_back(15 + i - R12Gap_neg);
      if (clust_RGaps_neg[R12Gap_neg] > 3.6)
      {
        hitMatches_neg[i] -= 1;
      }
      if (clust_RGaps_neg[R12Gap_neg] > 4.6)
      {
        hitMatches_neg[i] -= 1;
      }
      if (clust_RGaps_neg[R12Gap_neg] > 5.8)
      {
        hitMatches_neg[i] -= 1;
      }
    }
    else if (i < (R23Gap_neg + 1) && i >= (R12Gap_neg + 1))
    {
      hitMatches_neg.push_back(15 + 1 + i - (R12Gap_neg + 1));
    }
    else if (i >= R23Gap_neg + 1)
    {
      hitMatches_neg.push_back(23 + 1 + i - (R23Gap_neg + 1));
    }
  }

  // Match reco and truth positions
  // std::map<unsigned int, unsigned int> matched_pair;
  std::vector<std::pair<unsigned int, unsigned int>> matched_pair;
  std::vector<unsigned int> matched_nclus;

  std::vector<bool> hits_matched(m_truth_pos.size());
  std::vector<bool> clusts_matched(reco_pos.size());

  // do iterative matching m_nMatchIter times
  for (int matchIt = 0; matchIt < m_nMatchIter; matchIt++)
  {
    // loop over truth positions
    for (unsigned int i = 0; i < m_truth_pos.size(); ++i)
    {
      if (hits_matched[i])
      {
        continue;
      }

      const double z1 = m_truth_pos[i].Z();
      const double rad1 = get_r(m_truth_pos[i].X(), m_truth_pos[i].Y());
      const double phi1 = m_truth_pos[i].Phi();

      int hitRadIndex = -1;

      // get which hit radial index this it
      for (int k = 0; k < (int) hit_RPeaks.size(); k++)
      {
        if (std::abs(rad1 - hit_RPeaks[k]) < 0.5)
        {
          hitRadIndex = k;
          break;
        }
      }

      // for iterative looping: identify closest phi
      double prev_dphi = 1.1 * m_phi_cut;
      int matchJ = -1;

      // loop over cluster positions
      for (unsigned int j = 0; j < reco_pos.size(); ++j)
      {
        if (clusts_matched[j])
        {
          continue;
        }

        int angleR = -1;

        if (reco_pos[j].Perp() < 41)
        {
          angleR = 0;
        }
        else if (reco_pos[j].Perp() >= 41 && reco_pos[j].Perp() < 58)
        {
          angleR = 1;
        }
        if (reco_pos[j].Perp() >= 58)
        {
          angleR = 2;
        }

        // const auto& nclus = reco_nclusters[j];
        double phi2 = reco_pos[j].Phi();
        const double z2 = reco_pos[j].Z();
        const double rad2 = get_r(reco_pos[j].X(), reco_pos[j].Y());
        if (angleR != -1)
        {
          if (z2 > 0)
          {
            phi2 -= m_clustRotation_pos[angleR];
          }
          else
          {
            phi2 -= m_clustRotation_neg[angleR];
          }
        }

        int clustRMatchIndex = -1;
        if (z2 > 0)
        {
          clustRMatchIndex = getClusterRMatch(hitMatches_pos, clust_RPeaks_pos, rad2);
        }
        else
        {
          clustRMatchIndex = getClusterRMatch(hitMatches_neg, clust_RPeaks_neg, rad2);
        }

        if (clustRMatchIndex == -1)
        {
          continue;
        }

        // const double phi2 = reco_pos[j].Phi();

        // only match pairs that are on the same side of the TPC
        const bool accepted_z = ((z1 > 0) == (z2 > 0));
        if (!accepted_z)
        {
          continue;
        }

        const bool accepted_r = (hitRadIndex == clustRMatchIndex);

        const auto dphi = delta_phi(phi1 - phi2);
        const bool accepted_phi = std::abs(dphi) < m_phi_cut;

        if (!accepted_r || !accepted_phi)
        {
          continue;
        }

        if (fabs(dphi) < fabs(prev_dphi))
        {
          prev_dphi = dphi;
          matchJ = j;
          hits_matched[i] = true;
        }
      }  // end loop over reco_pos

      if (matchJ != -1)
      {
        clusts_matched[matchJ] = true;
        matched_pair.emplace_back(i, matchJ);
        matched_nclus.push_back(reco_nclusters[matchJ]);

        if (m_savehistograms)
        {
          const auto& nclus = reco_nclusters[matchJ];
          const double rad2 = get_r(reco_pos[matchJ].X(), reco_pos[matchJ].Y());
          const double phi2 = reco_pos[matchJ].Phi();

          const auto dr = rad1 - rad2;
          const auto dphi = delta_phi(phi1 - phi2);

          hnclus->Fill((float) nclus);

          double r = rad2;

          hdrphi->Fill(r * dphi);
          hdphi->Fill(dphi);
          hrdphi->Fill(r, dphi);
          hdrdphi->Fill(dr, dphi);
          hrdr->Fill(r, dr);
          if (nclus == 1)
          {
            if (r < 40.0)
            {
              hdr1_single->Fill(dr);
            }
            if (r >= 40.0 && r < 58.0)
            {
              hdr2_single->Fill(dr);
            }
            if (r >= 58.0)
            {
              hdr3_single->Fill(dr);
            }
          }
          else
          {
            if (r < 40.0)
            {
              hdr1_double->Fill(dr);
            }
            if (r >= 40.0 && r < 58.0)
            {
              hdr2_double->Fill(dr);
            }
            if (r >= 58.0)
            {
              hdr3_double->Fill(dr);
            }
          }
        }  // end save histos
      }    // end if match was found
    }      // end loop over truth pos
  }        // end loop over matching iterations

  // print some statistics:
  if (Verbosity())
  {
    const auto n_valid_truth = std::count_if(m_truth_pos.begin(), m_truth_pos.end(), [](const TVector3& pos)
                                             { return get_r(pos.x(), pos.y()) > 30; });
    const auto n_reco_size1 = std::count_if(reco_nclusters.begin(), reco_nclusters.end(), [](const unsigned int& value)
                                            { return value == 1; });
    const auto n_reco_size2 = std::count_if(reco_nclusters.begin(), reco_nclusters.end(), [](const unsigned int& value)
                                            { return value == 2; });
    std::cout << "PHTpcCentralMembraneMatcher::process_event - m_truth_pos size: " << m_truth_pos.size() << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_event - m_truth_pos size, r>30cm: " << n_valid_truth << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_event - reco_pos size: " << reco_pos.size() << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_event - reco_pos size (nclus==1): " << n_reco_size1 << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_event - reco_pos size (nclus==2): " << n_reco_size2 << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_event - matched_pair size: " << matched_pair.size() << std::endl;
  }

  for (unsigned int ip = 0; ip < matched_pair.size(); ++ip)
  {
    const std::pair<unsigned int, unsigned int>& p = matched_pair[ip];
    const unsigned int& nclus = matched_nclus[ip];

    // add to node tree
    unsigned int key = p.first;
    auto cmdiff = new CMFlashDifferencev1();
    cmdiff->setTruthPhi(m_truth_pos[p.first].Phi());
    cmdiff->setTruthR(m_truth_pos[p.first].Perp());
    cmdiff->setTruthZ(m_truth_pos[p.first].Z());

    cmdiff->setRecoPhi(reco_pos[p.second].Phi());
    cmdiff->setRecoR(reco_pos[p.second].Perp());
    cmdiff->setRecoZ(reco_pos[p.second].Z());

    cmdiff->setNclusters(nclus);

    m_cm_flash_diffs->addDifferenceSpecifyKey(key, cmdiff);

    if (m_savehistograms)
    {
      match_ntup->Fill(m_event_index, m_truth_pos[p.first].Perp(), m_truth_pos[p.first].Phi(), reco_pos[p.second].Perp(), reco_pos[p.second].Phi(), reco_pos[p.second].Z(), nclus, pos1[p.second].Perp(), pos1[p.second].Phi(), adc1[p.second], layer1[p.second], pos2[p.second].Perp(), pos2[p.second].Phi(), adc2[p.second], layer2[p.second]);
    }

    // store cluster position
    const double clus_r = reco_pos[p.second].Perp();
    double clus_phi = reco_pos[p.second].Phi();
    if (clus_phi < 0)
    {
      clus_phi += 2 * M_PI;
    }

    const double clus_z = reco_pos[p.second].z();
    const unsigned int side = (clus_z < 0) ? 0 : 1;

    // calculate residuals (cluster - truth)
    const double dr = reco_pos[p.second].Perp() - m_truth_pos[p.first].Perp();
    const double dphi = delta_phi(reco_pos[p.second].Phi() - m_truth_pos[p.first].Phi());
    const double rdphi = reco_pos[p.second].Perp() * dphi;
    const double dz = reco_pos[p.second].z() - m_truth_pos[p.first].z();

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
  }

  if (Verbosity())
  {
    std::cout << "PHTpcCentralMembraneMatcher::process_events - cmclusters: " << m_corrected_CMcluster_map->size() << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_events - matched pairs: " << matched_pair.size() << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_events - differences: " << m_cm_flash_diffs->size() << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_events - entries: " << m_dcc_out->m_hentries[0]->GetEntries() << ", " << m_dcc_out->m_hentries[1]->GetEntries() << std::endl;
  }

  // normalize per-event distortion correction histograms and fill guarding bins
  normalize_distortions(m_dcc_out);
  fill_guarding_bins(m_dcc_out);

  if (Verbosity())
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

  m_event_index++;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTpcCentralMembraneMatcher::End(PHCompositeNode* /*topNode*/)
{
  // write distortion corrections
  if (m_dcc_out_aggregated)
  {
    // normalize aggregated distortion correction histograms and fill guarding bins
    normalize_distortions(m_dcc_out_aggregated.get());
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

  if (m_savehistograms && fout2)
  {
    fout2->cd();

    match_ntup->Write();
    hit_r_phi->Write();
    clust_r_phi_pos->Write();
    clust_r_phi_neg->Write();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..

int PHTpcCentralMembraneMatcher::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  m_corrected_CMcluster_map = findNode::getClass<CMFlashClusterContainer>(topNode, "CORRECTED_CM_CLUSTER");
  if (!m_corrected_CMcluster_map)
  {
    std::cout << PHWHERE << "CORRECTED_CM_CLUSTER Node missing, abort." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // input tpc distortion correction static
  m_dcc_in_static = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerStatic");
  if (m_dcc_in_static)
  {
    std::cout << "PHTpcCentralMembraneMatcher:   found TPC distortion correction container static" << std::endl;
  }

  // input tpc distortion correction average
  m_dcc_in_average = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerAverage");
  if (m_dcc_in_average)
  {
    std::cout << "PHTpcCentralMembraneMatcher:   found TPC distortion correction container average" << std::endl;
  }

  // create node for results of matching
  std::cout << "Creating node CM_FLASH_DIFFERENCES" << std::endl;
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode* DetNode =
      dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode", "TRKR"));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode("TRKR");
    dstNode->addNode(DetNode);
  }

  m_cm_flash_diffs = new CMFlashDifferenceContainerv1;
  PHIODataNode<PHObject>* CMFlashDifferenceNode =
      new PHIODataNode<PHObject>(m_cm_flash_diffs, "CM_FLASH_DIFFERENCES", "PHObject");
  DetNode->addNode(CMFlashDifferenceNode);

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
      std::cout << "PHTpcCentralMembraneMatcher::InitRun - RUN Node missing, quitting" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    std::cout << "PHTpcCentralMembraneMatcher::GetNodes - creating TpcDistortionCorrectionContainer in node " << dcc_out_node_name << std::endl;
    m_dcc_out = new TpcDistortionCorrectionContainer;
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
      dcc->m_hDPint[i] = new TH2F(
          (boost::format("hIntDistortionP%s") % extension[i]).str().c_str(),
          (boost::format("hIntDistortionP%s") % extension[i]).str().c_str(),
          m_phibins + 2, phiMin, phiMax, m_rbins + 2, rMin, rMax);
      delete dcc->m_hDRint[i];
      dcc->m_hDRint[i] = new TH2F(
          (boost::format("hIntDistortionR%s") % extension[i]).str().c_str(),
          (boost::format("hIntDistortionR%s") % extension[i]).str().c_str(),
          m_phibins + 2, phiMin, phiMax, m_rbins + 2, rMin, rMax);
      delete dcc->m_hDZint[i];
      dcc->m_hDZint[i] = new TH2F(
          (boost::format("hIntDistortionZ%s") % extension[i]).str().c_str(),
          (boost::format("hIntDistortionZ%s") % extension[i]).str().c_str(),
          m_phibins + 2, phiMin, phiMax, m_rbins + 2, rMin, rMax);
      delete dcc->m_hentries[i];
      dcc->m_hentries[i] = new TH2I(
          (boost::format("hEntries%s") % extension[i]).str().c_str(),
          (boost::format("hEntries%s") % extension[i]).str().c_str(),
          m_phibins + 2, phiMin, phiMax, m_rbins + 2, rMin, rMax);

      // delete dcc->m_hDPint[i]; dcc->m_hDPint[i] = new TH2F( Form("hIntDistortionP%s", extension[i].c_str()), Form("hIntDistortionP%s", extension[i].c_str()), 205, phiBins, 68, r_bins );
      // delete dcc->m_hDRint[i]; dcc->m_hDRint[i] = new TH2F( Form("hIntDistortionR%s", extension[i].c_str()), Form("hIntDistortionR%s", extension[i].c_str()), 205, phiBins, 68, r_bins );
      // delete dcc->m_hDZint[i]; dcc->m_hDZint[i] = new TH2F( Form("hIntDistortionZ%s", extension[i].c_str()), Form("hIntDistortionZ%s", extension[i].c_str()), 205, phiBins, 68, r_bins );
      // delete dcc->m_hentries[i]; dcc->m_hentries[i] = new TH2I( Form("hEntries%s", extension[i].c_str()), Form("hEntries%s", extension[i].c_str()), 205, phiBins, 68, r_bins);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________
void PHTpcCentralMembraneMatcher::CalculateCenters(
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
