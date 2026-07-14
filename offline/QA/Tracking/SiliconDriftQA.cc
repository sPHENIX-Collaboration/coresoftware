#include "SiliconDriftQA.h"

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <TDirectory.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>

#include <cassert>
#include <climits>
#include <cmath>
#include <format>
#include <iostream>
#include <string>

namespace
{
  //! pt
  template <class T>
  T get_pt(const T& px, const T& py)
  {
    return std::sqrt(px * px + py * py);
  }

  //! piecewise fit function used for the drift velocity extraction
  //  par[0] = constrained slope
  //  par[1] = offset for eta < 0
  //  par[2] = offset for eta >= 0
  double fit_function_2d(double* x, double* par)
  {
    const int ieta = static_cast<int>(std::floor(x[0]));
    const double z = x[1];
    if (ieta < 0 || ieta > 1)
    {
      TF2::RejectPoint();
      return 0.;
    }
    return par[ieta + 1] + par[0] * z;
  }

  //! suffixes used in histogram names for the two eta bins
  const char* k_eta_suffix[2] = {"negeta", "poseta"};

  //! number of z bins of the dz vs z histograms
  constexpr int k_nzbins = 200;

}  // namespace

//____________________________________________________________________________..
SiliconDriftQA::SiliconDriftQA(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int SiliconDriftQA::InitRun(PHCompositeNode* /*topNode*/)
{
  createHistos();

  // reference histograms initialized in header file to histos in HistoManager
  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  for (int ieta = 0; ieta < 2; ieta++)
  {
    h_zsi_dz[ieta] = dynamic_cast<TH2*>(hm->getHisto(std::format("{}zsi_dz_{}", getHistoPrefix(), k_eta_suffix[ieta])));
  }
  h_dz = dynamic_cast<TH1*>(hm->getHisto(std::format("{}dz", getHistoPrefix())));
  h_ntracks = dynamic_cast<TH1*>(hm->getHisto(std::format("{}ntracks", getHistoPrefix())));
  h_driftSummary = dynamic_cast<TH1*>(hm->getHisto(std::format("{}driftSummary", getHistoPrefix())));

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SiliconDriftQA::process_event(PHCompositeNode* topNode)
{
  auto* track_map = findNode::getClass<SvtxTrackMap>(topNode, m_trackmapname);
  if (!track_map)
  {
    std::cout << PHWHERE << " " << m_trackmapname << " node missing, abort." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  int naccepted = 0;

  for (const auto& [track_id, track] : *track_map)
  {
    // require valid beam-crossing
    const auto crossing = track->get_crossing();
    if (crossing == SHRT_MAX)
    {
      if (Verbosity())
      {
        std::cout << PHWHERE << " invalid crossing, track ignored." << std::endl;
      }
      continue;
    }

    // require both seeds
    const auto* si_seed = track->get_silicon_seed();
    const auto* tpc_seed = track->get_tpc_seed();
    if (!si_seed || !tpc_seed)
    {
      continue;
    }

    // count clusters per subsystem
    unsigned int n_tpc = 0;
    unsigned int n_mvtx = 0;
    unsigned int n_intt = 0;

    for (const auto* seed : {si_seed, tpc_seed})
    {
      for (auto it = seed->begin_cluster_keys(); it != seed->end_cluster_keys(); ++it)
      {
        switch (TrkrDefs::getTrkrId(*it))
        {
        case TrkrDefs::tpcId:
          ++n_tpc;
          break;
        case TrkrDefs::mvtxId:
          ++n_mvtx;
          break;
        case TrkrDefs::inttId:
          ++n_intt;
          break;
        default:
          break;
        }
      }
    }

    // apply selection cuts
    if (n_tpc < m_min_nclusters_tpc)
    {
      continue;
    }
    if (n_mvtx < m_min_nclusters_mvtx)
    {
      continue;
    }
    if (n_intt < m_min_nclusters_intt)
    {
      continue;
    }

    const float eta = tpc_seed->get_eta();
    if (std::abs(eta) > m_max_eta)
    {
      continue;
    }

    const float pt = get_pt(track->get_px(), track->get_py());
    if (pt < m_min_pt)
    {
      continue;
    }

    // get seed z positions at POCA
    const auto si_pos = TrackSeedHelper::get_xyz(si_seed);
    const auto tpc_pos = TrackSeedHelper::get_xyz(tpc_seed);

    const float z_si = si_pos.z();
    const float z_tpc = tpc_pos.z();

    // dz = (z_tpc + sign(eta)*crossing*crossing_interval*dv) - z_si
    const double sign_eta = (eta >= 0) ? 1.0 : -1.0;
    const float z_tpc_corr = z_tpc + sign_eta * crossing * m_crossing_interval * m_drift_velocity;
    const float dz = z_tpc_corr - z_si;

    // fill histograms
    const int ieta = (eta >= 0) ? 1 : 0;
    h_zsi_dz[ieta]->Fill(z_si, dz);
    h_dz->Fill(dz);

    ++naccepted;
  }

  h_ntracks->Fill(naccepted);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SiliconDriftQA::End(PHCompositeNode* /*topNode*/)
{
  if (!(h_zsi_dz[0] && h_zsi_dz[1] && h_driftSummary))
  {
    std::cout << PHWHERE << " histograms not found, skipping drift velocity fit." << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  const int nEntries = static_cast<int>(h_zsi_dz[0]->GetEntries() + h_zsi_dz[1]->GetEntries());
  if (Verbosity())
  {
    std::cout << Name() << "::End - fitting " << nEntries << " entries" << std::endl;
  }

  // record input drift velocity even if the fit is skipped
  h_driftSummary->SetBinContent(3, m_drift_velocity);

  if (nEntries < 2 * m_min_slice_entries)
  {
    std::cout << Name() << "::End - not enough entries (" << nEntries << "), skipping drift velocity fit." << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // build mean-dz TH2 via FitSlicesY, one eta bin at a time
  // x = eta bin [0,2), y = z_si (cm), content = mean dz (cm)
  auto* h_fit = new TH2F("h_fit_silicon", "", 2, 0, 2, k_nzbins, -m_max_z, m_max_z);
  h_fit->SetDirectory(nullptr);

  for (int ieta = 0; ieta < 2; ++ieta)
  {
    auto* h2d = h_zsi_dz[ieta];

    // fit vertical slices; require a minimum of m_min_slice_entries per slice
    TObjArray slices;
    slices.SetOwner(kTRUE);
    h2d->FitSlicesY(nullptr, 0, -1, m_min_slice_entries, "QNR", &slices);
    auto* h_mean = dynamic_cast<TH1*>(slices.At(1));
    if (!h_mean)
    {
      continue;
    }

    for (int iz = 1; iz <= h_mean->GetNbinsX(); ++iz)
    {
      const double entries = h2d->Integral(iz, iz, 1, h2d->GetNbinsY());
      if (entries > 0)
      {
        h_fit->SetBinContent(ieta + 1, iz, h_mean->GetBinContent(iz));
      }
    }
  }

  // 2D piecewise fit: shared slope + per-eta offset
  auto* fit2d = new TF2("fit2d_silicon", fit_function_2d, 0, 2, -m_max_z, m_max_z, 3);
  for (int i = 0; i < 3; ++i)
  {
    fit2d->SetParameter(i, 0.0);
  }
  h_fit->Fit(fit2d, "0RQ");

  const double slope = fit2d->GetParameter(0);
  const double slope_err = fit2d->GetParError(0);
  const double off_neg = fit2d->GetParameter(1);  // ieta=0, eta<0
  const double off_pos = fit2d->GetParameter(2);  // ieta=1, eta>=0

  const double dv_new = m_drift_velocity / (1.0 + slope);
  const double dv_err = m_drift_velocity / std::pow(1.0 + slope, 2) * slope_err;
  const double t0_new = (off_pos - off_neg) / (2.0 * dv_new);

  std::cout << Name() << "::End"
            << " slope=" << slope
            << " dv_in=" << m_drift_velocity << " cm/ns"
            << " dv_new=" << dv_new << " +/- " << dv_err << " cm/ns"
            << " t0_new=" << t0_new << " ns"
            << std::endl;

  // store fit results in the summary histogram
  h_driftSummary->SetBinContent(1, slope);
  h_driftSummary->SetBinContent(2, slope_err);
  h_driftSummary->SetBinContent(4, dv_new);
  h_driftSummary->SetBinContent(5, dv_err);
  h_driftSummary->SetBinContent(6, t0_new);

  delete fit2d;
  delete h_fit;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
std::string SiliconDriftQA::getHistoPrefix() const
{
  // define prefix to all histos in HistoManager
  return std::string("h_") + Name() + std::string("_");
}

//____________________________________________________________________________..
void SiliconDriftQA::createHistos()
{
  // initialize HistoManager
  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // create and register histos in HistoManager
  for (int ieta = 0; ieta < 2; ieta++)
  {
    auto* h = new TH2F(std::format("{}zsi_dz_{}", getHistoPrefix(), k_eta_suffix[ieta]).c_str(),
                       std::format("{};z_{{silicon}} (cm);#Deltaz_{{TPC-silicon}} (cm)",
                                   (ieta == 0 ? "#eta_{TPC} < 0" : "#eta_{TPC} #geq 0"))
                           .c_str(),
                       k_nzbins, -m_max_z, m_max_z, 200, -m_max_dz, m_max_dz);
    hm->registerHisto(h);
  }

  {
    auto* h = new TH1F(std::format("{}dz", getHistoPrefix()).c_str(),
                       ";#Deltaz_{TPC-silicon} (cm);tracks", 200, -m_max_dz, m_max_dz);
    hm->registerHisto(h);
  }

  {
    auto* h = new TH1F(std::format("{}ntracks", getHistoPrefix()).c_str(),
                       ";accepted tracks per event;events", 50, -0.5, 49.5);
    hm->registerHisto(h);
  }

  {
    // summary of the drift velocity fit performed in End()
    auto* h = new TH1F(std::format("{}driftSummary", getHistoPrefix()).c_str(),
                       "drift velocity fit summary", 6, 0.5, 6.5);
    h->GetXaxis()->SetBinLabel(1, "slope");
    h->GetXaxis()->SetBinLabel(2, "slope_err");
    h->GetXaxis()->SetBinLabel(3, "v_{in} (cm/ns)");
    h->GetXaxis()->SetBinLabel(4, "v_{new} (cm/ns)");
    h->GetXaxis()->SetBinLabel(5, "v_{new} err (cm/ns)");
    h->GetXaxis()->SetBinLabel(6, "t_{0} (ns)");
    hm->registerHisto(h);
  }
}