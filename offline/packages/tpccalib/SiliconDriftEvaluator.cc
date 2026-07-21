#include "SiliconDriftEvaluator.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TParameter.h>

#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <memory>

//_____________________________________________________________________
namespace
{

  //! pt
  template <class T>
  T get_pt(const T& px, const T& py)
  {
    return std::sqrt(px * px + py * py);
  }

  //_____________________________________________________________________
  //  par[0] = constrained slope
  //  par[1] = offset for eta < 0
  //  par[2] = offset for eta >= 0
  //
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

  //! 1D version used to draw per-eta overlay lines on QA canvas
  double linear_function(double* x, double* par)
  {
    return par[0] * x[0] + par[1];
  }

  //! human-readable label for each eta bin
  const char* k_eta_labels[2] = {"#eta_{TPC} < 0", "#eta_{TPC} #geq 0"};

}  // namespace

//_____________________________________________________________________
SiliconDriftEvaluator::SiliconDriftEvaluator(const std::string& name)
  : SubsysReco(name)
{
}

//_____________________________________________________________________
int SiliconDriftEvaluator::Init(PHCompositeNode* topNode)
{
  // find DST node
  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "SiliconDriftEvaluator::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto* evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if (!evalNode)
  {
    // create
    std::cout << "SiliconDriftEvaluator::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode("EVAL");
    dstNode->addNode(evalNode);
  }

  // add container to output tree
  auto* newNode = new PHIODataNode<PHObject>(new Container, "SiliconDriftEvaluator::Container", "PHObject");

  // overwrite split level for easier offline browsing
  newNode->SplitLevel(99);
  evalNode->addNode(newNode);

  // book 3D accumulator histogram
  // x = eta bin: 0 = eta<0, 1 = eta>=0
  // y = z_si (cm)
  // z = dz   (cm)
  m_hist3D = new TH3F("SiliconDriftEval_hist3D", ";#eta bin;z_{silicon} (cm);#Deltaz_{TPC-silicon} (cm)", 2, 0, 2, 200, -m_max_z, m_max_z, 200, -m_max_dz, m_max_dz);
  m_hist3D->SetDirectory(nullptr);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SiliconDriftEvaluator::InitRun(PHCompositeNode* topNode)
{
  return load_nodes(topNode);
}

//_____________________________________________________________________
int SiliconDriftEvaluator::process_event(PHCompositeNode* topNode)
{
  // load nodes
  const auto res = load_nodes(topNode);
  if (res != Fun4AllReturnCodes::EVENT_OK) {return res;}

  // cleanup output
  if (m_container) {m_container->Reset();}

  evaluate_tracks();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SiliconDriftEvaluator::End(PHCompositeNode*)
{
  if (!m_hist3D)
  {
    std::cerr << Name() << "::End - histogram not found, skipping fit." << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  const int nEntries = static_cast<int>(m_hist3D->GetEntries());
  std::cout << Name() << "::End - fitting " << nEntries << " entries" << std::endl;

  // build mean-dz TH2F via FitSlicesY, one eta bin at a time
  // x = eta bin [0,2), y = z_si (cm), content = mean dz (cm)
  auto* h_fit = new TH2F("h_fit_silicon", "",
                        2, 0, 2,
                        200, -m_max_z, m_max_z);
  h_fit->SetDirectory(nullptr);

  for (int ieta = 0; ieta < 2; ++ieta)
  {
    m_hist3D->GetXaxis()->SetRange(ieta + 1, ieta + 1);
    auto* h2d = static_cast<TH2F*>(m_hist3D->Project3D("zy"));
    h2d->SetName(Form("h2d_etabin_%i", ieta));
    h2d->SetDirectory(nullptr);

    // fit vertical slices; require a minimum of m_min_slice_entries per slice
    h2d->FitSlicesY(nullptr, 0, -1, m_min_slice_entries);
    auto* h_mean = static_cast<TH1F*>(gDirectory->Get(Form("h2d_etabin_%i_1", ieta)));

    if (!h_mean)
    {
      delete h2d;
      continue;
    }

    for (int iz = 1; iz <= h_mean->GetNbinsX(); ++iz)
    {
      const double entries = h2d->Integral(iz, iz, 1, m_hist3D->GetNbinsZ());
      if (entries > 0)
      {
        h_fit->SetBinContent(ieta + 1, iz, h_mean->GetBinContent(iz));
      }
    }

    delete h2d;
  }

  m_hist3D->GetXaxis()->SetRange(0, 0);

  // 2D piecewise fit: shared slope + per-eta offset
  auto* fit2d = new TF2("fit2d_silicon", fit_function_2d, 0, 2, -m_max_z, m_max_z, 3);
  for (int i = 0; i < 3; ++i) {fit2d->SetParameter(i, 0.0);}
  h_fit->Fit(fit2d, "0R");

  const double slope = fit2d->GetParameter(0);
  const double slope_err = fit2d->GetParError(0);
  const double off_neg = fit2d->GetParameter(1);  // ieta=0, eta<0
  const double off_pos = fit2d->GetParameter(2);  // ieta=1, eta>=0

  const double dv_new = m_drift_velocity / (1.0 + slope);
  const double dv_err = m_drift_velocity / std::pow(1.0 + slope, 2) * slope_err;
  const double t0_new = (off_pos - off_neg) / (2.0 * dv_new);

  std::cout << Name() << "::End" << " slope=" << slope << " dv_in=" << m_drift_velocity << " cm/ns" << " dv_new=" << dv_new << " +/- " << dv_err << " cm/ns" << " t0_new=" << t0_new << " ns" << std::endl;

  // draw the plot
  auto* canvas = new TCanvas("silicon_drift_calib", "Silicon drift velocity calibration", 1400, 700);
  canvas->Divide(2, 1);

  for (int ieta = 0; ieta < 2; ++ieta)
  {
    canvas->cd(ieta + 1);
    gPad->SetTopMargin(0.13);
    gPad->SetRightMargin(0.18);

    // 2D distribution for this eta bin
    m_hist3D->GetXaxis()->SetRange(ieta + 1, ieta + 1);
    auto* h2d = static_cast<TH2F*>(m_hist3D->Project3D("zy"));
    h2d->SetName(Form("hplot_etabin_%i", ieta));
    h2d->SetTitle(";z_{silicon} (cm);#Deltaz_{TPC-silicon} (cm)");
    h2d->SetStats(0);
    h2d->Draw("COLZ");

    // mean-dz points from FitSlicesY
    auto* h_fit_proj = h_fit->ProjectionY(Form("h_fit_proj_%i", ieta), ieta + 1, ieta + 1);
    h_fit_proj->SetMarkerStyle(20);
    h_fit_proj->SetMarkerSize(0.6);
    h_fit_proj->SetMarkerColor(kRed);
    h_fit_proj->SetLineColor(kRed);
    h_fit_proj->Draw("same P");

    // 1D fit line for this eta bin
    auto* f1d = new TF1(Form("f1d_etabin_%i", ieta), linear_function, -m_max_z, m_max_z, 2);
    f1d->SetParameter(0, slope);
    f1d->SetParameter(1, (ieta == 0) ? off_neg : off_pos);
    f1d->SetLineColor(kGreen + 2);
    f1d->SetLineWidth(2);
    f1d->Draw("same");

    // reference line at dz = 0
    auto* zero = new TLine(-m_max_z, 0, m_max_z, 0);
    zero->SetLineStyle(2);
    zero->SetLineColor(kGray + 1);
    zero->Draw();

    auto* leg = new TLegend(0.13, 0.76, 0.82, 0.95);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.033);
    leg->SetHeader(Form("%s   entries: %i   v_{in}=%.4f cm/ns",
                        k_eta_labels[ieta], nEntries, m_drift_velocity),
                   "C");
    leg->AddEntry(h_fit_proj, "Gaussian slice mean", "p");
    leg->AddEntry(f1d, Form("slope=%.4f  v_{new}=%.4f#pm%.4f cm/ns  t_{0}=%.1f ns", slope, dv_new, dv_err, t0_new), "l");
    leg->Draw();
  }

  m_hist3D->GetXaxis()->SetRange(0, 0);

  canvas->SaveAs(m_plot_filename.c_str());
  std::cout << Name() << "::End - QA canvas saved to " << m_plot_filename << std::endl;

  // write histograms, fit and results to a ROOT file
  if (!m_root_filename.empty())
  {
    std::unique_ptr<TFile> outfile(TFile::Open(m_root_filename.c_str(), "RECREATE"));
    if (outfile && !outfile->IsZombie())
    {
      outfile->cd();
      m_hist3D->Write();
      h_fit->Write("h_fit_silicon");
      fit2d->Write();
      canvas->Write();
      TParameter<double>("slope", slope).Write();
      TParameter<double>("drift_velocity_in", m_drift_velocity).Write();
      TParameter<double>("drift_velocity_new", dv_new).Write();
      TParameter<double>("drift_velocity_err", dv_err).Write();
      TParameter<double>("t0_new", t0_new).Write();
      outfile->Close();
      std::cout << Name() << "::End - histograms and fit results saved to " << m_root_filename << std::endl;
    }
    else
    {
      std::cerr << Name() << "::End - could not open " << m_root_filename << " for writing." << std::endl;
    }
  }

  delete canvas;
  delete fit2d;
  delete h_fit;

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SiliconDriftEvaluator::load_nodes(PHCompositeNode* topNode)
{
  // track map
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, m_trackmapname);

  // local container
  m_container = findNode::getClass<Container>(topNode, "SiliconDriftEvaluator::Container");
  assert(m_container);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void SiliconDriftEvaluator::evaluate_tracks()
{
  if (!(m_track_map && m_container && m_hist3D)) {return;}

  // clear array
  m_container->clearTracks();

  for (const auto& [track_id, track] : *m_track_map)
  {
    // require valid beam-crossing
    const auto crossing = track->get_crossing();
    if (crossing == SHRT_MAX)
    {
      std::cout << "SiliconDriftEvaluator::evaluate_tracks - invalid crossing, track ignored." << std::endl;
      continue;
    }

    // require both seeds
    const auto* si_seed = track->get_silicon_seed();
    const auto* tpc_seed = track->get_tpc_seed();
    if (!si_seed || !tpc_seed) continue;

    // count clusters per subsystem
    unsigned int n_tpc = 0;
    unsigned int n_mvtx = 0;
    unsigned int n_intt = 0;

    for (const auto* seed : {track->get_silicon_seed(), track->get_tpc_seed()})
    {
      if (!seed) continue;
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
    if (n_tpc < m_min_nclusters_tpc) continue;
    if (n_mvtx < m_min_nclusters_mvtx) continue;
    if (n_intt < m_min_nclusters_intt) continue;

    const float eta = tpc_seed->get_eta();
    if (std::abs(eta) > m_max_eta) continue;

    const float pt = get_pt(track->get_px(), track->get_py());
    if (pt < m_min_pt) continue;

    // get seed z positions at POCA
    const auto si_pos = TrackSeedHelper::get_xyz(si_seed);
    const auto tpc_pos = TrackSeedHelper::get_xyz(tpc_seed);

    const float z_si = si_pos.z();
    const float z_tpc = tpc_pos.z();

    // dz = (z_tpc + sign(eta)*crossing*crossing_interval*dv) - z_si
    const double sign_eta = (eta >= 0) ? 1.0 : -1.0;
    const float z_tpc_corr = z_tpc + sign_eta * crossing * m_crossing_interval * m_drift_velocity;
    const float dz = z_tpc_corr - z_si;

    // fill track struct
    TrackStruct track_struct;
    track_struct._nclusters_mvtx = n_mvtx;
    track_struct._nclusters_intt = n_intt;
    track_struct._nclusters_tpc = n_tpc;
    track_struct._pt = pt;
    track_struct._eta = eta;
    track_struct._phi = tpc_seed->get_phi();
    track_struct._z_tpc = z_tpc;
    track_struct._z_si = z_si;
    track_struct._crossing = crossing;
    track_struct._dz = dz;

    // fill histogram
    // eta bin centre: 0.5 for eta<0, 1.5 for eta>=0
    const double eta_bin = (eta >= 0) ? 1.5 : 0.5;
    m_hist3D->Fill(eta_bin, z_si, dz);

    m_container->addTrack(track_struct);
  }
}