#include "INTTZvtx.h"

#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLine.h>
#include <TTree.h>

#include <boost/format.hpp>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

INTTZvtx::INTTZvtx(const std::string& runType,
                   const std::string& outFolderDirectory,
                   std::pair<double, double> beamOrigin,
                   double phiDiffCut,
                   std::pair<double, double> DCACut,
                   int NCluCutl,
                   int NCluCut,
                   unsigned int zvtxCalRequire,
                   std::pair<double, double> zvtxQAWidth,
                   bool drawEventDisplay,
                   bool enableQA,
                   bool printMessageOpt)
  : run_type(runType)
  , out_folder_directory(outFolderDirectory)
  , beam_origin(beamOrigin)
  , phi_diff_cut(phiDiffCut)
  , DCA_cut(DCACut)
  , N_clu_cut(NCluCut)
  , N_clu_cutl(NCluCutl)
  , zvtx_cal_require(zvtxCalRequire)
  , zvtx_QA_width(zvtxQAWidth)
  , draw_event_display(drawEventDisplay)
  , m_enable_qa(enableQA)
  , print_message_opt(printMessageOpt)
{
  // SetsPhenixStyle();
  gErrorIgnoreLevel = kWarning;  // note : To not print the "print plot info."

  temp_event_zvtx_info.clear();
  avg_event_zvtx_vec.clear();
  Z_resolution_vec.clear();
  good_comb_id = 0;
  MC_z_diff_peak = -777.;
  MC_z_diff_width = -777.;

  N_group_info_detail = {-1., -1., -1., -1.};

  N_comb.clear();
  N_comb_e.clear();
  z_mid.clear();
  z_range.clear();
  N_comb_phi.clear();
  // eff_N_comb.clear(); eff_z_mid.clear(); eff_N_comb_e.clear(); eff_z_range.clear(); // note : eff_sig

  inner_clu_phi_map.clear();
  outer_clu_phi_map.clear();
  inner_clu_phi_map = std::vector<std::vector<std::pair<bool, clu_info>>>(360);
  outer_clu_phi_map = std::vector<std::vector<std::pair<bool, clu_info>>>(360);

  ///////////
  // ntuple variables
  out_eID = 0;
  bco_full_out = 0;
  N_cluster_inner_out = -1;
  N_cluster_outer_out = -1;
  out_ES_zvtx = -1;
  out_ES_zvtxE = -1;
  out_ES_rangeL = -1;
  out_ES_rangeR = -1;
  out_ES_N_good = -1;
  out_ES_width_density = -1;

  out_LB_Gaus_Mean_mean = -1;
  out_LB_Gaus_Mean_meanE = -1;
  out_LB_Gaus_Mean_chi2 = -1;
  out_LB_Gaus_Mean_width = -1;

  out_LB_Gaus_Width_width = -1;
  out_LB_Gaus_Width_size_width = -1;
  out_LB_Gaus_Width_offset = -1;

  out_mid_cut_Ngroup = -1;
  out_mid_cut_peak_width = -1;
  out_mid_cut_peak_ratio = -1;

  out_LB_cut_Ngroup = -1;
  out_LB_cut_peak_width = -1;
  out_LB_cut_peak_ratio = -1;

  out_LB_geo_mean = -1;
  out_good_zvtx_tag = false;

  MC_true_zvtx = -9999.;

  out_centrality_bin = -1;

  out_N_cluster_north = 0;
  out_N_cluster_south = 0;

  ///////////
  //    Init();
}

INTTZvtx::~INTTZvtx()
{
  // histos for z-vertex calculation
  if (evt_possible_z != nullptr)
  {
    delete evt_possible_z;
  }
  if (line_breakdown_hist != nullptr)
  {
    delete line_breakdown_hist;
  }
  if (gaus_fit != nullptr)
  {
    delete gaus_fit;
  }
  if (zvtx_finder != nullptr)
  {
    delete zvtx_finder;
  }

  if (z_range_gr != nullptr)
  {
    delete z_range_gr;
  }

  if (draw_event_display)
  {
    // QA histograms
    // event by event
    delete evt_select_track_phi;
    delete evt_phi_diff_1D;
    delete evt_phi_diff_inner_phi;
    delete evt_inner_outer_phi;
    delete phi_diff_inner_phi;

    delete c2;
    // all the pads related to c2 are automatically deleted
    if (temp_event_xy != nullptr)
    {
      delete temp_event_xy;
    }
    if (temp_event_rz != nullptr)
    {
      delete temp_event_rz;
    }
    if (z_range_gr_draw != nullptr)
    {
      delete z_range_gr_draw;
    }
  }

  if (m_enable_qa)
  {
    delete c1;

    std::cout << "del : " << std::hex << (long) avg_event_zvtx << std::dec << std::endl;
    for (auto& itr : m_v_qahist)
    {
      // std::cout<<"del : "<<itr->GetTitle()<<std::endl;
      delete itr;
    }
  }

  if (draw_event_display)
  {
    delete final_fit_range_line;
    delete coord_line;
    delete ladder_line;
    delete bkg;
  }

  if (draw_event_display || m_enable_qa)
  {
    delete eff_sig_range_line;
    delete draw_text;
  }
}

void INTTZvtx::Characterize_Pad(TPad* pad, float left, float right, float top, float bottom, bool set_logY, int setgrid_bool)
{
  if (setgrid_bool == true)
  {
    pad->SetGrid(1, 1);
  }
  pad->SetLeftMargin(left);
  pad->SetRightMargin(right);
  pad->SetTopMargin(top);
  pad->SetBottomMargin(bottom);
  pad->SetTicks(1, 1);
  if (set_logY == true)
  {
    pad->SetLogy(1);
  }
}

void INTTZvtx::Init()
{
  if (!std::filesystem::exists(out_folder_directory))
  {
    std::filesystem::create_directory(out_folder_directory);
  }

  InitCanvas();
  InitHist();
  InitTreeOut();
  InitRest();

  if (draw_event_display)
  {
    c2->Print((boost::format("%s/temp_event_display.pdf") % out_folder_directory).str().c_str());
  }

  m_initialized = true;
}

void INTTZvtx::InitHist()
{
  // histos for z-vertex calculation
  evt_possible_z = new TH1F("evt_possible_z", "evt_possible_z", 50, evt_possible_z_range.first, evt_possible_z_range.second);
  evt_possible_z->SetLineWidth(1);
  evt_possible_z->GetXaxis()->SetTitle("Z [mm]");
  evt_possible_z->GetYaxis()->SetTitle("Entry");

  int N = 1200;        // note : N bins for each side, regardless the bin at zero
  double width = 0.5;  // note : bin width with the unit [mm]
  line_breakdown_hist = new TH1F("line_breakdown_hist", "line_breakdown_hist", 2 * N + 1, -1 * (width * N + width / 2.), width * N + width / 2.);
  line_breakdown_hist->SetLineWidth(1);
  line_breakdown_hist->GetXaxis()->SetTitle("Z [mm]");
  line_breakdown_hist->GetYaxis()->SetTitle("Entry");
  if (print_message_opt == true)
  {
    std::cout << "class INTTZvtx, Line breakdown hist, range : "
              << line_breakdown_hist->GetXaxis()->GetXmin() << " "
              << line_breakdown_hist->GetXaxis()->GetXmax() << " "
              << line_breakdown_hist->GetBinWidth(1) << std::endl;
  }

  // QA histograms
  // event by event
  if (draw_event_display)
  {
    evt_select_track_phi = new TH1F("evt_select_track_phi", "evt_select_track_phi", 361, 0, 361);
    evt_select_track_phi->GetXaxis()->SetTitle("Track phi [degree]");
    evt_select_track_phi->GetYaxis()->SetTitle("Entry");
    evt_select_track_phi->GetXaxis()->SetNdivisions(505);

    evt_phi_diff_inner_phi = new TH2F("evt_phi_diff_inner_phi", "evt_phi_diff_inner_phi", 361, 0, 361, 100, -1.5, 1.5);
    evt_phi_diff_inner_phi->GetXaxis()->SetTitle("Inner phi [degree]");
    evt_phi_diff_inner_phi->GetYaxis()->SetTitle("Inner - Outer [degree]");
    evt_phi_diff_inner_phi->GetXaxis()->SetNdivisions(505);

    evt_inner_outer_phi = new TH2F("evt_inner_outer_phi", "evt_inner_outer_phi", 120, 0, 360, 120, 0, 360);
    evt_inner_outer_phi->GetXaxis()->SetTitle("Inner phi [degree]");
    evt_inner_outer_phi->GetYaxis()->SetTitle("Outer phi [degree]");
    evt_inner_outer_phi->GetXaxis()->SetNdivisions(505);

    evt_phi_diff_1D = new TH1F("evt_phi_diff_1D", "evt_phi_diff_1D", 100, -10, 10);
    evt_phi_diff_1D->GetXaxis()->SetTitle("Inner - Outer [degree]");
    evt_phi_diff_1D->GetYaxis()->SetTitle("Entry");
    evt_phi_diff_1D->GetXaxis()->SetNdivisions(505);

    phi_diff_inner_phi = new TH2F("phi_diff_inner_phi", "All evt phi_diff_inner_phi", 361, 0, 361, 100, -1.5, 1.5);
    phi_diff_inner_phi->GetXaxis()->SetTitle("Inner phi [degree]");
    phi_diff_inner_phi->GetYaxis()->SetTitle("Inner - Outer [degree]");
    phi_diff_inner_phi->GetXaxis()->SetNdivisions(505);
  }

  if (m_enable_qa)
  {
    // whole run
    avg_event_zvtx = new TH1F("avg_event_zvtx", "avg_event_zvtx", 100, zvtx_hist_l, zvtx_hist_r);
    avg_event_zvtx->SetLineColor(TColor::GetColor("#1A3947"));
    avg_event_zvtx->SetLineWidth(2);
    avg_event_zvtx->GetYaxis()->SetTitle("Entry");
    avg_event_zvtx->GetXaxis()->SetTitle("Z vertex position [mm]");
    avg_event_zvtx->GetYaxis()->SetRangeUser(0, 50);
    avg_event_zvtx->SetTitleSize(0.06, "X");
    avg_event_zvtx->SetTitleSize(0.06, "Y");
    avg_event_zvtx->GetXaxis()->SetTitleOffset(0.82);
    avg_event_zvtx->GetYaxis()->SetTitleOffset(1.1);
    avg_event_zvtx->GetXaxis()->CenterTitle(true);
    avg_event_zvtx->GetYaxis()->CenterTitle(true);
    avg_event_zvtx->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(avg_event_zvtx);

    zvtx_evt_fitError = new TH1F("zvtx_evt_fitError", "zvtx_evt_fitError", 150, 0, 50);
    zvtx_evt_fitError->GetXaxis()->SetTitle(" mm ");
    zvtx_evt_fitError->GetYaxis()->SetTitle("entry");
    zvtx_evt_fitError->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(zvtx_evt_fitError);

    zvtx_evt_fitError_corre = new TH2F("zvtx_evt_fitError_corre", "zvtx_evt_fitError_corre", 200, 0, 10000, 200, 0, 20);
    zvtx_evt_fitError_corre->GetXaxis()->SetTitle(" # of clusters ");
    zvtx_evt_fitError_corre->GetYaxis()->SetTitle(" #pm mm ");
    zvtx_evt_fitError_corre->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(zvtx_evt_fitError_corre);

    zvtx_evt_width_corre = new TH2F("zvtx_evt_width_corre", "zvtx_evt_width_corre", 200, 0, 10000, 200, 0, 300);
    zvtx_evt_width_corre->GetXaxis()->SetTitle(" # of clusters ");
    zvtx_evt_width_corre->GetYaxis()->SetTitle(" mm ");
    zvtx_evt_width_corre->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(zvtx_evt_width_corre);

    zvtx_evt_nclu_corre = new TH2F("zvtx_evt_nclu_corre", "zvtx_evt_nclu_corre", 200, 0, 10000, 200, -1000, 1000);
    zvtx_evt_nclu_corre->GetXaxis()->SetTitle(" # of clusters ");
    zvtx_evt_nclu_corre->GetYaxis()->SetTitle(" zvtx [mm] ");
    zvtx_evt_nclu_corre->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(zvtx_evt_nclu_corre);

    width_density = new TH1F("width_density", "width_density", 200, 0, 0.006);  // note : N good hits / width
    width_density->GetXaxis()->SetTitle(" (N good track / width) / NClus ");
    width_density->GetYaxis()->SetTitle(" Entry ");
    width_density->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(width_density);

    ES_width = new TH1F("ES_width", "ES_width", 200, 0, 150);
    ES_width->GetXaxis()->SetTitle(" Width [mm] ");
    ES_width->GetYaxis()->SetTitle(" Entry ");
    ES_width->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(ES_width);

    ES_width_ratio = new TH1F("ES_width_ratio", "ES_width_ratio", 200, 0, 60);
    ES_width_ratio->GetXaxis()->SetTitle(" NClus / Width ");
    ES_width_ratio->GetYaxis()->SetTitle(" Entry ");
    ES_width_ratio->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(ES_width_ratio);

    Z_resolution_Nclu = new TH2F("Z_resolution_Nclu", "Z_resolution_Nclu", 200, 0, 10000, 200, -100, 100);
    Z_resolution_Nclu->GetXaxis()->SetTitle(" # of clusters ");
    Z_resolution_Nclu->GetYaxis()->SetTitle("#Delta Z (Reco - True) [mm]");
    Z_resolution_Nclu->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(Z_resolution_Nclu);

    Z_resolution_pos = new TH2F("Z_resolution_pos", "Z_resolution_pos", 200, -350, 350, 200, -50, 50);
    Z_resolution_pos->GetXaxis()->SetTitle("True Zvtx [mm]");
    Z_resolution_pos->GetYaxis()->SetTitle("#Delta Z (Reco - True) [mm]");
    Z_resolution_pos->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(Z_resolution_pos);

    Z_resolution_pos_cut = new TH2F("Z_resolution_pos_cut", "Z_resolution_pos_cut", 200, -350, 350, 200, -50, 50);
    Z_resolution_pos_cut->GetXaxis()->SetTitle("True Zvtx [mm]");
    Z_resolution_pos_cut->GetYaxis()->SetTitle("#Delta Z (Reco - True) [mm]");
    Z_resolution_pos_cut->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(Z_resolution_pos_cut);

    Z_resolution = new TH1F("Z_resolution", "Z_resolution", 200, -30, 30);
    Z_resolution->GetXaxis()->SetTitle("#Delta Z (Reco - True) [mm]");
    Z_resolution->GetYaxis()->SetTitle(" Entry ");
    Z_resolution->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(Z_resolution);

    line_breakdown_gaus_ratio_hist = new TH1F("line_breakdown_gaus_ratio_hist", "line_breakdown_gaus_ratio_hist", 200, 0, 0.0005);
    line_breakdown_gaus_ratio_hist->GetXaxis()->SetTitle("(Norm. size) / width");
    line_breakdown_gaus_ratio_hist->GetYaxis()->SetTitle(" Entry ");
    line_breakdown_gaus_ratio_hist->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(line_breakdown_gaus_ratio_hist);

    line_breakdown_gaus_width_hist = new TH1F("line_breakdown_gaus_width_hist", "line_breakdown_gaus_width_hist", 200, 0, 100);
    line_breakdown_gaus_width_hist->GetXaxis()->SetTitle("Width [mm]");
    line_breakdown_gaus_width_hist->GetYaxis()->SetTitle(" Entry ");
    line_breakdown_gaus_width_hist->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(line_breakdown_gaus_width_hist);

    gaus_width_Nclu = new TH2F("gaus_width_Nclu", "gaus_width_Nclu", 200, 0, 10000, 200, 0, 100);
    gaus_width_Nclu->GetXaxis()->SetTitle(" # of clusters ");
    gaus_width_Nclu->GetYaxis()->SetTitle("Gaus fit width [mm]");
    gaus_width_Nclu->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(gaus_width_Nclu);

    gaus_rchi2_Nclu = new TH2F("gaus_rchi2_Nclu", "gaus_rchi2_Nclu", 200, 0, 10000, 200, 0, 1);
    gaus_rchi2_Nclu->GetXaxis()->SetTitle(" # of clusters ");
    gaus_rchi2_Nclu->GetYaxis()->SetTitle("Gaus fit #chi2^{2}/NDF");
    gaus_rchi2_Nclu->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(gaus_rchi2_Nclu);

    final_fit_width = new TH2F("final_fit_width", "final_fit_width", 200, 0, 10000, 200, 0, 400);
    final_fit_width->GetXaxis()->SetTitle(" # of clusters ");
    final_fit_width->GetYaxis()->SetTitle("Fit width [mm]");
    final_fit_width->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(final_fit_width);

    N_track_candidate_Nclu = new TH2F("N_track_candidate_Nclu", "N_track_candidate_Nclu", 200, 0, 10000, 200, 0, 13000);
    N_track_candidate_Nclu->GetXaxis()->SetTitle(" # of clusters ");
    N_track_candidate_Nclu->GetYaxis()->SetTitle("N tracklet candidates");
    N_track_candidate_Nclu->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(N_track_candidate_Nclu);

    peak_group_width_hist = new TH1F("peak_group_width_hist", "peak_group_width_hist", 100, 0, 300);
    peak_group_width_hist->GetXaxis()->SetTitle("Width [mm]");
    peak_group_width_hist->GetYaxis()->SetTitle("Entry");
    peak_group_width_hist->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(peak_group_width_hist);

    peak_group_ratio_hist = new TH1F("peak_group_ratio_hist", "peak_group_ratio_hist", 110, 0, 1.1);
    peak_group_ratio_hist->GetXaxis()->SetTitle("Peak group ratio");
    peak_group_ratio_hist->GetYaxis()->SetTitle("Entry");
    peak_group_ratio_hist->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(peak_group_ratio_hist);

    N_group_hist = new TH1F("N_group_hist", "N_group_hist", 30, 0, 30);
    N_group_hist->GetXaxis()->SetTitle("N group post cut");
    N_group_hist->GetYaxis()->SetTitle("Entry");
    N_group_hist->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(N_group_hist);

    peak_group_detail_width_hist = new TH1F("peak_group_detail_width_hist", "peak_group_detail_width_hist", 200, 0, 300);
    peak_group_detail_width_hist->GetXaxis()->SetTitle("Width [mm]");
    peak_group_detail_width_hist->GetYaxis()->SetTitle("Entry");
    peak_group_detail_width_hist->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(peak_group_detail_width_hist);

    peak_group_detail_ratio_hist = new TH1F("peak_group_detail_ratio_hist", "peak_group_detail_ratio_hist", 110, 0, 1.1);
    peak_group_detail_ratio_hist->GetXaxis()->SetTitle("Peak group ratio");
    peak_group_detail_ratio_hist->GetYaxis()->SetTitle("Entry");
    peak_group_detail_ratio_hist->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(peak_group_detail_ratio_hist);

    N_group_detail_hist = new TH1F("N_group_detail_hist", "N_group_detail_hist", 30, 0, 30);
    N_group_detail_hist->GetXaxis()->SetTitle("N group post cut");
    N_group_detail_hist->GetYaxis()->SetTitle("Entry");
    N_group_detail_hist->GetXaxis()->SetNdivisions(505);
    m_v_qahist.push_back(N_group_detail_hist);

    dca_inner_phi = new TH2F("dca_inner_phi", "All dca_inner_phi", 90, 0, 360, 100, -10., 10);
    dca_inner_phi->GetXaxis()->SetTitle("Inner phi [degree]");
    dca_inner_phi->GetYaxis()->SetTitle("dca [mm]");
    // dca_inner_phi -> GetXaxis() -> SetNdivisions(505);
    m_v_qahist.push_back(dca_inner_phi);
  }
}

void INTTZvtx::InitCanvas()
{
  if (draw_event_display)
  {
    c2 = new TCanvas("", "", 4000, 1600);
    c2->cd();
    pad_xy = new TPad("pad_xy", "", 0.0, 0.5, 0.2, 1.0);
    Characterize_Pad(pad_xy, 0.15, 0.1, 0.1, 0.2, false, 0);
    pad_xy->Draw();

    pad_rz = new TPad("pad_rz", "", 0.2, 0.5, 0.40, 1.0);
    Characterize_Pad(pad_rz, 0.15, 0.1, 0.1, 0.2, false, 0);
    pad_rz->Draw();

    pad_z = new TPad("pad_z", "", 0.40, 0.5, 0.6, 1.0);
    Characterize_Pad(pad_z, 0.15, 0.1, 0.1, 0.2, false, 0);
    pad_z->Draw();

    pad_z_hist = new TPad("pad_z_hist", "", 0.6, 0.5, 0.8, 1.0);
    Characterize_Pad(pad_z_hist, 0.15, 0.1, 0.1, 0.2, false, 0);
    pad_z_hist->Draw();

    pad_z_line = new TPad("pad_z_line", "", 0.8, 0.5, 1, 1.0);
    Characterize_Pad(pad_z_line, 0.15, 0.1, 0.1, 0.2, false, 0);
    pad_z_line->Draw();

    pad_phi_diff = new TPad("pad_phi_diff", "", 0.0, 0.0, 0.2, 0.5);
    Characterize_Pad(pad_phi_diff, 0.15, 0.1, 0.1, 0.2, false, 0);
    pad_phi_diff->Draw();

    pad_track_phi = new TPad("pad_track_phi", "", 0.2, 0.0, 0.40, 0.5);
    Characterize_Pad(pad_track_phi, 0.15, 0.1, 0.1, 0.2, false, 0);
    pad_track_phi->Draw();

    pad_inner_outer_phi = new TPad("pad_inner_outer_phi", "", 0.4, 0.0, 0.60, 0.5);
    Characterize_Pad(pad_inner_outer_phi, 0.15, 0.1, 0.1, 0.2, false, 0);
    pad_inner_outer_phi->SetLogz(1);
    pad_inner_outer_phi->Draw();

    pad_phi_diff_1D = new TPad("pad_phi_diff_1D", "", 0.6, 0.0, 0.80, 0.5);
    Characterize_Pad(pad_phi_diff_1D, 0.15, 0.1, 0.1, 0.2, false, 0);
    // pad_phi_diff_1D -> SetLogz(1);
    pad_phi_diff_1D->Draw();
  }

  if (m_enable_qa)
  {
    c1 = new TCanvas("", "", 950, 800);
    c1->cd();
  }
}

void INTTZvtx::InitTreeOut()
{
  if (m_enable_qa)
  {
    out_file = new TFile((boost::format("%s/INTT_zvtx.root") % out_folder_directory).str().c_str(), "RECREATE");
    tree_out = new TTree("tree_Z", "INTT Z info.");

    tree_out->Branch("eID", &out_eID);
    tree_out->Branch("bco_full", &bco_full_out);
    tree_out->Branch("nclu_inner", &N_cluster_inner_out);
    tree_out->Branch("nclu_outer", &N_cluster_outer_out);
    tree_out->Branch("nclu_south", &out_N_cluster_south);
    tree_out->Branch("nclu_north", &out_N_cluster_north);
    tree_out->Branch("ES_zvtx", &out_ES_zvtx);                                    // note : effective sigma, pol0 fit z-vertex
    tree_out->Branch("ES_zvtxE", &out_ES_zvtxE);                                  // note : effective sigma, pol0 fit z-vertex error
    tree_out->Branch("ES_rangeL", &out_ES_rangeL);                                // note : effective sigma, selected range left
    tree_out->Branch("ES_rangeR", &out_ES_rangeR);                                // note : effective sigma, selected range right
    tree_out->Branch("ES_N_good", &out_ES_N_good);                                // note : effective sigma, number of z-vertex candidates in the range
    tree_out->Branch("ES_Width_density", &out_ES_width_density);                  // note : effective sigma, N good z-vertex candidates divided by width
    tree_out->Branch("LB_Gaus_Mean_mean", &out_LB_Gaus_Mean_mean);                // note : Line break loose offset - gaus mean
    tree_out->Branch("LB_Gaus_Mean_meanE", &out_LB_Gaus_Mean_meanE);              // note : Line break loose offset - gaus mean error
    tree_out->Branch("LB_Gaus_Mean_chi2", &out_LB_Gaus_Mean_chi2);                // note : Line break loose offset - reduce chi2
    tree_out->Branch("LB_Gaus_Mean_width", &out_LB_Gaus_Mean_width);              // note : Line break loose offset - width
    tree_out->Branch("LB_Gaus_Width_width", &out_LB_Gaus_Width_width);            // note : Line break tight offset - gaus width
    tree_out->Branch("LB_Gaus_Width_offset", &out_LB_Gaus_Width_offset);          // note : Line break tight offset - offset
    tree_out->Branch("LB_Gaus_Width_size_width", &out_LB_Gaus_Width_size_width);  // note : Line break tight offset - norm. height / width
    tree_out->Branch("LB_geo_mean", &out_LB_geo_mean);                            // note : Line break peak position directly from the distribution (with the bin width 0.5 mm)
    tree_out->Branch("good_zvtx_tag", &out_good_zvtx_tag);
    tree_out->Branch("mid_cut_Ngroup", &out_mid_cut_Ngroup);          // note : mid cut Ngroup
    tree_out->Branch("mid_cut_peak_width", &out_mid_cut_peak_width);  // note : mid cut peak width
    tree_out->Branch("mid_cut_peak_ratio", &out_mid_cut_peak_ratio);  // note : mid cut peak ratio
    tree_out->Branch("LB_cut_Ngroup", &out_LB_cut_Ngroup);            // note : LB cut Ngroup
    tree_out->Branch("LB_cut_peak_width", &out_LB_cut_peak_width);    // note : LB cut peak width
    tree_out->Branch("LB_cut_peak_ratio", &out_LB_cut_peak_ratio);    // note : LB cut peak ratio
    tree_out->Branch("MC_true_zvtx", &MC_true_zvtx);
    tree_out->Branch("Centrality_bin", &out_centrality_bin);
  }
}

void INTTZvtx::InitRest()
{
  gaus_fit = new TF1("gaus_fit", InttVertexUtil::gaus_func, evt_possible_z_range.first, evt_possible_z_range.second, 4);
  gaus_fit->SetLineColor(2);
  gaus_fit->SetLineWidth(1);
  gaus_fit->SetNpx(1000);

  zvtx_finder = new TF1("zvtx_finder", "pol0", -1, 20000);
  zvtx_finder->SetLineColor(2);
  zvtx_finder->SetLineWidth(1);

  if (draw_event_display)
  {
    final_fit_range_line = new TLine();

    coord_line = new TLine();

    ladder_line = new TLine();

    bkg = new TGraph(2);
  }

  if (draw_event_display || m_enable_qa)
  {
    eff_sig_range_line = new TLine();
    eff_sig_range_line->SetLineWidth(2);
    eff_sig_range_line->SetLineColor(TColor::GetColor("#A08144"));
    eff_sig_range_line->SetLineStyle(2);

    draw_text = new TLatex();
    draw_text->SetNDC();
    draw_text->SetTextSize(0.03);
  }
}

bool INTTZvtx::ProcessEvt(
    int event_i,
    std::vector<clu_info>& temp_sPH_inner_nocolumn_vec,
    std::vector<clu_info>& temp_sPH_outer_nocolumn_vec,
    std::vector<std::vector<double>>& temp_sPH_nocolumn_vec,
    std::vector<std::vector<double>>& temp_sPH_nocolumn_rz_vec,
    int NvtxMC,
    double TrigZvtxMC,
    bool PhiCheckTag,
    uint64_t bco_full,
    int centrality_bin)
{
  if (!m_initialized)
  {
    std::cout << "INTTZvtx is not initialized" << std::endl;
    exit(1);
  }

  ////////////////////////////////////////
  // initialize tree variables
  out_eID = event_i;
  bco_full_out = bco_full;
  N_cluster_inner_out = -1;
  N_cluster_outer_out = -1;

  out_ES_zvtx = -1;
  out_ES_zvtxE = -1;
  out_ES_rangeL = -1;
  out_ES_rangeR = -1;
  out_ES_N_good = -1;
  out_ES_width_density = -1;

  out_LB_Gaus_Mean_mean = -1;
  out_LB_Gaus_Mean_meanE = -1;
  out_LB_Gaus_Mean_chi2 = -1;
  out_LB_Gaus_Mean_width = -1;

  out_LB_Gaus_Width_width = -1;
  out_LB_Gaus_Width_size_width = -1;
  out_LB_Gaus_Width_offset = -1;

  out_mid_cut_Ngroup = -1;
  out_mid_cut_peak_width = -1;
  out_mid_cut_peak_ratio = -1;

  out_LB_cut_Ngroup = -1;
  out_LB_cut_peak_width = -1;
  out_LB_cut_peak_ratio = -1;

  out_LB_geo_mean = -1;
  out_good_zvtx_tag = false;

  MC_true_zvtx = TrigZvtxMC * 10.;

  out_centrality_bin = -1;

  out_N_cluster_north = 0;
  out_N_cluster_south = 0;

  m_zvtxinfo.clear();

  good_zvtx_tag = false;
  good_zvtx_tag_int = 0;

  loose_offset_peak = -9999;   // note : unit [mm]
  loose_offset_peakE = -9999;  // note : unit [mm]

  if (event_i % 1000 == 0 && print_message_opt == true)
  {
    std::cout << "In INTTZvtx class, running event : " << event_i << std::endl;
  }

  long total_NClus = temp_sPH_inner_nocolumn_vec.size() + temp_sPH_outer_nocolumn_vec.size();

  if (total_NClus < zvtx_cal_require)
  {
    if (m_enable_qa)
    {
      tree_out->Fill();
    }
    std::cout << "In INTTZvtx class, event : " << event_i << " return confirmation ntotal:" << total_NClus << std::endl;
    return false;
  }

  if (run_type == "MC" && NvtxMC != 1)
  {
    if (m_enable_qa)
    {
      tree_out->Fill();
    }
    std::cout << "In INTTZvtx class, event : " << event_i << " return Nvtx : " << NvtxMC << " Nvtx more than one " << std::endl;
    return false;
  }
  if (PhiCheckTag == false)
  {
    if (m_enable_qa)
    {
      tree_out->Fill();
    }
    std::cout << "In INTTZvtx class, event : " << event_i << " return Nvtx : " << NvtxMC << " Not full phi has hits " << std::endl;
    return false;
  }

  //--if (   temp_sPH_inner_nocolumn_vec.size() < 10
  //--    || temp_sPH_outer_nocolumn_vec.size() < 10
  //--    || total_NClus > N_clu_cut
  //--    || total_NClus < N_clu_cutl)
  if (temp_sPH_inner_nocolumn_vec.size() < 2     // originally 5
      || temp_sPH_outer_nocolumn_vec.size() < 2  // originally 5
      || total_NClus > N_clu_cut || total_NClus < N_clu_cutl)
  {
    if (m_enable_qa)
    {
      tree_out->Fill();
    }
    std::cout << (boost::format("In INTTZvtx class, event : %i, return low clu continue, NClus : %ld %lu %lu\n") % event_i % total_NClus % temp_sPH_inner_nocolumn_vec.size() % temp_sPH_outer_nocolumn_vec.size()).str();
    return false;
  }

  //--std::cout<<"--1--"<<std::endl;
  //-----------------
  // cluster pair
  // note : put the cluster into the phi map, the first bool is for the cluster usage.
  // note : false means the cluster is not used
  double Clus_InnerPhi_Offset = 0;  // note : the vertex in XY is not at zero, so the "offset" moves the offset back to the orign which is (0,0)
  double Clus_OuterPhi_Offset = 0;  // note : the vertex in XY is not at zero, so the "offset" moves the offset back to the orign which is (0,0)
  for (auto& inner_i : temp_sPH_inner_nocolumn_vec)
  {
    Clus_InnerPhi_Offset = (inner_i.y - beam_origin.second < 0)
                               ? atan2(inner_i.y - beam_origin.second,
                                       inner_i.x - beam_origin.first) *
                                         (180. / TMath::Pi()) +
                                     360
                               : atan2(inner_i.y - beam_origin.second,
                                       inner_i.x - beam_origin.first) *
                                     (180. / TMath::Pi());

    // std::cout<<"inner clu phi : "<<Clus_InnerPhi_Offset<<" origin: "<< temp_sPH_inner_nocolumn_vec[inner_i].phi <<std::endl;
    // std::cout<<" ("<<Clus_InnerPhi_Offset<<", "<< temp_sPH_inner_nocolumn_vec[inner_i].phi<<")" <<std::endl;
    //
    inner_clu_phi_map[int(Clus_InnerPhi_Offset)].push_back({false, inner_i});

    if (inner_i.z > 0)
    {
      out_N_cluster_north += 1;
    }
    else
    {
      out_N_cluster_south += 1;
    }
  }
  //--std::cout<<"--2--"<<std::endl;
  for (auto& outer_i : temp_sPH_outer_nocolumn_vec)
  {
    Clus_OuterPhi_Offset = (outer_i.y - beam_origin.second < 0)
                               ? atan2(outer_i.y - beam_origin.second,
                                       outer_i.x - beam_origin.first) *
                                         (180. / TMath::Pi()) +
                                     360
                               : atan2(outer_i.y - beam_origin.second,
                                       outer_i.x - beam_origin.first) *
                                     (180. / TMath::Pi());

    outer_clu_phi_map[int(Clus_OuterPhi_Offset)].push_back({false, outer_i});

    if (outer_i.z > 0)
    {
      out_N_cluster_north += 1;
    }
    else
    {
      out_N_cluster_south += 1;
    }
  }

  //--std::cout<<"--3--"<<std::endl;
  ////
  // tracklet reconstruction from inner and outer clusters
  int good_pair_count = 0;

  for (int inner_phi_i = 0; inner_phi_i < 360; inner_phi_i++)  // note : each phi cell (1 degree)
  {
    // note : N cluster in this phi cell
    for (unsigned int inner_phi_clu_i = 0; inner_phi_clu_i < inner_clu_phi_map[inner_phi_i].size(); inner_phi_clu_i++)
    {
      if (inner_clu_phi_map[inner_phi_i][inner_phi_clu_i].first == true)
      {
        continue;
      }

      Clus_InnerPhi_Offset = (inner_clu_phi_map[inner_phi_i][inner_phi_clu_i].second.y - beam_origin.second < 0)
                                 ? atan2(inner_clu_phi_map[inner_phi_i][inner_phi_clu_i].second.y - beam_origin.second,
                                         inner_clu_phi_map[inner_phi_i][inner_phi_clu_i].second.x - beam_origin.first) *
                                           (180. / TMath::Pi()) +
                                       360
                                 : atan2(inner_clu_phi_map[inner_phi_i][inner_phi_clu_i].second.y - beam_origin.second,
                                         inner_clu_phi_map[inner_phi_i][inner_phi_clu_i].second.x - beam_origin.first) *
                                       (180. / TMath::Pi());

      // todo: change the outer phi scan range
      // note : the outer phi index, -1, 0, 1
      for (int scan_i = -1; scan_i < 2; scan_i++)
      {
        int true_scan_i = ((inner_phi_i + scan_i) < 0)     ? 360 + (inner_phi_i + scan_i)
                          : ((inner_phi_i + scan_i) > 359) ? (inner_phi_i + scan_i) - 360
                                                           : inner_phi_i + scan_i;

        // note : N clusters in that outer phi cell
        for (unsigned int outer_phi_clu_i = 0; outer_phi_clu_i < outer_clu_phi_map[true_scan_i].size(); outer_phi_clu_i++)
        {
          if (outer_clu_phi_map[true_scan_i][outer_phi_clu_i].first == true)
          {
            continue;
          }

          Clus_OuterPhi_Offset = (outer_clu_phi_map[true_scan_i][outer_phi_clu_i].second.y - beam_origin.second < 0)
                                     ? atan2(outer_clu_phi_map[true_scan_i][outer_phi_clu_i].second.y - beam_origin.second,
                                             outer_clu_phi_map[true_scan_i][outer_phi_clu_i].second.x - beam_origin.first) *
                                               (180. / TMath::Pi()) +
                                           360
                                     : atan2(outer_clu_phi_map[true_scan_i][outer_phi_clu_i].second.y - beam_origin.second,
                                             outer_clu_phi_map[true_scan_i][outer_phi_clu_i].second.x - beam_origin.first) *
                                           (180. / TMath::Pi());

          double delta_phi = get_delta_phi(Clus_InnerPhi_Offset, Clus_OuterPhi_Offset);

          if (draw_event_display)
          {
            evt_phi_diff_1D->Fill(delta_phi);                                       // QA
            evt_phi_diff_inner_phi->Fill(Clus_InnerPhi_Offset, delta_phi);          // QA
            evt_inner_outer_phi->Fill(Clus_InnerPhi_Offset, Clus_OuterPhi_Offset);  // QA

            phi_diff_inner_phi->Fill(Clus_InnerPhi_Offset, delta_phi);  // QA
          }

          if (fabs(delta_phi) < phi_diff_cut)
          {
            double DCA_sign = calculateAngleBetweenVectors(
                outer_clu_phi_map[true_scan_i][outer_phi_clu_i].second.x, outer_clu_phi_map[true_scan_i][outer_phi_clu_i].second.y,
                inner_clu_phi_map[inner_phi_i][inner_phi_clu_i].second.x, inner_clu_phi_map[inner_phi_i][inner_phi_clu_i].second.y,
                beam_origin.first, beam_origin.second);
            if (m_enable_qa)
            {
              dca_inner_phi->Fill(Clus_InnerPhi_Offset, DCA_sign);
            }

            if (DCA_cut.first < DCA_sign && DCA_sign < DCA_cut.second)
            {
              good_pair_count += 1;
              // used_outer_check[outer_i] = 1; //note : this outer cluster was already used!

              // note : we basically transform the coordinate from cartesian to cylinder
              // note : we should set the offset first, otherwise it provides the bias
              // todo : which point should be used, DCA point or vertex xy ? Has to be studied
              std::pair<double, double> z_range_info = Get_possible_zvtx(
                  0.,  // get_radius(beam_origin.first,beam_origin.second),
                  {get_radius(inner_clu_phi_map[inner_phi_i][inner_phi_clu_i].second.x - beam_origin.first,
                              inner_clu_phi_map[inner_phi_i][inner_phi_clu_i].second.y - beam_origin.second),
                   inner_clu_phi_map[inner_phi_i][inner_phi_clu_i].second.z},  // note : unsign radius
                  {get_radius(outer_clu_phi_map[true_scan_i][outer_phi_clu_i].second.x - beam_origin.first,
                              outer_clu_phi_map[true_scan_i][outer_phi_clu_i].second.y - beam_origin.second),
                   outer_clu_phi_map[true_scan_i][outer_phi_clu_i].second.z}  // note : unsign radius
              );

              // note : try to remove some crazy background candidates. Can be a todo
              if (evt_possible_z_range.first < z_range_info.first && z_range_info.first < evt_possible_z_range.second)
              {
                N_comb.push_back(good_comb_id);
                N_comb_e.push_back(0);
                N_comb_phi.push_back(get_track_phi(Clus_InnerPhi_Offset, delta_phi));
                z_mid.push_back(z_range_info.first);
                z_range.push_back(z_range_info.second);

                evt_possible_z->Fill(z_range_info.first);  // used for calculation

                // note : fill the line_breakdwon histogram as well as a vector for the width determination
                line_breakdown(line_breakdown_hist,  // used for calculation
                               {z_range_info.first - z_range_info.second,
                                z_range_info.first + z_range_info.second});

                good_comb_id += 1;
              }
            }
          }
        }
      }  // note : end of outer clu loop
    }

  }  // note : end of inner clu loop
  //--std::cout<<"--4--"<<std::endl;

  // if (event_i == 906) {
  //     for (int hist_i = 0; hist_i < line_breakdown_hist->GetNbinsX(); hist_i++){
  //         std::cout<<line_breakdown_hist->GetBinContent(hist_i+1)<<","<<std::endl;
  //     }
  // }
  if (m_enable_qa)
  {
    //--std::cout<<"--4 0-- "<<total_NClus<<" "<< N_comb.size()<<std::endl;
    //--std::cout<<std::hex<<(long)N_track_candidate_Nclu<<std::dec<<std::endl;
    //--std::cout<<N_track_candidate_Nclu -> GetName()<<std::endl;
    N_track_candidate_Nclu->Fill(total_NClus, N_comb.size());  // QA
  }

  //--std::cout<<"--4 00--"<<std::endl;
  //-----------------------
  //
  // Fit histogram to determine Zvertex
  if (N_comb.size() > zvtx_cal_require)
  {
    //--std::cout<<"--4 1--"<<std::endl;
    N_group_info = find_Ngroup(evt_possible_z);
    N_group_info_detail = find_Ngroup(line_breakdown_hist);

    // note : first fit is for the width, so apply the constraints on the Gaussian offset
    gaus_fit->SetParameters(line_breakdown_hist->GetBinContent(line_breakdown_hist->GetMaximumBin()),
                            line_breakdown_hist->GetBinCenter(line_breakdown_hist->GetMaximumBin()),
                            40,
                            0);
    gaus_fit->SetParLimits(0, 0, 100000);  // note : size
    gaus_fit->SetParLimits(2, 5, 10000);   // note : Width
    gaus_fit->SetParLimits(3, 0, 10000);   // note : offset
    // todo : try to use single gaus to fit the distribution, and try to only fit the peak region (peak - 100 mm + peak + 100 mm)
    line_breakdown_hist->Fit(gaus_fit, "NQ", "",
                             line_breakdown_hist->GetBinCenter(line_breakdown_hist->GetMaximumBin()) - 90,
                             line_breakdown_hist->GetBinCenter(line_breakdown_hist->GetMaximumBin()) + 90);
    //----------------
    // 1st try z-vertex
    tight_offset_peak = gaus_fit->GetParameter(1);
    tight_offset_width = fabs(gaus_fit->GetParameter(2));

    double final_selection_widthU = (tight_offset_peak + tight_offset_width);
    double final_selection_widthD = (tight_offset_peak - tight_offset_width);

    double gaus_fit_offset = gaus_fit->GetParameter(3);
    gaus_fit->SetParameter(3, 0);  // note : in order to calculate the integration
    double gaus_ratio = (fabs(gaus_fit->GetParameter(0)) / gaus_fit->Integral(-600, 600)) / fabs(gaus_fit->GetParameter(2));
    gaus_fit->SetParameter(3, gaus_fit_offset);  // note : put the offset back to the function

    // note : second fit is for the peak position, therefore, loose the constraints on the Gaussian offset
    gaus_fit->SetParameters(line_breakdown_hist->GetBinContent(line_breakdown_hist->GetMaximumBin()),
                            line_breakdown_hist->GetBinCenter(line_breakdown_hist->GetMaximumBin()),
                            40,
                            0);

    gaus_fit->SetParLimits(0, 0, 100000);    // note : size
    gaus_fit->SetParLimits(2, 5, 10000);     // note : Width
    gaus_fit->SetParLimits(3, -200, 10000);  // note : offset
    // todo : try to use single gaus to fit the distribution, and try to only fit the peak region (peak - 100 mm + peak + 100 mm)
    line_breakdown_hist->Fit(gaus_fit, "NQ", "",
                             line_breakdown_hist->GetBinCenter(line_breakdown_hist->GetMaximumBin()) - 90,
                             line_breakdown_hist->GetBinCenter(line_breakdown_hist->GetMaximumBin()) + 90);

    // line_breakdown_hist -> Fit(gaus_fit, "NQ", "", N_group_info_detail[2]-10, N_group_info_detail[3]+10);

    //--std::cout<<"--5--"<<std::endl;
    //----------------
    // final z-vertex
    loose_offset_peak = gaus_fit->GetParameter(1);
    loose_offset_peakE = gaus_fit->GetParError(1);

    // additional QA below
    // note : eff sigma method, relatively sensitive to the background
    // note : use z-mid to do the effi_sig, because that line_breakdown takes too long time
    temp_event_zvtx_info = InttVertexUtil::sigmaEff_avg(z_mid, Integrate_portion);

    std::vector<double> eff_N_comb;    // QA
    std::vector<double> eff_N_comb_e;  // QA
    std::vector<double> eff_z_mid;     // QA
    std::vector<double> eff_z_range;   // QA note : eff_sig
    for (unsigned int track_i = 0; track_i < N_comb.size(); track_i++)
    {
      if (N_group_info[2] <= z_mid[track_i] && z_mid[track_i] <= N_group_info[3])
      {
        eff_N_comb.push_back(N_comb[track_i]);
        eff_N_comb_e.push_back(N_comb_e[track_i]);
        eff_z_mid.push_back(z_mid[track_i]);
        eff_z_range.push_back(z_range[track_i]);
      }

      if (draw_event_display)
      {
        if (final_selection_widthD <= z_mid[track_i] && z_mid[track_i] <= final_selection_widthU)
        {
          // note : for monitoring the the phi distribution that is used for the z vertex determination.
          // note : in principle, I expect it should be something uniform.
          evt_select_track_phi->Fill(N_comb_phi[track_i]);  // QA
        }
      }
    }

    //--std::cout<<"--6--"<<std::endl;
    if (z_range_gr != nullptr)
    {
      delete z_range_gr;
    }
    z_range_gr = new TGraphErrors(eff_N_comb.size(),
                                  &eff_N_comb[0], &eff_z_mid[0],
                                  &eff_N_comb_e[0], &eff_z_range[0]);

    z_range_gr->Fit(zvtx_finder, "NQ", "", 0, N_comb[N_comb.size() - 1]);
    double width_density_par = (double(eff_N_comb.size()) / fabs(temp_event_zvtx_info[2] - temp_event_zvtx_info[1]));

    if (zvtx_QA_width.first < tight_offset_width &&
        tight_offset_width < zvtx_QA_width.second &&
        100 < fabs(N_group_info_detail[3] - N_group_info_detail[2]) &&
        fabs(N_group_info_detail[3] - N_group_info_detail[2]) < 190)
    {
      if (N_group_info[0] < 4 &&
          N_group_info[1] >= 0.6 &&
          N_group_info_detail[0] < 7 &&
          N_group_info_detail[1] > 0.9)
      {
        good_zvtx_tag = true;
      }
      else
      {
        good_zvtx_tag = false;
      }
    }
    else
    {
      good_zvtx_tag = false;
    }

    good_zvtx_tag_int = (good_zvtx_tag == true) ? 1 : 0;  // note : so stupid way to convert bool to int...

    // note : final
    final_zvtx = loose_offset_peak;  // zvtx_finder -> GetParameter(0);

    ///////////////////////////////////

    //--std::cout<<"--7--"<<std::endl;
    if (m_enable_qa)
    {
      //
      // note : for the group information post the background cut
      peak_group_width_hist->Fill(fabs(N_group_info[3] - N_group_info[2]));
      peak_group_ratio_hist->Fill(N_group_info[1]);
      N_group_hist->Fill(N_group_info[0]);
      peak_group_detail_width_hist->Fill(fabs(N_group_info_detail[3] - N_group_info_detail[2]));
      peak_group_detail_ratio_hist->Fill(N_group_info_detail[1]);
      N_group_detail_hist->Fill(N_group_info_detail[0]);

      // note : gaus fit on the linebreak
      gaus_width_Nclu->Fill(total_NClus, tight_offset_width);
      gaus_rchi2_Nclu->Fill(total_NClus, gaus_fit->GetChisquare() / double(gaus_fit->GetNDF()));
      line_breakdown_gaus_ratio_hist->Fill(gaus_ratio);
      line_breakdown_gaus_width_hist->Fill(tight_offset_width);

      // note : the effi sig on the z_mid vector
      zvtx_evt_width_corre->Fill(total_NClus, fabs(temp_event_zvtx_info[2] - temp_event_zvtx_info[1]));
      width_density->Fill(width_density_par);
      ES_width->Fill(fabs(temp_event_zvtx_info[2] - temp_event_zvtx_info[1]));
      ES_width_ratio->Fill(total_NClus / fabs(temp_event_zvtx_info[2] - temp_event_zvtx_info[1]));

      // note : regarding to the zvtx determination, pol0 fit with error bar considered, fit range given by eff sigma method
      zvtx_evt_fitError->Fill(fabs(zvtx_finder->GetParError(0)));
      zvtx_evt_fitError_corre->Fill(total_NClus, fabs(zvtx_finder->GetParError(0)));

      if (good_zvtx_tag)
      {
        zvtx_evt_nclu_corre->Fill(total_NClus, final_zvtx);
        avg_event_zvtx->Fill(final_zvtx);
        avg_event_zvtx_vec.push_back(final_zvtx);
        if (run_type == "MC")
        {
          Z_resolution_vec.push_back(final_zvtx - (TrigZvtxMC * 10.));
          Z_resolution->Fill(final_zvtx - (TrigZvtxMC * 10.));
          Z_resolution_Nclu->Fill(total_NClus, final_zvtx - (TrigZvtxMC * 10.));
          Z_resolution_pos->Fill(TrigZvtxMC * 10., final_zvtx - (TrigZvtxMC * 10.));
          if (total_NClus > high_multi_line)
          {
            Z_resolution_pos_cut->Fill(TrigZvtxMC * 10., final_zvtx - (TrigZvtxMC * 10.));
          }
        }
      }

      final_fit_width->Fill(total_NClus, fabs(final_selection_widthU - final_selection_widthD));  // note : from LB gaus for the moment

      // note : output the root tree
      out_eID = event_i;
      N_cluster_inner_out = temp_sPH_inner_nocolumn_vec.size();
      N_cluster_outer_out = temp_sPH_outer_nocolumn_vec.size();

      out_ES_zvtx = zvtx_finder->GetParameter(0);
      out_ES_zvtxE = zvtx_finder->GetParError(0);
      out_ES_rangeL = temp_event_zvtx_info[1];
      out_ES_rangeR = temp_event_zvtx_info[2];
      out_ES_N_good = eff_N_comb.size();
      out_ES_width_density = width_density_par;

      out_LB_Gaus_Mean_mean = loose_offset_peak;
      out_LB_Gaus_Mean_meanE = gaus_fit->GetParError(1);                              // note : -> loose one
      out_LB_Gaus_Mean_chi2 = gaus_fit->GetChisquare() / double(gaus_fit->GetNDF());  // note : -> loose one
      out_LB_Gaus_Mean_width = fabs(gaus_fit->GetParameter(2));

      out_LB_Gaus_Width_width = tight_offset_width;
      out_LB_Gaus_Width_offset = gaus_fit_offset;
      out_LB_Gaus_Width_size_width = gaus_ratio;

      out_mid_cut_Ngroup = N_group_info[0];
      out_mid_cut_peak_ratio = N_group_info[1];
      out_mid_cut_peak_width = fabs(N_group_info[3] - N_group_info[2]) / 2.;

      out_LB_cut_Ngroup = N_group_info_detail[0];
      out_LB_cut_peak_ratio = N_group_info_detail[1];
      out_LB_cut_peak_width = fabs(N_group_info_detail[3] - N_group_info_detail[2]) / 2.;

      out_centrality_bin = centrality_bin;

      out_LB_geo_mean = LB_geo_mean(line_breakdown_hist,
                                    {(tight_offset_peak - tight_offset_width),
                                     (tight_offset_peak + tight_offset_width)},
                                    event_i);
      out_good_zvtx_tag = good_zvtx_tag;
      bco_full_out = bco_full;
      MC_true_zvtx = TrigZvtxMC * 10.;

      // note : if N good tracks in xy found > certain value
    }

    //--std::cout<<"--8--"<<std::endl;
    // drawing event display & QA histograms
    if (draw_event_display)
    {
      if (temp_event_xy != nullptr)
      {
        delete temp_event_xy;
      }
      temp_event_xy = new TGraph(temp_sPH_nocolumn_vec[0].size(),
                                 &temp_sPH_nocolumn_vec[0][0], &temp_sPH_nocolumn_vec[1][0]);
      temp_event_xy->SetTitle("INTT event display X-Y plane");
      temp_event_xy->GetXaxis()->SetLimits(-150, 150);
      temp_event_xy->GetYaxis()->SetRangeUser(-150, 150);
      temp_event_xy->GetXaxis()->SetTitle("X [mm]");
      temp_event_xy->GetYaxis()->SetTitle("Y [mm]");
      temp_event_xy->SetMarkerStyle(20);
      temp_event_xy->SetMarkerColor(2);
      temp_event_xy->SetMarkerSize(1);

      if (temp_event_rz != nullptr)
      {
        delete temp_event_rz;
      }
      temp_event_rz = new TGraph(temp_sPH_nocolumn_rz_vec[0].size(),
                                 &temp_sPH_nocolumn_rz_vec[0][0], &temp_sPH_nocolumn_rz_vec[1][0]);
      temp_event_rz->SetTitle("INTT event display r-Z plane");
      temp_event_rz->GetXaxis()->SetLimits(-500, 500);
      temp_event_rz->GetYaxis()->SetRangeUser(-150, 150);
      temp_event_rz->GetXaxis()->SetTitle("Z [mm]");
      temp_event_rz->GetYaxis()->SetTitle("Radius [mm]");
      temp_event_rz->SetMarkerStyle(20);
      temp_event_rz->SetMarkerColor(2);
      temp_event_rz->SetMarkerSize(1);

      // note : --------------------------------------------------------------------------------------------------------------------------

      // std::cout<<"pad "<<(long)pad_xy<<std::endl;
      // std::cout<<"bkg "<<(long)bkg<<std::endl;
      // std::cout<<"draw_text "<<(long)draw_text<<std::endl;
      pad_xy->cd();
      temp_event_xy->Draw("ap");
      // temp_event_xy -> Draw("p");
      draw_text->DrawLatex(0.2, 0.85,
                           (boost::format("eID : %i, inner Ncluster : %zu, outer Ncluster : %zu") % event_i % temp_sPH_inner_nocolumn_vec.size() % temp_sPH_outer_nocolumn_vec.size()).str().c_str());

      //--std::cout<<"--9--"<<std::endl;
      // note : --------------------------------------------------------------------------------------------------------------------------

      pad_rz->cd();
      temp_event_rz->Draw("ap");
      // eff_sig_range_line -> DrawLine(temp_event_zvtx_info[0],-150,temp_event_zvtx_info[0],150);
      coord_line->DrawLine(0, -150, 0, 150);
      coord_line->DrawLine(-500, 0, 500, 0);
      draw_text->DrawLatex(0.2, 0.85, "Negative radius : Clu_{outer} > 180^{0}");
      // draw_text -> DrawLatex(0.2, 0.81, (boost::format("EffSig avg : %.2f mm",temp_event_zvtx_info[0]));

      // note : --------------------------------------------------------------------------------------------------------------------------
      // std::cout<<"test tag 2-5"<<std::endl;
      pad_z->cd();
      if (z_range_gr_draw != nullptr)
      {
        delete z_range_gr_draw;
      }
      z_range_gr_draw = new TGraphErrors(N_comb.size(), &N_comb[0], &z_mid[0], &N_comb_e[0], &z_range[0]);
      z_range_gr_draw->GetYaxis()->SetRangeUser(-650, 650);
      z_range_gr_draw->GetXaxis()->SetTitle("Index");
      z_range_gr_draw->GetYaxis()->SetTitle("Z [mm]");
      z_range_gr_draw->SetMarkerStyle(20);
      z_range_gr_draw->Draw("ap");

      eff_sig_range_line->DrawLine(z_range_gr_draw->GetXaxis()->GetXmin(), N_group_info[2],
                                   z_range_gr_draw->GetXaxis()->GetXmax(), N_group_info[2]);
      eff_sig_range_line->DrawLine(z_range_gr_draw->GetXaxis()->GetXmin(), N_group_info[3],
                                   z_range_gr_draw->GetXaxis()->GetXmax(), N_group_info[3]);
      draw_text->DrawLatex(0.2, 0.82, (boost::format("#color[2]{Event Zvtx %.2f mm, error : #pm%.2f}") % zvtx_finder->GetParameter(0) % zvtx_finder->GetParError(0)).str().c_str());
      draw_text->DrawLatex(0.2, 0.78, (boost::format("#color[2]{Width density : %.6f}") % (width_density_par)).str().c_str());
      draw_text->DrawLatex(0.2, 0.74, (boost::format("#color[2]{Width (%.0f%%) : %.2f to %.2f mm #rightarrow %.2f mm}") % (Integrate_portion * 100.) % temp_event_zvtx_info[2] % temp_event_zvtx_info[1] % (std::fabs(temp_event_zvtx_info[2] - temp_event_zvtx_info[1]) / 2.)).str().c_str());

      // z_range_gr_draw -> Draw("p same");
      zvtx_finder->Draw("lsame");

      //--std::cout<<"--10--"<<std::endl;
      // note : --------------------------------------------------------------------------------------------------------------------------
      pad_z_hist->cd();
      evt_possible_z->SetMaximum(evt_possible_z->GetBinContent(evt_possible_z->GetMaximumBin()) * 1.4);
      evt_possible_z->Draw("hist");
      eff_sig_range_line->DrawLine(evt_possible_z->GetXaxis()->GetXmin(),
                                   evt_possible_z->GetBinContent(evt_possible_z->GetMaximumBin()) / 2.,
                                   evt_possible_z->GetXaxis()->GetXmax(),
                                   evt_possible_z->GetBinContent(evt_possible_z->GetMaximumBin()) / 2.);
      draw_text->DrawLatex(0.2, 0.82, (boost::format("N group : %.0f") % N_group_info[0]).str().c_str());
      draw_text->DrawLatex(0.2, 0.78, (boost::format("Main peak ratio : %.2f") % N_group_info[1]).str().c_str());
      draw_text->DrawLatex(0.2, 0.74, (boost::format("good z tag : %i") % good_zvtx_tag).str().c_str());

      // note : --------------------------------------------------------------------------------------------------------------------------
      pad_z_line->cd();
      line_breakdown_hist->SetMinimum(0);
      line_breakdown_hist->SetMaximum(line_breakdown_hist->GetBinContent(line_breakdown_hist->GetMaximumBin()) * 2);
      line_breakdown_hist->Draw("hist");
      gaus_fit->Draw("lsame");
      if (!good_zvtx_tag)
      {
        final_fit_range_line->DrawLine(line_breakdown_hist->GetXaxis()->GetXmin(), line_breakdown_hist->GetMaximum(),
                                       line_breakdown_hist->GetXaxis()->GetXmax(), line_breakdown_hist->GetMinimum());
      }
      final_fit_range_line->DrawLine(final_selection_widthD, line_breakdown_hist->GetMinimum(),
                                     final_selection_widthD, line_breakdown_hist->GetMaximum());
      final_fit_range_line->DrawLine(final_selection_widthU, line_breakdown_hist->GetMinimum(),
                                     final_selection_widthU, line_breakdown_hist->GetMaximum());
      draw_text->DrawLatex(0.2, 0.82, (boost::format("Gaus mean %.2f mm") % loose_offset_peak).str().c_str());
      draw_text->DrawLatex(0.2, 0.78, (boost::format("Width : %.2f mm") % tight_offset_width).str().c_str());
      draw_text->DrawLatex(0.2, 0.74, (boost::format("Reduced #chi2 : %.3f") % (gaus_fit->GetChisquare() / double(gaus_fit->GetNDF()))).str().c_str());
      draw_text->DrawLatex(0.2, 0.70, (boost::format("Norm. entry / Width : %.6f mm") % gaus_ratio).str().c_str());
      draw_text->DrawLatex(0.2, 0.66, (boost::format("LB Geo mean : %.3f mm") % out_LB_geo_mean).str().c_str());

      if (N_group_info_detail[0] != -1)
      {
        eff_sig_range_line->DrawLine(line_breakdown_hist->GetXaxis()->GetXmin(),
                                     line_breakdown_hist->GetBinContent(line_breakdown_hist->GetMaximumBin()) / 2.,
                                     line_breakdown_hist->GetXaxis()->GetXmax(),
                                     line_breakdown_hist->GetBinContent(line_breakdown_hist->GetMaximumBin()) / 2.);
        draw_text->DrawLatex(0.2, 0.62, (boost::format("N group : %.0f") % N_group_info_detail[0]).str().c_str());
        draw_text->DrawLatex(0.2, 0.58, (boost::format("Main peak ratio : %.2f") % N_group_info_detail[1]).str().c_str());
      }

      if (run_type == "MC")
      {
        draw_text->DrawLatex(0.2, 0.54, (boost::format("True MCz : %.3f mm") % (TrigZvtxMC * 10.)).str().c_str());
      }

      // note : --------------------------------------------------------------------------------------------------------------------------
      pad_phi_diff->cd();
      evt_phi_diff_inner_phi->Draw("colz0");

      // note : --------------------------------------------------------------------------------------------------------------------------
      pad_track_phi->cd();
      evt_select_track_phi->Draw("hist");

      // note : --------------------------------------------------------------------------------------------------------------------------
      pad_inner_outer_phi->cd();
      evt_inner_outer_phi->Draw("colz0");

      // note : --------------------------------------------------------------------------------------------------------------------------
      pad_phi_diff_1D->cd();
      evt_phi_diff_1D->Draw("hist");
      // std::cout<<"--11--"<<std::endl;

      if (draw_event_display && (event_i % print_rate) == 0 /*&& good_zvtx_tag == true*/)
      {
        c2->Print((boost::format("%s/temp_event_display.pdf") % out_folder_directory).str().c_str());
      }
      else if (draw_event_display && (final_zvtx > 0 || final_zvtx < -450))
      {
        std::cout << "In INTTZvtx class, event :" << event_i << " weird zvtx " << std::endl;
        c2->Print((boost::format("%s/temp_event_display.pdf") % out_folder_directory).str().c_str());
      }
      else if (draw_event_display && run_type == "MC" && fabs(final_zvtx - (TrigZvtxMC * 10.)) > 2. && total_NClus > 3000)
      {
        std::cout << "In INTTZvtx class, event :" << event_i << " High NClus, poor Z : " << fabs(final_zvtx - (TrigZvtxMC * 10.)) << std::endl;
        c2->Print((boost::format("%s/temp_event_display.pdf") % out_folder_directory).str().c_str());
      }
      else if (draw_event_display && run_type == "MC" && fabs(final_zvtx - (TrigZvtxMC * 10.)) > 2.)
      {
        std::cout << "In INTTZvtx class, event :" << event_i << " low NClus, poor Z : " << fabs(final_zvtx - (TrigZvtxMC * 10.)) << std::endl;
        c2->Print((boost::format("%s/temp_event_display.pdf") % out_folder_directory).str().c_str());
      }

      // std::cout<<"--12--"<<std::endl;
      pad_xy->Clear();
      pad_rz->Clear();
      pad_z->Clear();
      pad_z_hist->Clear();
      pad_z_line->Clear();
      pad_phi_diff->Clear();
      pad_track_phi->Clear();
      pad_inner_outer_phi->Clear();
      pad_phi_diff_1D->Clear();
      // std::cout<<"--13--"<<std::endl;

    }  // if (draw_event_display == true)

    m_zvtxinfo.zvtx = loose_offset_peak;
    m_zvtxinfo.zvtx_err = loose_offset_peakE;
    m_zvtxinfo.width = gaus_fit->GetParameter(2);
    m_zvtxinfo.chi2ndf = gaus_fit->GetChisquare() / double(gaus_fit->GetNDF());
    m_zvtxinfo.good = good_zvtx_tag;
    m_zvtxinfo.ngroup = N_group_info_detail[0];
    m_zvtxinfo.peakratio = N_group_info_detail[1];
    m_zvtxinfo.peakwidth = fabs(N_group_info_detail[3] - N_group_info_detail[2]) / 2.;

    std::cout << "chi2 : " << m_zvtxinfo.chi2ndf << " " << m_zvtxinfo.width << std::endl;

  }  // if (N_comb.size() > zvtx_cal_require)

  m_zvtxinfo.nclus = total_NClus;
  m_zvtxinfo.ntracklets = N_comb.size();

  ////////////////////////////////////////

  // std::cout<<"good pair count : "<<good_pair_count<<std::endl;
  std::cout << "evt : " << event_i << ", good pair count : " << N_comb.size() << " " << good_pair_count << std::endl;

  //--std::cout<<"--14 0--"<<std::endl;
  if (m_enable_qa)
  {
    tree_out->Fill();
  }
  //--std::cout<<"--14--"<<std::endl;

  return true;
}

void INTTZvtx::ClearEvt()
{
  if (!m_initialized)
  {
    std::cout << "INTTZvtx is not initialized" << std::endl;
    exit(1);
  }

  good_comb_id = 0;
  out_N_cluster_north = 0;
  out_N_cluster_south = 0;

  // note : ultra-stupid way to avoid the errors
  //--z_range_gr = new TGraphErrors(); z_range_gr -> Delete();

  // draw_event_display
  // z_range_gr_draw = new TGraphErrors(); z_range_gr_draw -> Delete();
  // temp_event_xy = new TGraph(); temp_event_xy -> Delete();
  // temp_event_rz = new TGraph(); temp_event_rz -> Delete();

  N_comb_phi.clear();
  N_comb.clear();
  N_comb_e.clear();
  z_mid.clear();
  z_range.clear();

  //--eff_N_comb.clear();
  //--eff_z_mid.clear();
  //--eff_N_comb_e.clear();
  //--eff_z_range.clear();
  temp_event_zvtx_info = {0, -1000, -999.99};

  N_group_info.clear();
  N_group_info_detail = {-1., -1., -1., -1.};

  evt_possible_z->Reset("ICESM");
  line_breakdown_hist->Reset("ICESM");

  if (draw_event_display)
  {
    evt_phi_diff_1D->Reset("ICESM");
    evt_inner_outer_phi->Reset("ICESM");
    evt_select_track_phi->Reset("ICESM");
    evt_phi_diff_inner_phi->Reset("ICESM");
  }

  inner_clu_phi_map.clear();
  outer_clu_phi_map.clear();
  inner_clu_phi_map = std::vector<std::vector<std::pair<bool, clu_info>>>(360);
  outer_clu_phi_map = std::vector<std::vector<std::pair<bool, clu_info>>>(360);

  // note : this is the distribution for full run
  // line_breakdown_gaus_ratio_hist -> Reset("ICESM");
}

void INTTZvtx::PrintPlots()
{
  if (!m_initialized)
  {
    std::cout << "INTTZvtx is not initialized" << std::endl;
    exit(1);
  }

  if (draw_event_display)
  {
    c2->Print((boost::format("%s/temp_event_display.pdf") % out_folder_directory).str().c_str());
    c2->Clear();
  }

  if (m_enable_qa)
  {
    c1->Clear();

    std::cout << "avg_event_zvtx_vec size : " << avg_event_zvtx_vec.size() << std::endl;

    TLatex* ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    ltx->SetTextAlign(31);

    std::string plot_text = (run_type == "MC") ? "Simulation" : "Work-in-progress";
    std::string inttlabel_text = (boost::format("#it{#bf{sPHENIX INTT}} %s") % plot_text).str();

    // note : ---------------------------------------------------------------------------------------

    c1->cd();
    std::vector<float> avg_event_zvtx_info = {0, 0, 0};
    if (avg_event_zvtx_vec.size() > 10)
    {
      avg_event_zvtx_info = InttVertexUtil::sigmaEff_avg(avg_event_zvtx_vec, Integrate_portion_final);
    }

    avg_event_zvtx->SetMinimum(0);
    avg_event_zvtx->SetMaximum(avg_event_zvtx->GetBinContent(avg_event_zvtx->GetMaximumBin()) * 1.5);
    avg_event_zvtx->Draw("hist");

    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());

    eff_sig_range_line->DrawLine(avg_event_zvtx_info[1], 0, avg_event_zvtx_info[1], avg_event_zvtx->GetMaximum());
    eff_sig_range_line->DrawLine(avg_event_zvtx_info[2], 0, avg_event_zvtx_info[2], avg_event_zvtx->GetMaximum());
    draw_text->DrawLatex(0.21, 0.87, (boost::format("EffSig min : %.2f mm") % avg_event_zvtx_info[1]).str().c_str());
    draw_text->DrawLatex(0.21, 0.83, (boost::format("EffSig max : %.2f mm") % avg_event_zvtx_info[2]).str().c_str());
    draw_text->DrawLatex(0.21, 0.79, (boost::format("EffSig avg : %.2f mm") % avg_event_zvtx_info[0]).str().c_str());
    c1->Print((boost::format("%s/avg_event_zvtx.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------

    width_density->Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/width_density.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------

    ES_width->Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/ES_width.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------

    ES_width_ratio->Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/ES_width_ratio.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------

    zvtx_evt_fitError->Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/zvtx_evt_fitError.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------
    if (run_type == "MC")
    {
      std::vector<float> Z_resolution_vec_info = InttVertexUtil::sigmaEff_avg(Z_resolution_vec, Integrate_portion_final);

      TF1* gaus_fit_2 = new TF1("gaus_fit_2", InttVertexUtil::gaus_func, evt_possible_z_range.first, evt_possible_z_range.second, 4);
      gaus_fit_2->SetLineColor(2);
      gaus_fit_2->SetLineWidth(2);
      gaus_fit_2->SetNpx(1000);

      gaus_fit_2->SetParameters(Z_resolution->GetBinContent(Z_resolution->GetMaximumBin()),
                                Z_resolution->GetBinCenter(Z_resolution->GetMaximumBin()),
                                3, 0);
      gaus_fit_2->SetParLimits(0, 0, 100000);  // note : size
      gaus_fit_2->SetParLimits(2, 0, 10000);   // note : Width
      gaus_fit_2->SetParLimits(3, 0, 10000);   // note : offset
      Z_resolution->Fit(gaus_fit_2, "NQ", "",
                        Z_resolution->GetBinCenter(Z_resolution->GetMaximumBin()) - (2 * Z_resolution->GetStdDev()),
                        Z_resolution->GetBinCenter(Z_resolution->GetMaximumBin()) + (2 * Z_resolution->GetStdDev()));

      Z_resolution->Draw("hist");
      gaus_fit_2->SetRange(gaus_fit_2->GetParameter(1) - gaus_fit_2->GetParameter(2) * 2.5,
                           gaus_fit_2->GetParameter(1) + gaus_fit_2->GetParameter(2) * 2.5);
      gaus_fit_2->Draw("lsame");
      eff_sig_range_line->DrawLine(Z_resolution_vec_info[1], 0, Z_resolution_vec_info[1], Z_resolution->GetMaximum());
      eff_sig_range_line->DrawLine(Z_resolution_vec_info[2], 0, Z_resolution_vec_info[2], Z_resolution->GetMaximum());
      draw_text->DrawLatex(0.21, 0.87, (boost::format("EffSig min : %.2f mm") % Z_resolution_vec_info[1]).str().c_str());
      draw_text->DrawLatex(0.21, 0.83, (boost::format("EffSig max : %.2f mm") % Z_resolution_vec_info[2]).str().c_str());
      draw_text->DrawLatex(0.21, 0.79, (boost::format("EffSig avg : %.2f mm") % Z_resolution_vec_info[0]).str().c_str());
      draw_text->DrawLatex(0.21, 0.71, (boost::format("Gaus mean  : %.2f mm") % gaus_fit_2->GetParameter(1)).str().c_str());
      draw_text->DrawLatex(0.21, 0.67, (boost::format("Gaus width : %.2f mm") % fabs(gaus_fit_2->GetParameter(2))).str().c_str());
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
      c1->Print((boost::format("%s/Z_resolution.pdf") % out_folder_directory).str().c_str());
      c1->Clear();

      MC_z_diff_peak = gaus_fit_2->GetParameter(1);
      MC_z_diff_width = fabs(gaus_fit_2->GetParameter(2));

      delete gaus_fit_2;
    }

    // note : ---------------------------------------------------------------------------------------

    if (run_type == "MC")
    {
      Z_resolution_Nclu->Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
      c1->Print((boost::format("%s/Z_resolution_Nclu.pdf") % out_folder_directory).str().c_str());
      c1->Clear();
    }
    // note : ---------------------------------------------------------------------------------------

    if (run_type == "MC")
    {
      Z_resolution_pos->Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
      c1->Print((boost::format("%s/Z_resolution_pos.pdf") % out_folder_directory).str().c_str());
      c1->Clear();
    }
    // note : ---------------------------------------------------------------------------------------

    if (run_type == "MC")
    {
      Z_resolution_pos_cut->Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
      c1->Print((boost::format("%s/Z_resolution_pos_cut.pdf") % out_folder_directory).str().c_str());
      c1->Clear();
    }
    // note : ---------------------------------------------------------------------------------------

    zvtx_evt_fitError_corre->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/zvtx_evt_fitError_corre.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------

    zvtx_evt_nclu_corre->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/zvtx_evt_nclu_corre.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------

    zvtx_evt_width_corre->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/zvtx_evt_width_corre.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------

    gaus_width_Nclu->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/gaus_width_Nclu.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------

    gaus_rchi2_Nclu->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/gaus_rchi2_Nclu.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------

    final_fit_width->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/final_fit_width.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------

    N_track_candidate_Nclu->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/N_track_candidate_Nclu.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------
    peak_group_width_hist->Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/peak_group_width_hist.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------
    peak_group_ratio_hist->Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/peak_group_ratio_hist.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------
    N_group_hist->Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/N_group_hist.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------
    peak_group_detail_width_hist->Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/peak_group_detail_width_hist.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------
    peak_group_detail_ratio_hist->Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/peak_group_detail_ratio_hist.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------
    N_group_detail_hist->Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());
    c1->Print((boost::format("%s/N_group_detail_hist.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------
    line_breakdown_gaus_ratio_hist->Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());

    c1->Print((boost::format("%s/line_breakdown_gaus_ratio_hist.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : ---------------------------------------------------------------------------------------
    gaus_fit->SetParameters(line_breakdown_gaus_width_hist->GetBinContent(line_breakdown_gaus_width_hist->GetMaximumBin()), line_breakdown_gaus_width_hist->GetBinCenter(line_breakdown_gaus_width_hist->GetMaximumBin()), 8, 0);
    gaus_fit->SetParLimits(0, 0, 100000);  // note : size
    gaus_fit->SetParLimits(2, 1, 10000);   // note : Width
    gaus_fit->SetParLimits(3, 0, 10000);   // note : offset

    // todo : try to use single gaus to fit the distribution, and try to only fit the peak region (peak - 100 mm + peak + 100 mm)
    line_breakdown_gaus_width_hist->Fit(gaus_fit, "NQ", "", line_breakdown_gaus_width_hist->GetBinCenter(line_breakdown_gaus_width_hist->GetMaximumBin()) - 100, line_breakdown_gaus_width_hist->GetBinCenter(line_breakdown_gaus_width_hist->GetMaximumBin()) + 100);
    line_breakdown_gaus_width_hist->Draw("hist");
    gaus_fit->Draw("lsame");

    draw_text->DrawLatex(0.2, 0.82, (boost::format("Gaus mean %.2f mm") % gaus_fit->GetParameter(1)).str().c_str());
    draw_text->DrawLatex(0.2, 0.78, (boost::format("Width : %.2f mm") % fabs(gaus_fit->GetParameter(2))).str().c_str());
    draw_text->DrawLatex(0.2, 0.74, (boost::format("Reduced #chi2 : %.3f") % (gaus_fit->GetChisquare() / double(gaus_fit->GetNDF()))).str().c_str());
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, inttlabel_text.c_str());

    c1->Print((boost::format("%s/line_breakdown_gaus_width_hist.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    delete ltx;
  }
}

void INTTZvtx::EndRun()
{
  if (!m_initialized)
  {
    std::cout << "INTTZvtx is not initialized" << std::endl;
    exit(1);
  }

  if (m_enable_qa)
  {
    out_file->cd();
    tree_out->SetDirectory(out_file);
    tree_out->Write("", TObject::kOverwrite);

    for (auto& itr : m_v_qahist)
    {
      itr->Write();
    }

    out_file->Close();
  }
}

double INTTZvtx::GetZdiffPeakMC()
{
  if (run_type == "MC" && MC_z_diff_peak != -777.)
  {
    return MC_z_diff_peak;
  }
  else
  {
    std::cout << "In INTTZvtx. Are you playing with data? The MC_z_diff_peak wasn't assigned, the value is still -777. Pleak check" << std::endl;
    return -777.;
  }
}

double INTTZvtx::GetZdiffWidthMC()
{
  if (run_type == "MC" && MC_z_diff_width != -777.)
  {
    return MC_z_diff_width;
  }
  else
  {
    std::cout << "In INTTZvtx. Are you playing with data? The MC_z_diff_width wasn't assigned, the value is still -777. Pleak check" << std::endl;
    return -777.;
  }
}

std::vector<double> INTTZvtx::GetEvtZPeak()
{
  return {good_zvtx_tag_int, loose_offset_peak, loose_offset_peakE};
}

double INTTZvtx::get_radius(double x, double y)
{
  return sqrt(pow(x, 2) + pow(y, 2));
}

// note : Function to calculate the angle between two vectors in degrees using the cross product
double INTTZvtx::calculateAngleBetweenVectors(double x1, double y1, double x2, double y2, double targetX, double targetY)
{
  // Calculate the vectors vector_1 (point_1 to point_2) and vector_2 (point_1 to target)
  double vector1X = x2 - x1;
  double vector1Y = y2 - y1;

  double vector2X = targetX - x1;
  double vector2Y = targetY - y1;

  // Calculate the cross product of vector_1 and vector_2 (z-component)
  double crossProduct = vector1X * vector2Y - vector1Y * vector2X;

  // std::cout<<" crossProduct : "<<crossProduct<<std::endl;

  // Calculate the magnitudes of vector_1 and vector_2
  double magnitude1 = std::sqrt(vector1X * vector1X + vector1Y * vector1Y);
  double magnitude2 = std::sqrt(vector2X * vector2X + vector2Y * vector2Y);

  // Calculate the angle in radians using the inverse tangent of the cross product and dot product
  //--double dotProduct = vector1X * vector2X + vector1Y * vector2Y;

  //--double angleInRadians = std::atan2(std::abs(crossProduct), dotProduct);
  // Convert the angle from radians to degrees and return it
  //--double angleInDegrees = angleInRadians * 180.0 / M_PI;

  double angleInRadians_new = std::asin(crossProduct / (magnitude1 * magnitude2));
  //--double angleInDegrees_new = angleInRadians_new * 180.0 / M_PI;

  // std::cout<<"angle : "<<angleInDegrees_new<<std::endl;

  double DCA_distance = sin(angleInRadians_new) * magnitude2;

  return DCA_distance;
}

double INTTZvtx::Get_extrapolation(double given_y, double p0x, double p0y, double p1x, double p1y)  // note : x : z, y : r
{
  if (fabs(p0x - p1x) < 0.00001)
  {  // note : the line is vertical (if z is along the x axis)
    return p0x;
  }
  else
  {
    double slope = (p1y - p0y) / (p1x - p0x);
    double yIntercept = p0y - slope * p0x;
    double xCoordinate = (given_y - yIntercept) / slope;
    return xCoordinate;
  }
}

std::pair<double, double> INTTZvtx::Get_possible_zvtx(double rvtx, std::vector<double> p0, std::vector<double> p1)  // note : inner p0, outer p1, vector {r,z}, -> {y,x}
{
  std::vector<double> p0_z_edge = {(fabs(p0[1]) < 130) ? p0[1] - 8. : p0[1] - 10., (fabs(p0[1]) < 130) ? p0[1] + 8. : p0[1] + 10.};  // note : vector {left edge, right edge}
  std::vector<double> p1_z_edge = {(fabs(p1[1]) < 130) ? p1[1] - 8. : p1[1] - 10., (fabs(p1[1]) < 130) ? p1[1] + 8. : p1[1] + 10.};  // note : vector {left edge, right edge}

  double edge_first = Get_extrapolation(rvtx, p0_z_edge[0], p0[0], p1_z_edge[1], p1[0]);
  double edge_second = Get_extrapolation(rvtx, p0_z_edge[1], p0[0], p1_z_edge[0], p1[0]);

  double mid_point = (edge_first + edge_second) / 2.;
  double possible_width = fabs(edge_first - edge_second) / 2.;

  return {mid_point, possible_width};  // note : first : mid point, second : width
}

void INTTZvtx::line_breakdown(TH1* hist_in, std::pair<double, double> line_range)
{
  int first_bin = int((line_range.first - hist_in->GetXaxis()->GetXmin()) / hist_in->GetBinWidth(1)) + 1;
  int last_bin = int((line_range.second - hist_in->GetXaxis()->GetXmin()) / hist_in->GetBinWidth(1)) + 1;

  if (first_bin < 1)
  {
    first_bin = 0;
  }
  else if (first_bin > hist_in->GetNbinsX())
  {
    first_bin = hist_in->GetNbinsX() + 1;
  }

  if (last_bin < 1)
  {
    last_bin = 0;
  }
  else if (last_bin > hist_in->GetNbinsX())
  {
    last_bin = hist_in->GetNbinsX() + 1;
  }

  // std::cout<<"Digitize the bin : "<<first_bin<<" "<<last_bin<<std::endl;

  // note : if first:last = (0:0) or (N+1:N+1) -> the subtraction of them euqals to zero.
  for (int i = 0; i < (last_bin - first_bin) + 1; i++)
  {
    hist_in->SetBinContent(first_bin + i, hist_in->GetBinContent(first_bin + i) + 1);
  }
}

// note : search_range : should be the gaus fit range
double INTTZvtx::LB_geo_mean(TH1* hist_in, std::pair<double, double> search_range, int /*event_i*/)
{
  int Highest_bin_index = hist_in->GetMaximumBin();
  double Highest_bin_center = hist_in->GetBinCenter(Highest_bin_index);
  double Highest_bin_content = hist_in->GetBinContent(Highest_bin_index);
  if (Highest_bin_center < search_range.first || search_range.second < Highest_bin_center)
  {
    // std::cout<<"In INTTZvtx class, event : "<<event_i<<", interesting event, different group was fit"<<std::endl;
    return -999.;
  }

  std::vector<float> same_height_bin;
  same_height_bin.clear();
  same_height_bin.push_back(Highest_bin_center);

  int search_index = 1;
  while ((hist_in->GetBinCenter(Highest_bin_index + search_index)) < search_range.second)
  {
    if (hist_in->GetBinContent(Highest_bin_index + search_index) == Highest_bin_content)
    {
      same_height_bin.push_back(hist_in->GetBinCenter(Highest_bin_index + search_index));
    }
    search_index++;
  }

  // std::cout<<"test, search_index right : "<<search_index<<std::endl;

  search_index = 1;
  while (search_range.first < (hist_in->GetBinCenter(Highest_bin_index - search_index)))
  {
    if (hist_in->GetBinContent(Highest_bin_index - search_index) == Highest_bin_content)
    {
      same_height_bin.push_back(hist_in->GetBinCenter(Highest_bin_index - search_index));
    }
    search_index++;
  }

  // std::cout<<"test, search_index left : "<<search_index<<std::endl;
  // std::cout<<"test, same_height_bin.size(): "<<same_height_bin.size()<<std::endl;

  return accumulate(same_height_bin.begin(), same_height_bin.end(), 0.0) / double(same_height_bin.size());
}

//  note : N group, group size, group tag, group width ?  // note : {group size, group entry, group tag, group widthL, group widthR}
// note : {N_group, ratio (if two), peak widthL, peak widthR}
std::vector<double> INTTZvtx::find_Ngroup(TH1* hist_in)
{
  double Highest_bin_Content = hist_in->GetBinContent(hist_in->GetMaximumBin());
  double Highest_bin_Center = hist_in->GetBinCenter(hist_in->GetMaximumBin());

  int group_Nbin = 0;
  int peak_group_ID = 0;  // =0 added by TH 20240418
  double group_entry = 0;
  double peak_group_ratio;
  std::vector<int> group_Nbin_vec;
  group_Nbin_vec.clear();
  std::vector<double> group_entry_vec;
  group_entry_vec.clear();
  std::vector<double> group_widthL_vec;
  group_widthL_vec.clear();
  std::vector<double> group_widthR_vec;
  group_widthR_vec.clear();

  for (int i = 0; i < hist_in->GetNbinsX(); i++)
  {
    // todo : the background rejection is here : Highest_bin_Content/2. for the time being
    double bin_content = (hist_in->GetBinContent(i + 1) <= Highest_bin_Content / 2.) ? 0. : (hist_in->GetBinContent(i + 1) - Highest_bin_Content / 2.);

    if (bin_content != 0)
    {
      if (group_Nbin == 0)
      {
        group_widthL_vec.push_back(hist_in->GetBinCenter(i + 1) - (hist_in->GetBinWidth(i + 1) / 2.));
      }

      group_Nbin += 1;
      group_entry += bin_content;
    }
    else if (bin_content == 0 && group_Nbin != 0)
    {
      group_widthR_vec.push_back(hist_in->GetBinCenter(i + 1) - (hist_in->GetBinWidth(i + 1) / 2.));
      group_Nbin_vec.push_back(group_Nbin);
      group_entry_vec.push_back(group_entry);
      group_Nbin = 0;
      group_entry = 0;
    }
  }
  if (group_Nbin != 0)
  {
    group_Nbin_vec.push_back(group_Nbin);
    group_entry_vec.push_back(group_entry);
    group_widthR_vec.push_back(hist_in->GetXaxis()->GetXmax());
  }  // note : the last group at the edge

  // note : find the peak group
  for (unsigned int i = 0; i < group_Nbin_vec.size(); i++)
  {
    if (group_widthL_vec[i] < Highest_bin_Center && Highest_bin_Center < group_widthR_vec[i])
    {
      peak_group_ID = i;
      break;
    }
  }

  if (group_entry_vec.size() > 0)
  {
    peak_group_ratio = group_entry_vec[peak_group_ID] / (accumulate(group_entry_vec.begin(), group_entry_vec.end(), 0.0));
  }
  else
  {
    peak_group_ratio = 0.0;
  }

  // for (int i = 0; i < group_Nbin_vec.size(); i++)
  // {
  //     std::cout<<" "<<std::endl;
  //     std::cout<<"group size : "<<group_Nbin_vec[i]<<std::endl;
  //     std::cout<<"group entry : "<<group_entry_vec[i]<<std::endl;
  //     std::cout<<group_widthL_vec[i]<<" "<<group_widthR_vec[i]<<std::endl;
  // }

  // std::cout<<" "<<std::endl;
  // std::cout<<"N group : "<<group_Nbin_vec.size()<<std::endl;
  // std::cout<<"Peak group ID : "<<peak_group_ID<<std::endl;
  // std::cout<<"peak group width : "<<group_widthL_vec[peak_group_ID]<<" "<<group_widthR_vec[peak_group_ID]<<std::endl;
  // std::cout<<"ratio : "<<peak_group_ratio<<std::endl;

  // for the case that all bin content in the for statemene above is 0
  if (int(group_widthL_vec.size()) <= peak_group_ID || int(group_widthR_vec.size()) <= peak_group_ID)
  {  // added by Genki (Jan 2025)
    return {double(group_Nbin_vec.size()), peak_group_ratio, -9999, -9999};
  }

  // note : {N_group, ratio (if two), peak widthL, peak widthR}
  return {double(group_Nbin_vec.size()), peak_group_ratio, group_widthL_vec[peak_group_ID], group_widthR_vec[peak_group_ID]};
}

double INTTZvtx::get_delta_phi(double angle_1, double angle_2)
{
  std::vector<double> vec_abs = {fabs(angle_1 - angle_2), fabs(angle_1 - angle_2 + 360), fabs(angle_1 - angle_2 - 360)};
  std::vector<double> vec = {(angle_1 - angle_2), (angle_1 - angle_2 + 360), (angle_1 - angle_2 - 360)};
  return vec[std::distance(vec_abs.begin(), std::min_element(vec_abs.begin(), vec_abs.end()))];
}

double INTTZvtx::get_track_phi(double inner_clu_phi_in, double delta_phi_in)
{
  double track_phi = inner_clu_phi_in - (delta_phi_in / 2.);
  if (track_phi < 0)
  {
    track_phi += 360;
  }
  else if (track_phi > 360)
  {
    track_phi -= 360;
  }
  else if (track_phi == 360)
  {
    track_phi = 0;
  }
  // else {track_phi = track_phi;}
  return track_phi;
}
