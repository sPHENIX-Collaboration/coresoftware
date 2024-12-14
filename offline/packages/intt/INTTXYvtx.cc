#include "INTTXYvtx.h"

#include <TF1.h>

#include <boost/format.hpp>

#include <cmath>

double cos_func(double* x, double* par)
{
  return -1 * par[0] * cos(par[1] * (x[0] - par[2])) + par[3];
}

INTTXYvtx::INTTXYvtx(const std::string& runType,
                     const std::string& outFolderDirectory,
                     std::pair<double, double> beamOrigin,
                     double phiDiffCut,
                     std::pair<double, double> DCACut,
                     int NCluCutl,
                     int NCluCut,
                     double angleDiffNew_l,
                     double angleDiffNew_r,
                     double peekCut,
                     bool printMessageOpt)
  : run_type(runType)
  , out_folder_directory(outFolderDirectory)
  , beam_origin(beamOrigin)
  , phi_diff_cut(phiDiffCut)
  , DCA_cut(DCACut)
  , N_clu_cutl(NCluCutl)
  , N_clu_cut(NCluCut)
  , angle_diff_new_l(angleDiffNew_l)
  , angle_diff_new_r(angleDiffNew_r)
  , peek(peekCut)
  , print_message_opt(printMessageOpt)
{
  gErrorIgnoreLevel = kWarning;  // note : To not print the "print plot info."

  // Init();
  plot_text = (run_type == "MC") ? "Simulation" : "Work-in-progress";

  cluster_pair_vec.clear();

  current_vtxX = beam_origin.first;
  current_vtxY = beam_origin.second;
}

INTTXYvtx::~INTTXYvtx()
{
  // InitHist
  for (auto& itr : m_v_hist)
  {
    //--std::cout<<"del : "<<itr->GetTitle()<<" "<<std::hex<<(long)itr<<std::hex<<std::endl;
    delete itr;
  }
  //--for(auto& itr: m_v_hist){
  //--  std::cout<<"after del : "<<std::hex<<(long)itr<<std::hex<<std::endl;
  //--}

  // InitGraph
  delete angle_diff_inner_phi_peak_profile_graph;
  delete angle_diff_outer_phi_peak_profile_graph;
  delete DCA_distance_inner_phi_peak_profile_graph;
  delete DCA_distance_outer_phi_peak_profile_graph;

  // created in FillLine_FindVertex
  delete xy_hist;
  delete xy_hist_bkgrm;

  // InitRest
  delete cos_fit;
  delete gaus_fit;
  delete horizontal_fit_inner;
  delete horizontal_fit_angle_diff_inner;
  delete horizontal_fit_outer;
  delete horizontal_fit_angle_diff_outer;

  delete c1;
  delete ltx;
  delete draw_text;
}

void INTTXYvtx::Init()
{
  if (m_enable_drawhist || m_savehist)
  {
    if (!std::filesystem::exists(out_folder_directory))
    {
      std::filesystem::create_directory(out_folder_directory);
    }
  }

  InitHist();
  InitGraph();
  InitRest();

  //--    InitTreeOut();

  m_initialized = true;
}

void INTTXYvtx::InitGraph()
{
  angle_diff_inner_phi_peak_profile_graph = new TGraph();
  angle_diff_outer_phi_peak_profile_graph = new TGraph();
  DCA_distance_inner_phi_peak_profile_graph = new TGraph();
  DCA_distance_outer_phi_peak_profile_graph = new TGraph();
}

void INTTXYvtx::InitRest()
{
  horizontal_fit_inner = new TF1("horizontal_fit_inner", "pol0", 0, 360);
  horizontal_fit_inner->SetLineWidth(2);
  horizontal_fit_inner->SetLineColor(2);

  horizontal_fit_angle_diff_inner = new TF1("horizontal_fit_angle_diff_inner", "pol0", 0, 360);
  horizontal_fit_angle_diff_inner->SetLineWidth(2);
  horizontal_fit_angle_diff_inner->SetLineColor(2);

  horizontal_fit_outer = new TF1("horizontal_fit_outer", "pol0", 0, 360);
  horizontal_fit_outer->SetLineWidth(2);
  horizontal_fit_outer->SetLineColor(2);

  horizontal_fit_angle_diff_outer = new TF1("horizontal_fit_angle_diff_outer", "pol0", 0, 360);
  horizontal_fit_angle_diff_outer->SetLineWidth(2);
  horizontal_fit_angle_diff_outer->SetLineColor(2);

  if (m_enable_qa)
  {
    cos_fit = new TF1("cos_fit", cos_func, 0, 360, 4);
    cos_fit->SetParNames("[A]", "[B]", "[C]", "[D]");
    cos_fit->SetLineColor(2);

    gaus_fit = new TF1("gaus_fit", InttVertexUtil::gaus_func, 0, 360, 4);
    gaus_fit->SetLineColor(4);
    gaus_fit->SetLineWidth(1);
    gaus_fit->SetParNames("size", "mean", "width", "offset");
    gaus_fit->SetNpx(1000);
  }

  if (m_enable_drawhist)
  {
    c1 = new TCanvas("", "", 950, 800);
    c1->cd();

    ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    ltx->SetTextAlign(31);

    draw_text = new TLatex();
    draw_text->SetNDC();
    draw_text->SetTextSize(0.03);
  }
}

void INTTXYvtx::InitHist()
{
  ///////////////////////////////////////
  // QA histograms in process_evt
  if (m_enable_qa)
  {
    N_cluster_correlation = new TH2F("N_cluster_correlation", "N_cluster_correlation", 100, 0, 4000, 100, 0, 4000);
    N_cluster_correlation->SetStats(false);
    N_cluster_correlation->GetXaxis()->SetTitle("inner N_cluster");
    N_cluster_correlation->GetYaxis()->SetTitle("Outer N_cluster");
    N_cluster_correlation->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(N_cluster_correlation);

    N_cluster_correlation_close = new TH2F("N_cluster_correlation_close", "N_cluster_correlation_close", 100, 0, 500, 100, 0, 500);
    N_cluster_correlation_close->SetStats(false);
    N_cluster_correlation_close->GetXaxis()->SetTitle("inner N_cluster");
    N_cluster_correlation_close->GetYaxis()->SetTitle("Outer N_cluster");
    N_cluster_correlation_close->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(N_cluster_correlation_close);

    inner_pos_xy = new TH2F("inner_pos_xy", "inner_pos_xy", 360, -100, 100, 360, -100, 100);
    inner_pos_xy->SetStats(false);
    inner_pos_xy->GetXaxis()->SetTitle("X axis [mm]");
    inner_pos_xy->GetYaxis()->SetTitle("Y axis [mm]");
    inner_pos_xy->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(inner_pos_xy);

    outer_pos_xy = new TH2F("outer_pos_xy", "outer_pos_xy", 360, -150, 150, 360, -150, 150);
    outer_pos_xy->SetStats(false);
    outer_pos_xy->GetXaxis()->SetTitle("X axis [mm]");
    outer_pos_xy->GetYaxis()->SetTitle("Y axis [mm]");
    outer_pos_xy->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(outer_pos_xy);

    // inner_outer_pos_xy = new TH2F("inner_outer_pos_xy","inner_outer_pos_xy",360,-150,150,360,-150,150);
    inner_outer_pos_xy = new TH2F("inner_outer_pos_xy", "inner_outer_pos_xy", 360, -150, 150, 360, -150, 150);
    inner_outer_pos_xy->SetStats(false);
    inner_outer_pos_xy->GetXaxis()->SetTitle("X axis [mm]");
    inner_outer_pos_xy->GetYaxis()->SetTitle("Y axis [mm]");
    inner_outer_pos_xy->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(inner_outer_pos_xy);
  }

  ///////////////////////////////////////
  // histograms for quadrant method
  DCA_distance_inner_phi = new TH2F("DCA_distance_inner_phi", "DCA_distance_inner_phi", 100, 0, 360, 100, -10, 10);
  DCA_distance_inner_phi->SetStats(false);
  DCA_distance_inner_phi->GetXaxis()->SetTitle("Inner phi [degree]");
  DCA_distance_inner_phi->GetYaxis()->SetTitle("DCA [mm]");
  DCA_distance_inner_phi->GetXaxis()->SetNdivisions(505);
  m_v_hist.push_back(DCA_distance_inner_phi);

  angle_diff_inner_phi = new TH2F("angle_diff_inner_phi", "angle_diff_inner_phi", 361, 0, 361, 100, -1.5, 1.5);
  angle_diff_inner_phi->SetStats(false);
  angle_diff_inner_phi->GetXaxis()->SetTitle("Inner phi [degree]");
  angle_diff_inner_phi->GetYaxis()->SetTitle("Inner - Outer [degree]");
  angle_diff_inner_phi->GetXaxis()->SetNdivisions(505);
  m_v_hist.push_back(angle_diff_inner_phi);

  ///////////////////////////////////////
  // note : it's for the geometry correction
  angle_diff_inner_phi_peak_final = new TH2F("angle_diff_inner_phi_peak_final",
                                             "angle_diff_inner_phi_peak_final", 361, 0, 361, 100, -1.5, 1.5);
  angle_diff_inner_phi_peak_final->SetStats(false);
  angle_diff_inner_phi_peak_final->GetXaxis()->SetTitle("Inner phi [degree]");
  angle_diff_inner_phi_peak_final->GetYaxis()->SetTitle("Inner - Outer [degree]");
  angle_diff_inner_phi_peak_final->GetXaxis()->SetNdivisions(505);
  m_v_hist.push_back(angle_diff_inner_phi_peak_final);

  DCA_distance_inner_phi_peak_final = new TH2F("DCA_distance_inner_phi_peak_final",
                                               "DCA_distance_inner_phi_peak_final", 100, 0, 360, 100, -10, 10);
  DCA_distance_inner_phi_peak_final->SetStats(false);
  DCA_distance_inner_phi_peak_final->GetXaxis()->SetTitle("Inner phi [degree]");
  DCA_distance_inner_phi_peak_final->GetYaxis()->SetTitle("DCA [mm]");
  DCA_distance_inner_phi_peak_final->GetXaxis()->SetNdivisions(505);
  m_v_hist.push_back(DCA_distance_inner_phi_peak_final);

  // QA
  angle_diff = new TH1F("angle_diff", "angle_diff", 100, 0, 5);
  angle_diff->SetStats(false);
  angle_diff->GetXaxis()->SetTitle("|Inner - Outer| [degree]");
  angle_diff->GetYaxis()->SetTitle("Entry");
  angle_diff->GetXaxis()->SetNdivisions(505);
  m_v_hist.push_back(angle_diff);

  angle_diff_new = new TH1F("angle_diff_new", "angle_diff_new", 100, angle_diff_new_l, angle_diff_new_r);
  angle_diff_new->SetStats(false);
  angle_diff_new->GetXaxis()->SetTitle("|Inner - Outer| [degree]");
  angle_diff_new->GetYaxis()->SetTitle("Entry");
  angle_diff_new->GetXaxis()->SetNdivisions(505);
  m_v_hist.push_back(angle_diff_new);

  DCA_distance = new TH1F("DCA_distance", "DCA_distance", 100, -10, 10);
  DCA_distance->SetStats(false);
  DCA_distance->GetXaxis()->SetTitle("DCA [mm]");
  DCA_distance->GetYaxis()->SetTitle("Entry");
  DCA_distance->GetXaxis()->SetNdivisions(505);
  m_v_hist.push_back(DCA_distance);

  if (m_enable_qa)
  {
    DCA_distance_outer_phi = new TH2F("DCA_distance_outer_phi", "DCA_distance_outer_phi", 100, 0, 360, 100, -10, 10);
    DCA_distance_outer_phi->SetStats(false);
    DCA_distance_outer_phi->GetXaxis()->SetTitle("Outer phi [degree]");
    DCA_distance_outer_phi->GetYaxis()->SetTitle("DCA [mm]");
    DCA_distance_outer_phi->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(DCA_distance_outer_phi);

    angle_diff_outer_phi = new TH2F("angle_diff_outer_phi", "angle_diff_outer_phi", 361, 0, 361, 100, -1.5, 1.5);
    angle_diff_outer_phi->SetStats(false);
    angle_diff_outer_phi->GetXaxis()->SetTitle("Outer phi [degree]");
    angle_diff_outer_phi->GetYaxis()->SetTitle("Inner - Outer [degree]");
    angle_diff_outer_phi->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(angle_diff_outer_phi);

    angle_correlation = new TH2F("angle_correlation", "angle_correlation", 361, 0, 361, 361, 0, 361);
    angle_correlation->SetStats(false);
    angle_correlation->GetXaxis()->SetTitle("Inner Phi [degree]");
    angle_correlation->GetYaxis()->SetTitle("Outer Phi [degree]");
    angle_correlation->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(angle_correlation);

    angle_diff_DCA_dist = new TH2F("angle_diff_DCA_dist", "angle_diff_DCA_dist", 100, -1.5, 1.5, 100, -3.5, 3.5);
    angle_diff_DCA_dist->SetStats(false);
    angle_diff_DCA_dist->GetXaxis()->SetTitle("Inner - Outer [degree]");
    angle_diff_DCA_dist->GetYaxis()->SetTitle("DCA distance [mm]");
    angle_diff_DCA_dist->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(angle_diff_DCA_dist);

    DCA_point = new TH2F("DCA_point", "DCA_point", 100, -10, 10, 100, -10, 10);
    DCA_point->SetStats(false);
    DCA_point->GetXaxis()->SetTitle("X pos [mm]");
    DCA_point->GetYaxis()->SetTitle("Y pos [mm]");
    DCA_point->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(DCA_point);

    DCA_distance_inner_X = new TH2F("DCA_distance_inner_X", "DCA_distance_inner_X", 100, -100, 100, 100, -10, 10);
    DCA_distance_inner_X->SetStats(false);
    DCA_distance_inner_X->GetXaxis()->SetTitle("Inner cluster X [mm]");
    DCA_distance_inner_X->GetYaxis()->SetTitle("DCA [mm]");
    DCA_distance_inner_X->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(DCA_distance_inner_X);

    DCA_distance_inner_Y = new TH2F("DCA_distance_inner_Y", "DCA_distance_inner_Y", 100, -100, 100, 100, -10, 10);
    DCA_distance_inner_Y->SetStats(false);
    DCA_distance_inner_Y->GetXaxis()->SetTitle("Inner cluster Y [mm]");
    DCA_distance_inner_Y->GetYaxis()->SetTitle("DCA [mm]");
    DCA_distance_inner_Y->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(DCA_distance_inner_Y);

    DCA_distance_outer_X = new TH2F("DCA_distance_outer_X", "DCA_distance_outer_X", 100, -130, 130, 100, -10, 10);
    DCA_distance_outer_X->SetStats(false);
    DCA_distance_outer_X->GetXaxis()->SetTitle("Outer cluster X [mm]");
    DCA_distance_outer_X->GetYaxis()->SetTitle("DCA [mm]");
    DCA_distance_outer_X->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(DCA_distance_outer_X);

    DCA_distance_outer_Y = new TH2F("DCA_distance_outer_Y", "DCA_distance_outer_Y", 100, -130, 130, 100, -10, 10);
    DCA_distance_outer_Y->SetStats(false);
    DCA_distance_outer_Y->GetXaxis()->SetTitle("Outer cluster Y [mm]");
    DCA_distance_outer_Y->GetYaxis()->SetTitle("DCA [mm]");
    DCA_distance_outer_Y->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(DCA_distance_outer_Y);

    ///////////////////////////////////////
    // note : it's for the geometry correction
    angle_diff_outer_phi_peak_final = new TH2F("angle_diff_outer_phi_peak_final",
                                               "angle_diff_outer_phi_peak_final", 361, 0, 361, 100, -1.5, 1.5);
    angle_diff_outer_phi_peak_final->SetStats(false);
    angle_diff_outer_phi_peak_final->GetXaxis()->SetTitle("Outer phi [degree]");
    angle_diff_outer_phi_peak_final->GetYaxis()->SetTitle("Inner - Outer [degree]");
    angle_diff_outer_phi_peak_final->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(angle_diff_outer_phi_peak_final);

    DCA_distance_outer_phi_peak_final = new TH2F("DCA_distance_outer_phi_peak_final",
                                                 "DCA_distance_outer_phi_peak_final", 100, 0, 360, 100, -10, 10);
    DCA_distance_outer_phi_peak_final->SetStats(false);
    DCA_distance_outer_phi_peak_final->GetXaxis()->SetTitle("Outer phi [degree]");
    DCA_distance_outer_phi_peak_final->GetYaxis()->SetTitle("DCA [mm]");
    DCA_distance_outer_phi_peak_final->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(DCA_distance_outer_phi_peak_final);

    angle_diff_new_bkg_remove_final = new TH1F("angle_diff_new_bkg_remove_final",
                                               "angle_diff_new_bkg_remove_final", 100, angle_diff_new_l, angle_diff_new_r);
    angle_diff_new_bkg_remove_final->SetStats(false);
    angle_diff_new_bkg_remove_final->GetXaxis()->SetTitle("|Inner - Outer| [degree]");
    angle_diff_new_bkg_remove_final->GetYaxis()->SetTitle("Entry");
    angle_diff_new_bkg_remove_final->GetXaxis()->SetNdivisions(505);
    m_v_hist.push_back(angle_diff_new_bkg_remove_final);
  }
}

// note : this function only prepare the pairs for the vertex XY calculation, it's like a general vertex for the whole run
void INTTXYvtx::ProcessEvt(
    int event_i,
    const std::vector<clu_info>& temp_sPH_inner_nocolumn_vec,
    const std::vector<clu_info>& temp_sPH_outer_nocolumn_vec,
    const std::vector<std::vector<double>>& /*temp_sPH_nocolumn_vec*/,
    const std::vector<std::vector<double>>& /*temp_sPH_nocolumn_rz_vec*/,
    int NvtxMC,
    double /*TrigZvtxMC*/,
    bool PhiCheckTag,
    Long64_t /*bco_full*/
)
{
  if (!m_initialized)
  {
    std::cout << "INTTXYvtx is not initialized, abort in ProcessEvt" << std::endl;
    exit(1);
  }

  if (print_message_opt && event_i % 10000 == 0)
  {
    std::cout << "In INTTXYvtx class, running event : " << event_i << std::endl;
  }

  total_NClus = temp_sPH_inner_nocolumn_vec.size() + temp_sPH_outer_nocolumn_vec.size();

  // note : the Move these two in the beginning of the function, in order to avoid those event-reject cuts
  if (m_enable_qa)
  {
    N_cluster_correlation->Fill(temp_sPH_inner_nocolumn_vec.size(), temp_sPH_outer_nocolumn_vec.size());
    N_cluster_correlation_close->Fill(temp_sPH_inner_nocolumn_vec.size(), temp_sPH_outer_nocolumn_vec.size());
  }

  if (total_NClus < zvtx_cal_require)
  {
    return;
    std::cout << "return confirmation" << std::endl;
  }

  if (run_type == "MC" && NvtxMC != 1)
  {
    return;
    std::cout << "In INTTXYvtx class, event : " << event_i
              << " Nvtx : " << NvtxMC << " Nvtx more than one " << std::endl;
  }
  if (PhiCheckTag == false)
  {
    return;
    std::cout << "In INTTXYvtx class, event : " << event_i
              << " Nvtx : " << NvtxMC << " Not full phi has hits " << std::endl;
  }

  if (temp_sPH_inner_nocolumn_vec.size() < 10 || temp_sPH_outer_nocolumn_vec.size() < 10 || total_NClus > N_clu_cut || total_NClus < N_clu_cutl)
  {
    return;
  }

  //-------------------------------
  // tracklet reconstruction accumulated multiple events
  for (auto& inner_i : temp_sPH_inner_nocolumn_vec)
  {
    for (auto& outer_i : temp_sPH_outer_nocolumn_vec)
    {
      // note : try to ease the analysis and also make it quick.
      if (fabs(inner_i.phi - outer_i.phi) < 7)  // todo : the pre phi cut is here, can be optimized
      {
        cluster_pair_vec.push_back({{inner_i.x,
                                     inner_i.y},
                                    {outer_i.x,
                                     outer_i.y}});
      }
    }
  }

  //-------------------------------
  // QA histogram
  if (m_enable_qa)
  {
    for (auto& inner_i : temp_sPH_inner_nocolumn_vec)
    {
      inner_pos_xy->Fill(inner_i.x, inner_i.y);
      inner_outer_pos_xy->Fill(inner_i.x, inner_i.y);
    }

    for (auto& outer_i : temp_sPH_outer_nocolumn_vec)
    {
      outer_pos_xy->Fill(outer_i.x, outer_i.y);
      inner_outer_pos_xy->Fill(outer_i.x, outer_i.y);
    }
  }

  //--  std::cout<<"  "<<event_i<<" clusterpair:size : "<<cluster_pair_vec.size()<<std::endl;
}

void INTTXYvtx::ClearEvt()
{
  return;
}

unsigned long INTTXYvtx::GetVecNele()
{
  return cluster_pair_vec.size();
}

std::pair<double, double> INTTXYvtx::GetFinalVTXxy()
{
  return {current_vtxX, current_vtxY};
}

std::pair<std::vector<TH2F*>, std::vector<TH1F*>> INTTXYvtx::GetHistFinal()
{
  return {
      {DCA_distance_inner_phi_peak_final,
       angle_diff_inner_phi_peak_final,
       DCA_distance_outer_phi_peak_final,
       angle_diff_outer_phi_peak_final,
       xy_hist,
       xy_hist_bkgrm},
      {angle_diff_new_bkg_remove_final}};
}

void INTTXYvtx::PrintPlots()
{
  if (!m_initialized)
  {
    std::cout << "INTTXYvtx is not initialized, abort in PrintPlots" << std::endl;
    exit(1);
  }

  if (m_enable_drawhist && m_enable_qa)
  {
    std::string s_inttlabel = (boost::format("#it{#bf{sPHENIX INTT}} %s") % plot_text).str().c_str();
    // note : -----------------------------------------------------------------------------------------
    inner_outer_pos_xy->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
    c1->Print((boost::format("%s/xyvtx_qa.pdf(") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : -----------------------------------------------------------------------------------------
    inner_pos_xy->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
    c1->Print((boost::format("%s/xyvtx_qa.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : -----------------------------------------------------------------------------------------
    outer_pos_xy->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
    c1->Print((boost::format("%s/xyvtx_qa.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : -----------------------------------------------------------------------------------------
    N_cluster_correlation->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
    c1->Print((boost::format("%s/xyvtx_qa.pdf") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : -----------------------------------------------------------------------------------------
    N_cluster_correlation_close->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
    c1->Print((boost::format("%s/xyvtx_qa.pdf)") % out_folder_directory).str().c_str());
    c1->Clear();
  }
}

//--------------------------------------------------------------------------
// quadrant method
std::vector<std::pair<double, double>> INTTXYvtx::MacroVTXSquare(double length, int N_trial)
{
  if (!m_initialized)
  {
    std::cout << "INTTXYvtx is not initialized, abort in MacroVTXSquare" << std::endl;
    exit(1);
  }

  bool draw_plot_opt = m_enable_drawhist;

  const double original_length = length;
  std::pair<double, double> origin = {0, 0};
  std::vector<std::pair<double, double>> vtx_vec = Get4vtx(origin, length);  // vtx_vec.push_back(origin);

  int small_index{0};
  std::vector<double> small_info_vec(18, -999);
  std::vector<double> grr_x{};
  std::vector<double> grr_E{};
  std::vector<double> grr_y{};

  std::vector<double> All_FitError_DCA_Y{};
  std::vector<double> All_FitError_DCA_X{};
  std::vector<double> All_FitError_angle_Y{};
  std::vector<double> All_FitError_angle_X{};

  std::vector<double> Winner_FitError_DCA_Y{};
  std::vector<double> Winner_FitError_DCA_X{};
  std::vector<double> Winner_FitError_angle_Y{};
  std::vector<double> Winner_FitError_angle_X{};

  if (print_message_opt == true)
  {
    std::cout << "In INTTXYvtx::MacroVTXSquare, N pairs : " << cluster_pair_vec.size() << std::endl;

    std::cout << N_trial << " runs, smart. which gives you the resolution down to "
              << length / pow(2, N_trial) << " mm" << std::endl;
  }

  if (cluster_pair_vec.size() == 0)
  {  // minimum tracklet cut. need to be tuned
    return {
        beam_origin,  // note : the best vertex
        origin,       // note : the origin in that trial
        {0, 0},       // note : horizontal_fit_inner -> GetParError(0),  horizontal_fit_angle_diff_inner -> GetParError(0)
        {0, 0},       // note : horizontal_fit_inner -> GetParameter(0), horizontal_fit_angle_diff_inner -> GetParameter(0)
        {0, 0},       // note : horizontal_fit_outer -> GetParError(0),  horizontal_fit_angle_diff_outer -> GetParError(0)
        {0, 0},       // note : horizontal_fit_outer -> GetParameter(0), horizontal_fit_angle_diff_outer -> GetParameter(0)
        {0, 0},       // note : the mean and stddev of angle_diff
        {0, 0},       // note : the mean and stddev of DCA_distance
        {0, 0},       // note : the mean and stddev of angle_diff, but with the background removed
    };
  }

  //-------------------------------
  // info_vec contents
  //   0 : angle_diff_new_bkg_remove       -> Mean,
  //   1 : angle_diff_new_bkg_remove       -> StdDev,      // note : angle diff stddev and error (1D histogram)
  //   2 : horizontal_fit_inner            -> Chisq / NDF),
  //   3 : horizontal_fit_inner            -> ParErr(0),   // note : inner DCA, pol0
  //   4 : horizontal_fit_angle_diff_inner -> Chisq / NDF),
  //   5 : horizontal_fit_angle_diff_inner -> ParErr(0),   // note : inner angle diff, pol0
  //   6 : horizontal_fit_outer            -> Chisq / NDF),
  //   7 : horizontal_fit_outer            -> ParErr(0),   // note : outer DCA, pol0
  //   8 : horizontal_fit_angle_diff_outer -> Chisq / NDF),
  //   9 : horizontal_fit_angle_diff_outer -> ParErr(0),   // note : outer angle diff, pol0
  //  10 : horizontal_fit_inner            -> Par(0),
  //  11 : horizontal_fit_angle_diff_inner -> Par(0),
  //  12 : horizontal_fit_outer            -> Par(0),
  //  13 : horizontal_fit_angle_diff_outer -> Par(0),
  //  14 : angle_diff                      -> Mean,
  //  15 : angle_diff                      -> StdDev,
  //  16 : DCA_distance                    -> Mean,
  //  17 : DCA_distance                    -> StdDev
  //

  if (m_enable_drawhist)
  {
    c1->cd();
    c1->Range(0, 0, 1, 1);
    ltx->DrawLatex(0.5, 0.5, "QA plots for quadrant method");
    c1->Print((boost::format("%s/%s(") % out_folder_directory % m_quad_pdfname).str().c_str());
    c1->Clear();
  }

  // current algorithm uses info_vec[3] and [5], others are for QA
  for (int i = 0; i < N_trial; i++)
  {
    if (print_message_opt == true)
    {
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~"
                << " step " << i << " "
                << "~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    }
    for (unsigned int i1 = 0; i1 < vtx_vec.size(); i1++)
    {
      if (print_message_opt == true)
      {
        std::cout << "tested vertex : " << vtx_vec[i1].first << " " << vtx_vec[i1].second << std::endl;
      }

      if (m_enable_drawhist)
      {
        c1->cd();
        c1->Range(0, 0, 1, 1);
        ltx->DrawLatex(0.5, 0.5, (boost::format("New_trial_square_%i_%i") % i % i1).str().c_str());
        c1->Print((boost::format("%s/%s") % out_folder_directory.c_str() % m_quad_pdfname).str().c_str());
        c1->Clear();
      }

      current_vtxX = vtx_vec[i1].first;
      current_vtxY = vtx_vec[i1].second;

      std::vector<double> info_vec = subMacroVTXxyCorrection(i, i1, draw_plot_opt);

      if (print_message_opt == true)
      {
        std::cout << "trial : " << i
                  << " vertex : " << i1
                  << " DCA fit error : " << info_vec[3]
                  << " angle diff fit error : " << info_vec[5] << std::endl;
      }

      All_FitError_DCA_Y.push_back(info_vec[3]);
      All_FitError_DCA_X.push_back(i);
      All_FitError_angle_Y.push_back(info_vec[5]);
      All_FitError_angle_X.push_back(i);

      if ((i1 == 0) ||
          (info_vec[3] < small_info_vec[3] &&
           info_vec[5] < small_info_vec[5])  // note : the fit error of the pol0 fit
      )
      {
        small_info_vec = info_vec;
        small_index = i1;

        TH2F_FakeClone(DCA_distance_inner_phi_peak, DCA_distance_inner_phi_peak_final);
        TH2F_FakeClone(angle_diff_inner_phi_peak, angle_diff_inner_phi_peak_final);
        if (m_enable_qa)
        {
          TH2F_FakeClone(DCA_distance_outer_phi_peak, DCA_distance_outer_phi_peak_final);
          TH2F_FakeClone(angle_diff_outer_phi_peak, angle_diff_outer_phi_peak_final);
          TH1F_FakeClone(angle_diff_new_bkg_remove, angle_diff_new_bkg_remove_final);
        }
      }
      if (print_message_opt == true)
      {
        std::cout << " " << std::endl;
      }

      ClearHist(1);
    }

    if (print_message_opt == true)
    {
      std::cout << "the Quadrant " << small_index << " won the competition" << std::endl;
    }

    Winner_FitError_DCA_Y.push_back(small_info_vec[3]);
    Winner_FitError_DCA_X.push_back(i);
    Winner_FitError_angle_Y.push_back(small_info_vec[5]);
    Winner_FitError_angle_X.push_back(i);

    grr_x.push_back(i);
    grr_y.push_back(small_index);
    grr_E.push_back(0);

    // note : generating the new 4 vertex for the next comparison
    // note : start to shrink the square
    if (i != N_trial - 1)
    {
      origin = {(vtx_vec[small_index].first + origin.first) / 2.,
                (vtx_vec[small_index].second + origin.second) / 2.};

      // std::cout<<"test : "<<origin.first<<" "<<origin.second<<" length: "<<length<<std::endl;
      // if (small_index == 4) {length /= 1.5;}
      // else {length /= 2.;}
      length /= 2.;
      vtx_vec = Get4vtx(origin, length);  // vtx_vec.push_back(origin);
    }
  }

  if (m_enable_drawhist)
  {
    c1->cd();
    c1->Print((boost::format("%s/%s)") % out_folder_directory.c_str() % m_quad_pdfname).str().c_str());
    c1->Clear();
  }

  if (draw_plot_opt == true)
  {
    DrawTGraphErrors(grr_x, grr_y, grr_E, grr_E, out_folder_directory,
                     {
                         (boost::format("Square_scan_history_%.1fmm_%iTrials") % original_length % N_trial).str().c_str()  // title
                         ,
                         "nth scan"  // x_title
                         ,
                         "Winner index"  // y_title
                         ,
                         "APL"  // draw option
                         ,
                         "quadorant_qa.pdf("  // pdf name
                     });
    Draw2TGraph(All_FitError_angle_X, All_FitError_angle_Y, Winner_FitError_angle_X, Winner_FitError_angle_Y, out_folder_directory,
                {
                    (boost::format("Angle_diff_fit_error_%iTrials") % N_trial).str().c_str()  // title
                    ,
                    "n iteration"  // x_title
                    ,
                    "#Delta#phi fit error [degree]"  // y_title
                    ,
                    ""  // draw option    // draw option
                    ,
                    "quadorant_qa.pdf"  // pdf name   // pdf name
                });
    Draw2TGraph(All_FitError_DCA_X, All_FitError_DCA_Y, Winner_FitError_DCA_X, Winner_FitError_DCA_Y, out_folder_directory,
                {
                    (boost::format("DCA_fit_error_%iTrials") % N_trial).str().c_str()  // title
                    ,
                    "n iteration"  // x_title
                    ,
                    "DCA fit error [mm]"  // y_title
                    ,
                    ""  // draw option
                    ,
                    "quadorant_qa.pdf)"  // pdf name
                });
  }

  return {
      vtx_vec[small_index],                      // note : the best vertex
      origin,                                    // note : the origin in that trial
      {small_info_vec[3], small_info_vec[5]},    // note : horizontal_fit_inner -> GetParError(0),  horizontal_fit_angle_diff_inner -> GetParError(0)
      {small_info_vec[10], small_info_vec[11]},  // note : horizontal_fit_inner -> GetParameter(0), horizontal_fit_angle_diff_inner -> GetParameter(0)
      {small_info_vec[7], small_info_vec[9]},    // note : horizontal_fit_outer -> GetParError(0),  horizontal_fit_angle_diff_outer -> GetParError(0)
      {small_info_vec[12], small_info_vec[13]},  // note : horizontal_fit_outer -> GetParameter(0), horizontal_fit_angle_diff_outer -> GetParameter(0)
      {small_info_vec[14], small_info_vec[15]},  // note : the mean and stddev of angle_diff
      {small_info_vec[16], small_info_vec[17]},  // note : the mean and stddev of DCA_distance
      {small_info_vec[0], small_info_vec[1]},    // note : the mean and stddev of angle_diff, but with the background removed
  };
}

std::vector<double> INTTXYvtx::subMacroVTXxyCorrection(int test_index, int trial_index, bool draw_plot_opt)
{
  int true_trial_index = test_index * 4 + trial_index;
  std::vector<double> out_vec = GetVTXxyCorrection_new(true_trial_index);

  std::string sub_out_folder_name{};
  if (draw_plot_opt == true)
  {
    sub_out_folder_name = (boost::format("%s/New_trial_square_%i_%i") % out_folder_directory % test_index % trial_index).str();

    //--if (std::filesystem::exists(sub_out_folder_name.c_str()) == false) {
    //--  system((boost::format("mkdir %s",sub_out_folder_name.c_str()));
    //--}

    // PrintPlotsVTXxy(sub_out_folder_name);
    PrintPlotsVTXxy();
  }

  return out_vec;
}

// note : {circle radius, possible correction angle, the chi2/NDF of pol0 fit}
std::vector<double> INTTXYvtx::GetVTXxyCorrection_new(int trial_index)
{
  if (print_message_opt == true)
  {
    std::cout << "Trial : " << trial_index
              << "---------------------------- ---------------------------- ----------------------------" << std::endl;
    std::cout << "Given vertex: " << current_vtxX << " " << current_vtxY << std::endl;
  }

  if (m_enable_qa)
  {
    cos_fit->SetParameters(4, 1.49143e-02, 185, 0.3);  // todo : we may have to apply more constaints on the fitting
    cos_fit->SetParLimits(2, 0, 360);                  // note : the peak location has to be positive

    // note : here is the test with a gaus fitting to find the peak
    gaus_fit->SetParameters(-4.5, 197, 50, 0);
    gaus_fit->SetParLimits(0, -100, 0);  // note : the gaus distribution points down
    // DCA_distance_inner_phi_peak_profile -> Fit(gaus_fit, "N","",100, 260);
    // std::cout<<"test, gaus fit range : "<<gaus_fit->GetParameter(1) - 25<<" "<<gaus_fit->GetParameter(1) + 25<<std::endl;
  }

  subMacroPlotWorking(true, 100, 260, 25);

  return {
      angle_diff_new_bkg_remove->GetMean(),
      angle_diff_new_bkg_remove->GetStdDev(),  // note : angle diff stddev and error (1D histogram)
      horizontal_fit_inner->GetChisquare() / double(horizontal_fit_inner->GetNDF()),
      horizontal_fit_inner->GetParError(0),  // note : inner DCA, pol0
      horizontal_fit_angle_diff_inner->GetChisquare() / double(horizontal_fit_angle_diff_inner->GetNDF()),
      horizontal_fit_angle_diff_inner->GetParError(0),  // note : inner angle diff, pol0
      horizontal_fit_outer->GetChisquare() / double(horizontal_fit_outer->GetNDF()),
      horizontal_fit_outer->GetParError(0),  // note : outer DCA, pol0
      horizontal_fit_angle_diff_outer->GetChisquare() / double(horizontal_fit_angle_diff_outer->GetNDF()),
      horizontal_fit_angle_diff_outer->GetParError(0),  // note : outer angle diff, pol0
      horizontal_fit_inner->GetParameter(0),
      horizontal_fit_angle_diff_inner->GetParameter(0),  // note : 10, 11
      horizontal_fit_outer->GetParameter(0),
      horizontal_fit_angle_diff_outer->GetParameter(0),  // note : 12, 13
      angle_diff->GetMean(),
      angle_diff->GetStdDev(),  // note : 14, 15
      DCA_distance->GetMean(),
      DCA_distance->GetStdDev(),  // note : 16, 17
  };
}

void INTTXYvtx::subMacroPlotWorking(
    bool phi_correction,
    double cos_fit_rangel,
    double cos_fit_ranger,
    double guas_fit_range)
{
  //   3 : horizontal_fit_inner            -> ParErr(0),   // note : inner DCA, pol0
  //   5 : horizontal_fit_angle_diff_inner -> ParErr(0),   // note : inner angle diff, pol0

  for (auto& i : cluster_pair_vec)
  {
    std::vector<double> DCA_info_vec = calculateDistanceAndClosestPoint(
        i.first.x, i.first.y,
        i.second.x, i.second.y,
        current_vtxX, current_vtxY);

    double DCA_sign = calculateAngleBetweenVectors(
        i.second.x, i.second.y,
        i.first.x, i.first.y,
        current_vtxX, current_vtxY);

    if (phi_correction == true)
    {
      // std::cout<<"option selected "<<std::endl;
      Clus_InnerPhi_Offset = (i.first.y - current_vtxY < 0)
                                 ? atan2(i.first.y - current_vtxY, i.first.x - current_vtxX) * (180. / M_PI) + 360
                                 : atan2(i.first.y - current_vtxY, i.first.x - current_vtxX) * (180. / M_PI);
      Clus_OuterPhi_Offset = (i.second.y - current_vtxY < 0)
                                 ? atan2(i.second.y - current_vtxY, i.second.x - current_vtxX) * (180. / M_PI) + 360
                                 : atan2(i.second.y - current_vtxY, i.second.x - current_vtxX) * (180. / M_PI);
    }
    else  // note : phi_correction == false
    {
      Clus_InnerPhi_Offset = (i.first.y < 0)
                                 ? atan2(i.first.y, i.first.x) * (180. / M_PI) + 360
                                 : atan2(i.first.y, i.first.x) * (180. / M_PI);
      Clus_OuterPhi_Offset = (i.second.y < 0)
                                 ? atan2(i.second.y, i.second.x) * (180. / M_PI) + 360
                                 : atan2(i.second.y, i.second.x) * (180. / M_PI);
    }

    //----------------------
    // this is used for quadrant method
    DCA_distance_inner_phi->Fill(Clus_InnerPhi_Offset, DCA_sign);
    angle_diff_inner_phi->Fill(Clus_InnerPhi_Offset, Clus_InnerPhi_Offset - Clus_OuterPhi_Offset);

    //----------------------
    if (m_enable_qa)
    {
      DCA_distance_outer_phi->Fill(Clus_OuterPhi_Offset, DCA_sign);
      angle_diff_outer_phi->Fill(Clus_OuterPhi_Offset, Clus_InnerPhi_Offset - Clus_OuterPhi_Offset);

      angle_diff->Fill(std::abs(Clus_InnerPhi_Offset - Clus_OuterPhi_Offset));
      angle_diff_new->Fill(std::abs(Clus_InnerPhi_Offset - Clus_OuterPhi_Offset));
      DCA_distance->Fill(DCA_sign);

      // draw only
      angle_correlation->Fill(Clus_InnerPhi_Offset, Clus_OuterPhi_Offset);
      angle_diff_DCA_dist->Fill(Clus_InnerPhi_Offset - Clus_OuterPhi_Offset, DCA_sign);
      DCA_point->Fill(DCA_info_vec[1], DCA_info_vec[2]);

      DCA_distance_inner_X->Fill(i.first.x, DCA_sign);
      DCA_distance_inner_Y->Fill(i.first.y, DCA_sign);
      DCA_distance_outer_X->Fill(i.second.x, DCA_sign);
      DCA_distance_outer_Y->Fill(i.second.y, DCA_sign);
    }

  }  // note : end of the loop for the cluster pair

  // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
  delete DCA_distance_inner_phi_peak;
  DCA_distance_inner_phi_peak = (TH2F*) DCA_distance_inner_phi->Clone("DCA_distance_inner_phi_peak");
  TH2F_threshold(DCA_distance_inner_phi_peak, 0.5);  // todo : the background cut can be modified, the ratio 0.5
  delete DCA_distance_inner_phi_peak_profile;
  DCA_distance_inner_phi_peak_profile = DCA_distance_inner_phi_peak->ProfileX("DCA_distance_inner_phi_peak_profile");

  double point_index = 0;
  std::vector<double> hist_column_content = SumTH2FColumnContent(DCA_distance_inner_phi_peak);
  for (int i = 0; i < DCA_distance_inner_phi_peak_profile->GetNbinsX(); i++)
  {
    if (hist_column_content[i] < 5)
    {
      continue;
    }  // note : in order to remove some remaining background

    DCA_distance_inner_phi_peak_profile_graph->SetPoint(point_index,
                                                        DCA_distance_inner_phi_peak_profile->GetBinCenter(i + 1),
                                                        DCA_distance_inner_phi_peak_profile->GetBinContent(i + 1));
    // std::cout<<"("<<DCA_distance_inner_phi_peak_profile->GetBinCenter(i+1)<<", "<< DCA_distance_inner_phi_peak_profile->GetBinContent(i+1)<<")"<<std::endl;
    point_index += 1;
  }

  //------------------------------------------------------------------
  // this is used to constrain the quadrant
  // info_vec[3];
  // todo : the fit range of the gaussian fit can be modified here
  //
  horizontal_fit_inner->SetParameter(0, 0);
  DCA_distance_inner_phi_peak_profile_graph->Fit(horizontal_fit_inner, "NQ", "", 0, 360);
  //------------------------------------------------------------------

  if (m_enable_qa)
  {
    DCA_distance_inner_phi_peak_profile_graph->Fit(gaus_fit, "NQ", "",
                                                   cos_fit->GetParameter(2) - guas_fit_range,   // note : what we want and need is the peak position,
                                                   cos_fit->GetParameter(2) + guas_fit_range);  // so we fit the peak again
    DCA_distance_inner_phi_peak_profile_graph->Fit(cos_fit, "NQ", "", cos_fit_rangel, cos_fit_ranger);
  }

  // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
  delete angle_diff_inner_phi_peak;
  angle_diff_inner_phi_peak = (TH2F*) angle_diff_inner_phi->Clone("angle_diff_inner_phi_peak");
  TH2F_threshold_advanced_2(angle_diff_inner_phi_peak, 0.5);  // todo : threshold ratio can be modified here
  hist_column_content = SumTH2FColumnContent(angle_diff_inner_phi_peak);
  angle_diff_inner_phi_peak_profile = angle_diff_inner_phi_peak->ProfileX("angle_diff_inner_phi_peak_profile");
  point_index = 0;
  for (int i = 0; i < angle_diff_inner_phi_peak_profile->GetNbinsX(); i++)
  {
    if (hist_column_content[i] < 5)
    {
      continue;
    }  // note : in order to remove some remaining background

    angle_diff_inner_phi_peak_profile_graph->SetPoint(point_index,
                                                      angle_diff_inner_phi_peak_profile->GetBinCenter(i + 1),
                                                      angle_diff_inner_phi_peak_profile->GetBinContent(i + 1));
    // std::cout<<"("<<angle_diff_inner_phi_peak_profile->GetBinCenter(i+1)<<", "<< angle_diff_inner_phi_peak_profile->GetBinContent(i+1)<<")"<<std::endl;
    point_index += 1;
  }

  //------------------------------------------------------------------
  // this is used to constrain the quadrant
  // info_vec[5];
  horizontal_fit_angle_diff_inner->SetParameter(0, 0);
  angle_diff_inner_phi_peak_profile_graph->Fit(horizontal_fit_angle_diff_inner, "NQ", "", 0, 360);
  //------------------------------------------------------------------

  angle_diff_inner_phi_peak_profile_graph->Set(0);
  DCA_distance_outer_phi_peak_profile_graph->Set(0);

  //------------------------------------------------------------------------
  // all others are not used for XY vertex calculation. for QA
  //------------------------------------------------------------------------

  // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
  delete angle_diff_new_bkg_remove;
  angle_diff_new_bkg_remove = (TH1F*) angle_diff_new->Clone("angle_diff_new_bkg_remove");
  angle_diff_new_bkg_remove->SetLineColor(2);
  Hist_1D_bkg_remove(angle_diff_new_bkg_remove, 1.5);

  if (m_enable_qa)
  {
    // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
    delete DCA_distance_outer_phi_peak;
    DCA_distance_outer_phi_peak = (TH2F*) DCA_distance_outer_phi->Clone("DCA_distance_outer_phi_peak");
    TH2F_threshold(DCA_distance_outer_phi_peak, 0.5);  // todo : the background cut can be modified, the ratio 0.5
    DCA_distance_outer_phi_peak_profile = DCA_distance_outer_phi_peak->ProfileX("DCA_distance_outer_phi_peak_profile");
    point_index = 0;
    hist_column_content = SumTH2FColumnContent(DCA_distance_outer_phi_peak);
    for (int i = 0; i < DCA_distance_outer_phi_peak_profile->GetNbinsX(); i++)
    {
      if (hist_column_content[i] < 5)
      {
        continue;
      }  // note : in order to remove some remaining background

      DCA_distance_outer_phi_peak_profile_graph->SetPoint(point_index,
                                                          DCA_distance_outer_phi_peak_profile->GetBinCenter(i + 1),
                                                          DCA_distance_outer_phi_peak_profile->GetBinContent(i + 1));
      // std::cout<<"("<<DCA_distance_outer_phi_peak_profile->GetBinCenter(i+1)<<", "<< DCA_distance_outer_phi_peak_profile->GetBinContent(i+1)<<")"<<std::endl;
      point_index += 1;
    }

    horizontal_fit_outer->SetParameter(0, 0);
    // todo : the fit range of the gaussian fit can be modified here
    DCA_distance_outer_phi_peak_profile_graph->Fit(horizontal_fit_outer, "NQ", "", 0, 360);

    // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
    delete angle_diff_outer_phi_peak;
    angle_diff_outer_phi_peak = (TH2F*) angle_diff_outer_phi->Clone("angle_diff_outer_phi_peak");
    TH2F_threshold_advanced_2(angle_diff_outer_phi_peak, 0.5);  // todo : threshold ratio can be modified here
    hist_column_content = SumTH2FColumnContent(angle_diff_outer_phi_peak);
    angle_diff_outer_phi_peak_profile = angle_diff_outer_phi_peak->ProfileX("angle_diff_outer_phi_peak_profile");
    point_index = 0;
    for (int i = 0; i < angle_diff_outer_phi_peak_profile->GetNbinsX(); i++)
    {
      if (hist_column_content[i] < 5)
      {
        continue;
      }  // note : in order to remove some remaining background

      angle_diff_outer_phi_peak_profile_graph->SetPoint(point_index,
                                                        angle_diff_outer_phi_peak_profile->GetBinCenter(i + 1),
                                                        angle_diff_outer_phi_peak_profile->GetBinContent(i + 1));
      // std::cout<<"("<<angle_diff_outer_phi_peak_profile->GetBinCenter(i+1)<<", "<< angle_diff_outer_phi_peak_profile->GetBinContent(i+1)<<")"<<std::endl;
      point_index += 1;
    }

    horizontal_fit_angle_diff_outer->SetParameter(0, 0);
    angle_diff_outer_phi_peak_profile_graph->Fit(horizontal_fit_angle_diff_outer, "NQ", "", 0, 360);

    // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------

    angle_diff_outer_phi_peak_profile_graph->Set(0);
    DCA_distance_inner_phi_peak_profile_graph->Set(0);
  }

  if (m_enable_qa && print_message_opt == true)
  {
    std::cout << "circle radius : " << std::abs(gaus_fit->GetParameter(0) + gaus_fit->GetParameter(3))
              << " possible correction angle : " << gaus_fit->GetParameter(1) << std::endl;
  }
}

// void INTTXYvtx::PrintPlotsVTXxy(std::string sub_out_folder_name)
void INTTXYvtx::PrintPlotsVTXxy()
{
  if (m_enable_drawhist)
  {
    std::string s_inttlabel = (boost::format("#it{#bf{sPHENIX INTT}} %s") % plot_text).str().c_str();
    std::string s_pdfname = out_folder_directory + "/" + m_quad_pdfname;
    std::cout << s_pdfname << std::endl;

    // note : -----------------------------------------------------------------------------------------
    DCA_distance_inner_phi->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
    c1->Print(s_pdfname.c_str());
    c1->Clear();

    // note : -----------------------------------------------------------------------------------------
    angle_diff_inner_phi->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
    c1->Print(s_pdfname.c_str());
    c1->Clear();

    // note : -----------------------------------------------------------------------------------------
    DCA_distance_inner_phi_peak->SetStats(false);
    DCA_distance_inner_phi_peak->GetXaxis()->SetTitle("Inner phi [degree]");
    DCA_distance_inner_phi_peak->GetYaxis()->SetTitle("DCA [mm]");
    DCA_distance_inner_phi_peak->Draw("colz0");
    DCA_distance_inner_phi_peak_profile->Draw("same");
    // cos_fit -> Draw("l same");
    // gaus_fit -> Draw("l same");
    horizontal_fit_inner->Draw("l same");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
    draw_text->DrawLatex(0.25, 0.84, (boost::format("#color[2]{Assumed vertex: %.3f mm, %.3f mm}") % current_vtxX % current_vtxY).str().c_str());
    draw_text->DrawLatex(0.25, 0.80, (boost::format("#color[2]{Pol0 fit chi2/NDF: %.3f, fit error: %.3f}") % (horizontal_fit_inner->GetChisquare() / double(horizontal_fit_inner->GetNDF())) % horizontal_fit_inner->GetParError(0)).str().c_str());
    c1->Print(s_pdfname.c_str());
    c1->Clear();

    // note : -----------------------------------------------------------------------------------------
    angle_diff_inner_phi_peak->SetStats(false);
    angle_diff_inner_phi_peak->GetXaxis()->SetTitle("Inner phi [degree]");
    angle_diff_inner_phi_peak->GetYaxis()->SetTitle("Inner - Outer [degree]");
    angle_diff_inner_phi_peak->Draw("colz0");
    angle_diff_inner_phi_peak_profile->Draw("same");
    horizontal_fit_angle_diff_inner->Draw("l same");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
    draw_text->DrawLatex(0.25, 0.84, (boost::format("#color[2]{Assumed vertex: %.3f mm, %.3f mm}") % current_vtxX % current_vtxY).str().c_str());
    c1->Print(s_pdfname.c_str());
    c1->Clear();

    //----------------------
    if (m_enable_qa)
    {
      // note : -----------------------------------------------------------------------------------------
      DCA_distance_outer_phi->Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      c1->Print(s_pdfname.c_str());
      c1->Clear();

      // note : -----------------------------------------------------------------------------------------
      angle_diff_outer_phi->Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      c1->Print(s_pdfname.c_str());
      c1->Clear();

      // note : -----------------------------------------------------------------------------------------
      angle_diff->SetMinimum(0);
      angle_diff->Draw("hist");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      c1->Print(s_pdfname.c_str());
      c1->Clear();

      // note : -----------------------------------------------------------------------------------------
      angle_diff_new->SetMinimum(0);
      angle_diff_new->Draw("hist");
      angle_diff_new_bkg_remove->Draw("hist same");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      draw_text->DrawLatex(0.4, 0.80, (boost::format("#color[2]{Dist. StdDev: %.4f}") % angle_diff_new_bkg_remove->GetStdDev()).str().c_str());
      c1->Print(s_pdfname.c_str());
      c1->Clear();

      // note : -----------------------------------------------------------------------------------------
      DCA_distance->Draw("hist");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      c1->Print(s_pdfname.c_str());
      c1->Clear();

      // note : -----------------------------------------------------------------------------------------
      angle_correlation->Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      c1->Print(s_pdfname.c_str());
      c1->Clear();

      // note : -----------------------------------------------------------------------------------------
      angle_diff_DCA_dist->Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      c1->Print(s_pdfname.c_str());
      c1->Clear();

      // note : -----------------------------------------------------------------------------------------
      DCA_point->Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      c1->Print(s_pdfname.c_str());
      c1->Clear();

      // note : -----------------------------------------------------------------------------------------
      DCA_distance_inner_X->Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      c1->Print(s_pdfname.c_str());
      c1->Clear();

      // note : -----------------------------------------------------------------------------------------
      DCA_distance_inner_Y->Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      c1->Print(s_pdfname.c_str());
      c1->Clear();

      // note : -----------------------------------------------------------------------------------------
      DCA_distance_outer_X->Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      c1->Print(s_pdfname.c_str());
      c1->Clear();

      // note : -----------------------------------------------------------------------------------------
      DCA_distance_outer_Y->Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      c1->Print(s_pdfname.c_str());
      c1->Clear();

      // note : -----------------------------------------------------------------------------------------
      DCA_distance_outer_phi_peak->SetStats(false);
      DCA_distance_outer_phi_peak->GetXaxis()->SetTitle("Outer phi [degree]");
      DCA_distance_outer_phi_peak->GetYaxis()->SetTitle("DCA [mm]");
      DCA_distance_outer_phi_peak->Draw("colz0");
      DCA_distance_outer_phi_peak_profile->Draw("same");
      horizontal_fit_outer->Draw("l same");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      draw_text->DrawLatex(0.25, 0.80, (boost::format("#color[2]{Assumed vertex: %.3f mm, %.3f mm}") % current_vtxX % current_vtxY).str().c_str());
      c1->Print(s_pdfname.c_str());
      c1->Clear();

      // note : -----------------------------------------------------------------------------------------
      angle_diff_outer_phi_peak->SetStats(false);
      angle_diff_outer_phi_peak->GetXaxis()->SetTitle("Outer phi [degree]");
      angle_diff_outer_phi_peak->GetYaxis()->SetTitle("Inner - Outer [degree]");
      angle_diff_outer_phi_peak->Draw("colz0");
      angle_diff_outer_phi_peak_profile->Draw("same");
      horizontal_fit_angle_diff_outer->Draw("l same");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, (boost::format("%s, peak : %f") % s_inttlabel.c_str() % peek).str().c_str());
      draw_text->DrawLatex(0.25, 0.84, (boost::format("#color[2]{Assumed vertex: %.3f mm, %.3f mm}") % current_vtxX % current_vtxY).str().c_str());
      c1->Print(s_pdfname.c_str());
      c1->Clear();
    }
  }
}

void INTTXYvtx::ClearHist(int /*print_option*/)
{
  // clear histograms for quadrant method
  DCA_distance_inner_phi->Reset("ICESM");
  angle_diff_inner_phi->Reset("ICESM");

  if (m_enable_qa)
  {
    DCA_distance_outer_phi->Reset("ICESM");
    angle_diff_outer_phi->Reset("ICESM");

    angle_correlation->Reset("ICESM");
    angle_diff_DCA_dist->Reset("ICESM");
    angle_diff->Reset("ICESM");
    DCA_point->Reset("ICESM");
    DCA_distance->Reset("ICESM");
    DCA_distance_inner_X->Reset("ICESM");
    DCA_distance_inner_Y->Reset("ICESM");
    DCA_distance_outer_X->Reset("ICESM");
    DCA_distance_outer_Y->Reset("ICESM");

    angle_diff_new->Reset("ICESM");
    //
    // DCA_distance_inner_phi_peak -> Reset("ICESM");
    // DCA_distance_outer_phi_peak -> Reset("ICESM");
    // angle_diff_inner_phi_peak -> Reset("ICESM");
    // angle_diff_outer_phi_peak -> Reset("ICESM");
    // angle_diff_new_bkg_remove -> Reset("ICESM");
  }

  delete angle_diff_new_bkg_remove;
  angle_diff_new_bkg_remove = nullptr;
  delete angle_diff_inner_phi_peak;
  angle_diff_inner_phi_peak = nullptr;
  delete angle_diff_outer_phi_peak;
  angle_diff_outer_phi_peak = nullptr;
  delete DCA_distance_inner_phi_peak;
  DCA_distance_inner_phi_peak = nullptr;
  delete DCA_distance_outer_phi_peak;
  DCA_distance_outer_phi_peak = nullptr;
}

void INTTXYvtx::EndRun()
{
  if (!m_initialized)
  {
    std::cout << "INTTXYvtx is not initialized, abort in EndRun" << std::endl;
    exit(1);
  }

  //------------------------
  // write histograms
  if (m_savehist)
  {
    if (!std::filesystem::exists(out_folder_directory))
    {
      std::filesystem::create_directory(out_folder_directory);
    }

    TFile* file_out = new TFile((boost::format("%s/run_XY_histo.root") % out_folder_directory).str().c_str(), "RECREATE");

    for (auto& itr : m_v_hist)
    {
      itr->Write();
    }

    if (xy_hist != nullptr)
    {
      xy_hist->Write();
    }
    if (xy_hist_bkgrm != nullptr)
    {
      xy_hist_bkgrm->Write();
    }

    file_out->Close();
    delete file_out;
  }

  // ProcessEvt QA histos
  //--    inner_pos_xy -> Reset("ICESM");
  //--    outer_pos_xy -> Reset("ICESM");
  //--    inner_outer_pos_xy -> Reset("ICESM");
  //--    N_cluster_correlation -> Reset("ICESM");
  //--    N_cluster_correlation_close -> Reset("ICESM");

  // file_out -> cd();
  // tree_out -> SetDirectory(file_out);
  // tree_out -> Write();
  // file_out -> Close();

  // c1 -> Delete();
  // ltx -> Delete();
  // draw_text -> Delete();
  // cos_fit -> Delete();
  // gaus_fit -> Delete();

  // horizontal_fit_inner -> Delete();
  // horizontal_fit_angle_diff_inner -> Delete();
  // horizontal_fit_outer -> Delete();
  // horizontal_fit_angle_diff_outer -> Delete();

  //--delete gROOT->FindObject("angle_diff");
  //--delete gROOT->FindObject("angle_diff_new");

  // angle_correlation -> Delete();
  // angle_diff_DCA_dist -> Delete();
  // angle_diff -> Delete();
  // angle_diff_new -> Delete();
  // angle_diff_new_bkg_remove -> Delete();
  // inner_pos_xy -> Delete();
  // outer_pos_xy -> Delete();
  // inner_outer_pos_xy -> Delete();
  // DCA_point -> Delete();
  // DCA_distance -> Delete();
  // N_cluster_correlation -> Delete();
  // N_cluster_correlation_close -> Delete();

  // DCA_distance_inner_X -> Delete();
  // DCA_distance_inner_Y -> Delete();
  // DCA_distance_outer_X -> Delete();
  // DCA_distance_outer_Y -> Delete();

  // DCA_distance_inner_phi -> Delete();
  // DCA_distance_inner_phi_peak -> Delete();
  // DCA_distance_inner_phi_peak_profile -> Delete();
  // DCA_distance_outer_phi -> Delete();
  // DCA_distance_outer_phi_peak -> Delete();
  // DCA_distance_outer_phi_peak_profile -> Delete();

  // angle_diff_inner_phi -> Delete();
  // angle_diff_inner_phi_peak -> Delete();
  // angle_diff_inner_phi_peak_profile -> Delete();
  // angle_diff_outer_phi -> Delete();
  // angle_diff_outer_phi_peak -> Delete();
  // angle_diff_outer_phi_peak_profile -> Delete();

  // DCA_distance_inner_phi_peak_final -> Delete();
  // angle_diff_inner_phi_peak_final -> Delete();

  return;
}

void INTTXYvtx::TH2F_threshold(TH2F* hist, double threshold)
{
  double max_cut = hist->GetMaximum() * threshold;

  for (int xi = 0; xi < hist->GetNbinsX(); xi++)
  {
    for (int yi = 0; yi < hist->GetNbinsY(); yi++)
    {
      if (hist->GetBinContent(xi + 1, yi + 1) < max_cut)
      {
        hist->SetBinContent(xi + 1, yi + 1, 0);
      }
    }
  }
}

void INTTXYvtx::TH2F_threshold_advanced_2(TH2F* hist, double threshold)
{
  // note : this function is to remove the background of the 2D histogram
  // note : but the threshold is given by average of the contents of the top "chosen_bin" bins and timing the threshold
  double max_cut = 0;
  int chosen_bin = 7;

  std::vector<float> all_bin_content_vec{};  //--all_bin_content_vec.clear();
  for (int xi = 0; xi < hist->GetNbinsX(); xi++)
  {
    for (int yi = 0; yi < hist->GetNbinsY(); yi++)
    {
      all_bin_content_vec.push_back(hist->GetBinContent(xi + 1, yi + 1));
    }
  }
  std::vector<unsigned long> ind(all_bin_content_vec.size(), 0);
  TMath::Sort(all_bin_content_vec.size(), &all_bin_content_vec[0], &ind[0]);
  for (int i = 0; i < chosen_bin; i++)
  {
    max_cut += all_bin_content_vec[ind[i]];
    /*std::cout<<"test : "<<all_bin_content_vec[ind[i]]<<std::endl;*/
  }

  max_cut = (max_cut / double(chosen_bin)) * threshold;

  for (int xi = 0; xi < hist->GetNbinsX(); xi++)
  {
    for (int yi = 0; yi < hist->GetNbinsY(); yi++)
    {
      if (hist->GetBinContent(xi + 1, yi + 1) < max_cut)
      {
        hist->SetBinContent(xi + 1, yi + 1, 0);
      }
    }
  }
}

std::vector<double>
INTTXYvtx::calculateDistanceAndClosestPoint(
    double x1,
    double y1,
    double x2,
    double y2,
    double target_x,
    double target_y)
{
  if (x1 != x2)
  {
    // Calculate the slope and intercept of the line passing through (x1, y1) and (x2, y2)
    double a = (y2 - y1) / (x2 - x1);
    double b = y1 - a * x1;

    // std::cout<<"slope : y="<<a<<"x+"<<b<<std::endl;

    // Calculate the closest distance from (target_x, target_y) to the line y = ax + b
    double closest_distance = std::abs(a * target_x - target_y + b) / std::sqrt(a * a + 1);

    // Calculate the coordinates of the closest point (Xc, Yc) on the line y = ax + b
    double Xc = (target_x + a * target_y - a * b) / (a * a + 1);
    double Yc = a * Xc + b;

    return {closest_distance, Xc, Yc};
  }
  else
  {
    double closest_distance = std::abs(x1 - target_x);
    double Xc = x1;
    double Yc = target_y;

    return {closest_distance, Xc, Yc};
  }
}

// note : Function to calculate the angle between two vectors in degrees using the cross product
double INTTXYvtx::calculateAngleBetweenVectors(
    double x1,
    double y1,
    double x2,
    double y2,
    double targetX,
    double targetY)
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

  double DCA_value = sin(angleInRadians_new) * magnitude2;  // DCA_value insread of DCA_distance

  return DCA_value;
}

void INTTXYvtx::Hist_1D_bkg_remove(TH1F* hist_in, double factor)
{
  // todo : N bins considered to be used in the background quantification
  std::vector<double> Nbin_content_vec{};
  for (int i = hist_in->GetNbinsX() - 5; i < hist_in->GetNbinsX(); i++)
  {
    Nbin_content_vec.push_back(hist_in->GetBinContent(i + 1));
  }

  double bkg_level = accumulate(Nbin_content_vec.begin(), Nbin_content_vec.end(), 0.0) / Nbin_content_vec.size();
  // std::cout<<"test, bkg cut : "<<bkg_level * factor<<std::endl;

  for (int i = 0; i < hist_in->GetNbinsX(); i++)
  {
    // note : the background rejection is here : bkg_level * 1.5 for the time being
    double bin_content = (hist_in->GetBinContent(i + 1) <= bkg_level * factor)
                             ? 0.
                             : (hist_in->GetBinContent(i + 1) - bkg_level * factor);

    hist_in->SetBinContent(i + 1, bin_content);
  }
}

void INTTXYvtx::DrawTGraphErrors(
    std::vector<double> x_vec,
    std::vector<double> y_vec,
    std::vector<double> xE_vec,
    std::vector<double> yE_vec,
    const std::string& output_directory,
    std::vector<std::string> plot_name)
{
  if (m_enable_drawhist)
  {
    c1->cd();

    TGraphErrors* g = new TGraphErrors(x_vec.size(), &x_vec[0], &y_vec[0], &xE_vec[0], &yE_vec[0]);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.5);
    g->SetMarkerColor(1);
    g->SetLineColor(1);
    g->GetXaxis()->SetTitle(plot_name[1].c_str());
    g->GetYaxis()->SetTitle(plot_name[2].c_str());
    g->SetTitle(plot_name[0].c_str());
    if (plot_name.size() == 4)
    {
      g->Draw(plot_name[3].c_str());
    }
    else
    {
      g->Draw("AP");
    }

    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, (boost::format("#it{#bf{sPHENIX INTT}} %s") % plot_text).str().c_str());
    c1->Print((boost::format("%s/%s") % output_directory.c_str() % plot_name[4]).str().c_str());
    c1->Clear();

    delete g;
  }
}

void INTTXYvtx::Draw2TGraph(
    std::vector<double> x1_vec,
    std::vector<double> y1_vec,
    std::vector<double> x2_vec,
    std::vector<double> y2_vec,
    const std::string& output_directory,
    std::vector<std::string> plot_name)
{
  if (m_enable_drawhist)
  {
    c1->cd();
    c1->SetLogy(1);

    TGraph* g1 = new TGraph(x1_vec.size(), &x1_vec[0], &y1_vec[0]);
    g1->SetMarkerStyle(5);
    g1->SetMarkerSize(1);
    g1->SetMarkerColor(1);
    g1->SetLineColor(1);
    g1->GetXaxis()->SetTitle(plot_name[1].c_str());
    g1->GetYaxis()->SetTitle(plot_name[2].c_str());
    g1->GetXaxis()->SetNdivisions(505);
    g1->GetXaxis()->SetLimits(-1, x1_vec[x1_vec.size() - 1] + 1);
    g1->SetTitle(plot_name[0].c_str());
    g1->Draw("AP");

    TGraph* g2 = new TGraph(x2_vec.size(), &x2_vec[0], &y2_vec[0]);
    g2->SetMarkerStyle(5);
    g2->SetMarkerSize(1);
    g2->SetMarkerColor(2);
    g2->SetLineColor(2);
    g2->Draw("PL same");

    TLegend* legend = new TLegend(0.4, 0.75, 0.65, 0.9);
    // legend -> SetMargin(0);
    legend->SetTextSize(0.03);
    legend->AddEntry(g1, "Tested vertex candidates", "p");
    legend->AddEntry(g2, "Better performed candidates", "p");
    legend->Draw("same");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, (boost::format("#it{#bf{sPHENIX INTT}} %s") % plot_text).str().c_str());
    c1->Print((boost::format("%s/%s") % output_directory.c_str() % plot_name[4]).str().c_str());
    c1->Clear();
    c1->SetLogy(0);

    delete g1;
    delete g2;
    delete legend;
  }
}

std::vector<double> INTTXYvtx::SumTH2FColumnContent(TH2F* hist_in)
{
  std::vector<double> sum_vec;
  sum_vec.clear();
  for (int i = 0; i < hist_in->GetNbinsX(); i++)
  {
    double sum = 0;
    for (int j = 0; j < hist_in->GetNbinsY(); j++)
    {
      sum += hist_in->GetBinContent(i + 1, j + 1);
    }
    sum_vec.push_back(sum);
  }
  return sum_vec;
}

std::vector<std::pair<double, double>> INTTXYvtx::Get4vtx(std::pair<double, double> origin, double length)
{
  std::vector<std::pair<double, double>> unit_vtx = {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
  std::vector<std::pair<double, double>> vec_out{};  //-- vec_out.clear();

  vec_out.reserve(unit_vtx.size());
  for (std::pair i1 : unit_vtx)
  {
    vec_out.emplace_back(i1.first * length + origin.first, i1.second * length + origin.second);
  }

  return vec_out;
}

void INTTXYvtx::TH2F_FakeClone(TH2F* hist_in, TH2F* hist_out)
{
  if (hist_in->GetNbinsX() != hist_out->GetNbinsX() ||
      hist_in->GetNbinsY() != hist_out->GetNbinsY())
  {
    std::cout << "In INTTXYvtx::TH2F_FakeClone, the input and output histogram have different binning!" << std::endl;
    return;
  }

  for (int i = 0; i < hist_in->GetNbinsX(); i++)
  {
    for (int j = 0; j < hist_in->GetNbinsY(); j++)
    {
      hist_out->SetBinContent(i + 1, j + 1, hist_in->GetBinContent(i + 1, j + 1));
    }
  }
}

void INTTXYvtx::TH1F_FakeClone(TH1F* hist_in, TH1F* hist_out)
{
  if (hist_in->GetNbinsX() != hist_out->GetNbinsX())
  {
    std::cout << "In INTTXYvtx::TH1F_FakeClone, the input and output histogram have different binning!" << std::endl;
    return;
  }

  for (int i = 0; i < hist_in->GetNbinsX(); i++)
  {
    hist_out->SetBinContent(i + 1, hist_in->GetBinContent(i + 1));
  }
}

void INTTXYvtx::TH2FSampleLineFill(
    TH2F* hist_in,
    double segmentation,
    std::pair<double, double> inner_clu,
    std::pair<double, double> outer_clu)
{
  if (!m_initialized)
  {
    std::cout << "INTTXYvtx is not initialized, abort in MacroVTXSquare" << std::endl;
    exit(1);
  }

  double x_min = hist_in->GetXaxis()->GetXmin();
  double x_max = hist_in->GetXaxis()->GetXmax();
  double y_min = hist_in->GetYaxis()->GetXmin();
  double y_max = hist_in->GetYaxis()->GetXmax();

  double seg_x, seg_y;
  double angle;
  int n_seg = 0;

  while (true)
  {
    angle = atan2(inner_clu.second - outer_clu.second, inner_clu.first - outer_clu.first);
    seg_x = (n_seg * segmentation) * cos(angle) + outer_clu.first;  // note : atan2(y,x), point.first is the radius
    seg_y = (n_seg * segmentation) * sin(angle) + outer_clu.second;

    if ((seg_x > x_min && seg_x < x_max && seg_y > y_min && seg_y < y_max) != true)
    {
      break;
    }

    hist_in->Fill(seg_x, seg_y);
    n_seg += 1;
  }

  n_seg = 1;
  while (true)
  {
    angle = atan2(inner_clu.second - outer_clu.second, inner_clu.first - outer_clu.first);
    seg_x = (-1 * n_seg * segmentation) * cos(angle) + outer_clu.first;  // note : atan2(y,x), point.first is the radius
    seg_y = (-1 * n_seg * segmentation) * sin(angle) + outer_clu.second;

    if ((seg_x > x_min && seg_x < x_max && seg_y > y_min && seg_y < y_max) != true)
    {
      break;
    }
    hist_in->Fill(seg_x, seg_y);
    n_seg += 1;
  }
}

std::vector<std::pair<double, double>>
INTTXYvtx::FillLine_FindVertex(
    std::pair<double, double> window_center,
    double segmentation,
    double window_width,
    int N_bins)
{
  bool draw_plot = m_enable_drawhist;

  if (cluster_pair_vec.size() == 0)
  {  // minimum tracklet cut. should be tuned
    return {beam_origin,
            {0, 0},
            {0, 0}};
  }

  delete xy_hist;  // if xy_hist is nullptr, nothing happen;
  xy_hist = new TH2F("xy_hist", "xy_hist",
                     N_bins, -1 * window_width / 2. + window_center.first,
                     window_width / 2. + window_center.first,
                     N_bins, -1 * window_width / 2. + window_center.second,
                     window_width / 2. + window_center.second);
  xy_hist->SetStats(false);
  xy_hist->GetXaxis()->SetTitle("X axis [mm]");
  xy_hist->GetYaxis()->SetTitle("Y axis [mm]");
  xy_hist->GetXaxis()->SetNdivisions(505);

  delete xy_hist_bkgrm;
  xy_hist_bkgrm = new TH2F("xy_hist_bkgrm", "xy_hist_bkgrm",
                           N_bins, -1 * window_width / 2. + window_center.first,
                           window_width / 2. + window_center.first,
                           N_bins, -1 * window_width / 2. + window_center.second,
                           window_width / 2. + window_center.second);
  xy_hist_bkgrm->SetStats(false);
  xy_hist_bkgrm->GetXaxis()->SetTitle("X axis [mm]");
  xy_hist_bkgrm->GetYaxis()->SetTitle("Y axis [mm]");
  xy_hist_bkgrm->GetXaxis()->SetNdivisions(505);

  // std::cout<<"test test size and bin of the hist xy_hist : "<<xy_hist -> GetNbinsX()<<" "<<xy_hist -> GetNbinsY()<<std::endl;
  // std::cout<<"test test bin width of the hist xy_hist : "<<xy_hist -> GetXaxis() -> GetBinWidth(1)<<" "<<xy_hist -> GetYaxis() -> GetBinWidth(1)<<std::endl;
  // std::cout<<"draw_plot status : "<<draw_plot<<std::endl;

  for (auto& i : cluster_pair_vec)
  {
    std::vector<double> DCA_info_vec = calculateDistanceAndClosestPoint(
        i.first.x, i.first.y,
        i.second.x, i.second.y,
        window_center.first, window_center.second);

    double DCA_sign = calculateAngleBetweenVectors(
        i.second.x, i.second.y,
        i.first.x, i.first.y,
        window_center.first, window_center.second);

    if (DCA_info_vec[0] != fabs(DCA_sign) && fabs(DCA_info_vec[0] - fabs(DCA_sign)) > 0.1)
    {
      std::cout << "different DCA : " << DCA_info_vec[0] << " " << DCA_sign << " diff : " << DCA_info_vec[0] - fabs(DCA_sign) << std::endl;
    }

    Clus_InnerPhi_Offset = (i.first.y - window_center.second < 0)
                               ? atan2(i.first.y - window_center.second,
                                       i.first.x - window_center.first) *
                                         (180. / M_PI) +
                                     360
                               : atan2(i.first.y - window_center.second,
                                       i.first.x - window_center.first) *
                                     (180. / M_PI);
    Clus_OuterPhi_Offset = (i.second.y - window_center.second < 0)
                               ? atan2(i.second.y - window_center.second,
                                       i.second.x - window_center.first) *
                                         (180. / M_PI) +
                                     360
                               : atan2(i.second.y - window_center.second,
                                       i.second.x - window_center.first) *
                                     (180. / M_PI);

    if (fabs(Clus_InnerPhi_Offset - Clus_OuterPhi_Offset) < 5)
    {
      TH2FSampleLineFill(xy_hist,
                         segmentation,
                         {i.first.x, i.first.y},
                         {DCA_info_vec[1], DCA_info_vec[2]});
      // note : the DCA cut may be biased since the DCA was calculated with respect to the vertex calculated by the quadrant method
      // if (DCA_cut.first < DCA_sign && DCA_sign < DCA_cut.second)
      // {
      //     TH2FSampleLineFill(xy_hist,
      //                        segmentation,
      //                        {cluster_pair_vec[i].first.x, cluster_pair_vec[i].first.y},
      //                        {DCA_info_vec[1], DCA_info_vec[2]}
      //                       );
      // }
    }
  }

  // note : try to remove the background
  TH2F_FakeClone(xy_hist, xy_hist_bkgrm);
  TH2F_threshold_advanced_2(xy_hist_bkgrm, 0.7);

  double reco_vtx_x = xy_hist_bkgrm->GetMean(1);  // note : the TH2F calculate the GetMean based on the bin center, no need to apply additional offset
  double reco_vtx_y = xy_hist_bkgrm->GetMean(2);  // note : the TH2F calculate the GetMean based on the bin center, no need to apply additional offset

  // std::cout<<"test : in the line filled, the process is almost done"<<std::endl;

  if (draw_plot)
  {
    TGraph* reco_vertex_gr = new TGraph();
    reco_vertex_gr->SetMarkerStyle(50);
    reco_vertex_gr->SetMarkerColor(2);
    reco_vertex_gr->SetMarkerSize(1);
    reco_vertex_gr->SetPoint(reco_vertex_gr->GetN(), reco_vtx_x, reco_vtx_y);

    std::string s_inttlabel = (boost::format("#it{#bf{sPHENIX INTT}} %s") % plot_text).str().c_str();
    // note : -----------------------------------------------------------------------------------------
    xy_hist->Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
    c1->Print((boost::format("%s/linefill_qa.pdf(") % out_folder_directory).str().c_str());
    c1->Clear();

    // note : -----------------------------------------------------------------------------------------
    xy_hist_bkgrm->Draw("colz0");
    draw_text->DrawLatex(0.21, 0.71 + 0.13, (boost::format("Vertex of the Run: %.4f mm, %.4f mm") % reco_vtx_x % reco_vtx_y).str().c_str());
    draw_text->DrawLatex(0.21, 0.67 + 0.13, (boost::format("Vertex error: %.4f mm, %.4f mm") % xy_hist_bkgrm->GetMeanError(1) % xy_hist_bkgrm->GetMeanError(2)).str().c_str());
    reco_vertex_gr->Draw("p same");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
    c1->Print((boost::format("%s/linefill_qa.pdf)") % out_folder_directory).str().c_str());
    c1->Clear();

    // std::cout<<"test : hello, can you see me ?"<<std::endl;
  }

  return {{reco_vtx_x, reco_vtx_y},
          {xy_hist_bkgrm->GetMeanError(1), xy_hist_bkgrm->GetMeanError(2)},
          {xy_hist_bkgrm->GetStdDev(1), xy_hist_bkgrm->GetStdDev(2)}};
}
