#ifndef INTT_INTTZVTX_H
#define INTT_INTTZVTX_H

#include "InttVertexUtil.h"

#include <cstdint>
#include <string>
#include <vector>

class TCanvas;
class TFile;
class TF1;
class TGraph;
class TGraphErrors;
class TH1;
class TH2;
class TLatex;
class TLine;
class TPad;
class TTree;

class INTTZvtx
{
 public:
  struct clu_info
  {
    int column;
    // int chip_id;
    double avg_chan;
    int sum_adc;
    int sum_adc_conv;
    int size;

    double x;
    double y;
    double z;
    int layer;
    double phi;
    // std::vector<double> bco_diff_vec; // note : for the multi-hit cluster, more than one hit was included. so more than one bco_diff
  };

  struct ZvtxInfo
  {
   public:
    double zvtx{-9999.};         // gaussian fit
    double zvtx_err{-1};         // gaussian fit
    double chi2ndf{-1};          // gaussian fit
    double width{-1};            // gaussian fit
    bool good{false};            // gaussian fit
    long nclus{0};               // Total_nclus
    unsigned int ntracklets{0};  // ngood tracklet N_comb
    unsigned int ngroup{0};      // ngroup in likebreak histogram
    double peakratio{-1};        // peak ratio in likebreak histogram
    double peakwidth{-1};        // peak width in likebreak histogram

    void clear()
    {
      zvtx = -9999.;   // gaussian fit
      zvtx_err = -1;   // gaussian fit
      chi2ndf = -1;    // gaussian fit
      width = -1;      // gaussian fit
      good = false;    // gaussian fit
      nclus = 0;       // Total_nclus
      ntracklets = 0;  // ngood tracklet N_comb
      ngroup = 0;      // ngroup in likebreak histogram
      peakratio = -1;  // peak ratio in likebreak histogram
      peakwidth = -1;  // peak width in likebreak histogram
    };
  };

 public:
  INTTZvtx(const std::string& runType,
           const std::string& outFolderDirectory,
           std::pair<double, double> beamOrigin,
           double phiDiffCut = 0.11,
           std::pair<double, double> DCACut = {-1, 1},
           int NCluCutl = 20,
           int NCluCut = 10000,
           unsigned int zvtxCal_require = 15,
           std::pair<double, double> zvtxQAWidth = {39.62, 65.36},
           bool drawEventDisplay = true,
           bool enableQA = true,
           bool printMessageOpt = true);

  virtual ~INTTZvtx();

  void Init();  // initialize all histograms and others

  bool ProcessEvt(int event_i,
                  std::vector<clu_info>& temp_sPH_inner_nocolumn_vec,
                  std::vector<clu_info>& temp_sPH_outer_nocolumn_vec,
                  std::vector<std::vector<double>>& temp_sPH_nocolumn_vec,
                  std::vector<std::vector<double>>& temp_sPH_nocolumn_rz_vec,
                  int NvtxMC,
                  double TrigZvtxMC,
                  bool PhiCheckTag,
                  uint64_t bco_full,
                  int centrality_bin);

  void ClearEvt();
  void PrintPlots();
  void EndRun();

  void EnableEventDisplay(const bool enableEvtDisp) { draw_event_display = enableEvtDisp; }
  void EnableQA(const bool enableQA) { m_enable_qa = enableQA; }

  double GetZdiffPeakMC();
  double GetZdiffWidthMC();

  std::vector<double> GetEvtZPeak();
  std::pair<double, double> GetBeamOrigin() { return beam_origin; }
  ZvtxInfo& GetZvtxInfo() { return m_zvtxinfo; }
  void SetBeamOrigin(double beamx, double beamy) { beam_origin = std::make_pair(beamx, beamy); }
  void SetPrintMessageOpt(const bool opt) { print_message_opt = opt; }
  void SetOutDirectory(const std::string& sOutDirectory) { out_folder_directory = sOutDirectory; }

 private:
  bool m_initialized{false};

  std::string run_type;
  std::string out_folder_directory;
  std::pair<double, double> beam_origin;
  double phi_diff_cut;                      // note : if (< phi_diff_cut)      -> pass      unit degree
  std::pair<double, double> DCA_cut;        // note : if (< DCA_cut)           -> pass      unit mm
  int N_clu_cut;                            // note : if (> N_clu_cut)         -> continue  unit number
  int N_clu_cutl;                           // note : if (< N_clu_cutl)        -> continue  unit number
  unsigned int zvtx_cal_require;            // note : if (> zvtx_cal_require)  -> pass
  std::pair<double, double> zvtx_QA_width;  // note : for the zvtx range Quality check, check the width
  bool draw_event_display{false};
  bool m_enable_qa{false};
  bool print_message_opt;

  std::pair<double, double> evt_possible_z_range = {-700, 700};

  std::vector<std::string> conversion_mode_BD = {"ideal", "survey_1_XYAlpha_Peek", "full_survey_3.32"};
  double Integrate_portion_final = 0.68;  // cut in effSig, PrintPlots
  double Integrate_portion = 0.35;        // cut in effSig, processEvt, todo : Width selection per event (the range finder for fitting)
  double high_multi_line = 1000;          // cut in ProcessEvt, todo : the cut to classify the high-low multiplicity, which are fit with different method
  double zvtx_hist_l = -500;              // histogram range for QA
  double zvtx_hist_r = 500;               // histogram range for QA
  int print_rate = 50;                    // if_print in processEvt, todo : the print rate is here

  std::vector<std::vector<std::pair<bool, clu_info>>> inner_clu_phi_map{};  // note: phi
  std::vector<std::vector<std::pair<bool, clu_info>>> outer_clu_phi_map{};  // note: phi

  ZvtxInfo m_zvtxinfo;

  TH1* evt_possible_z{nullptr};
  TH1* line_breakdown_hist{nullptr};  // note : try to fill the line into the histogram
  TF1* gaus_fit{nullptr};
  TF1* zvtx_finder{nullptr};
  TGraphErrors* z_range_gr{nullptr};  // ana // memory leak

  //////////////////////////////////////
  // for event display
  TH1* evt_select_track_phi{nullptr};    // ProcessEvt
  TH1* evt_phi_diff_1D{nullptr};         // ProcessEvt
  TH2* evt_phi_diff_inner_phi{nullptr};  // ProcessEvt
  TH2* evt_inner_outer_phi{nullptr};     // ProcessEvt

  // for QA
  TH1* avg_event_zvtx{nullptr};                  // Fill: ProcessEvt, Draw: PrintPlot
  TH1* zvtx_evt_fitError{nullptr};               // Fill: ProcessEvt, Draw: PrintPlot
  TH2* zvtx_evt_fitError_corre{nullptr};         // Fill: ProcessEvt, Draw: PrintPlot
  TH2* zvtx_evt_width_corre{nullptr};            // Fill: ProcessEvt, Draw: PrintPlot
  TH2* zvtx_evt_nclu_corre{nullptr};             // Fill: ProcessEvt, Draw: PrintPlot
  TH1* width_density{nullptr};                   // note : N good hits / width // Fill: ProcessEvt, Draw: PrintPlot
  TH1* ES_width{nullptr};                        // Fill: ProcessEvt, Draw: PrintPlot
  TH1* ES_width_ratio{nullptr};                  // Fill: ProcessEvt, Draw: PrintPlot
  TH2* Z_resolution_Nclu{nullptr};               // Fill: ProcessEvt, Draw: PrintPlot
  TH2* Z_resolution_pos{nullptr};                // Fill: ProcessEvt, Draw: PrintPlot
  TH2* Z_resolution_pos_cut{nullptr};            // Fill: ProcessEvt, Draw: PrintPlot
  TH1* Z_resolution{nullptr};                    // Fill: ProcessEvt, Draw: PrintPlot
  TH1* line_breakdown_gaus_ratio_hist{nullptr};  // note : the distribution of the entry/width of gaus fit // Fill: ProcessEvt, Draw: PrintPlot
  TH1* line_breakdown_gaus_width_hist{nullptr};  // note : the distribution of the gaus fit width // Fill: ProcessEvt, Draw: PrintPlot
  TH2* gaus_width_Nclu{nullptr};                 // Fill: ProcessEvt, Draw: PrintPlot
  TH2* gaus_rchi2_Nclu{nullptr};                 // Fill: ProcessEvt, Draw: PrintPlot
  TH2* final_fit_width{nullptr};                 // Fill: ProcessEvt, Draw: PrintPlot
  TH2* N_track_candidate_Nclu{nullptr};          // note : Number of tracklet candidate (In xy plane) vs number of clusters// Fill: ProcessEvt, Draw: PrintPlot
  TH1* peak_group_width_hist{nullptr};           // Fill: ProcessEvt, Draw: PrintPlot
  TH1* peak_group_ratio_hist{nullptr};           // Fill: ProcessEvt, Draw: PrintPlot
  TH1* N_group_hist{nullptr};                    // Fill: ProcessEvt, Draw: PrintPlot
  TH1* peak_group_detail_width_hist{nullptr};    // Fill: ProcessEvt, Draw: PrintPlot
  TH1* peak_group_detail_ratio_hist{nullptr};    // Fill: ProcessEvt, Draw: PrintPlot
  TH1* N_group_detail_hist{nullptr};             // Fill: ProcessEvt, Draw: PrintPlot

  TH2* phi_diff_inner_phi{nullptr};  // ProcessEvt
  TH2* dca_inner_phi{nullptr};       // ProcessEvt

  std::vector<TH1*> m_v_qahist{};

  TCanvas* c1{nullptr};  // PrintPlots

  TCanvas* c2{nullptr};                // processEvt
  TPad* pad_xy{nullptr};               // ProcessEvt
  TPad* pad_rz{nullptr};               // ProcessEvt
  TPad* pad_z{nullptr};                // ProcessEvt
  TPad* pad_z_hist{nullptr};           // ProcessEvt
  TPad* pad_z_line{nullptr};           // ProcessEvt
  TPad* pad_phi_diff{nullptr};         // ProcessEvt
  TPad* pad_track_phi{nullptr};        // ProcessEvt
  TPad* pad_inner_outer_phi{nullptr};  // ProcessEvt
  TPad* pad_phi_diff_1D{nullptr};      // ProcessEvt

  TLine* ladder_line{nullptr};           // tempbkg
  TLine* final_fit_range_line{nullptr};  // ProcessEvt
  TLine* coord_line{nullptr};            // ProcessEvt
  TLatex* draw_text{nullptr};            // ProcessEvt, PrintPlots
  TLine* eff_sig_range_line{nullptr};    // ProcessEvt, PrintPlots

  TFile* out_file{nullptr};
  TTree* tree_out{nullptr};

  // note : for tree_out
  double out_ES_zvtx, out_ES_zvtxE, out_ES_rangeL, out_ES_rangeR, out_ES_width_density, MC_true_zvtx;
  double out_LB_Gaus_Mean_mean, out_LB_Gaus_Mean_meanE, out_LB_Gaus_Mean_width, out_LB_Gaus_Mean_chi2;
  double out_LB_Gaus_Width_width, out_LB_Gaus_Width_size_width, out_LB_Gaus_Width_offset, out_LB_geo_mean;
  double out_mid_cut_peak_width, out_mid_cut_peak_ratio, out_LB_cut_peak_width, out_LB_cut_peak_ratio;
  bool out_good_zvtx_tag;
  int out_eID, N_cluster_outer_out, N_cluster_inner_out, out_ES_N_good, out_mid_cut_Ngroup, out_LB_cut_Ngroup, out_centrality_bin;
  int out_N_cluster_north, out_N_cluster_south;
  uint64_t bco_full_out;

  // note : for out parameters
  double MC_z_diff_peak, MC_z_diff_width;  // note : the comparison between Reco - true in MC. Values are from fitting foucsing on the peak region.

  void InitHist();
  void InitCanvas();
  void InitTreeOut();
  void InitRest();

  std::vector<float> temp_event_zvtx_info;  // effSig method
  std::vector<float> avg_event_zvtx_vec;    // for QA
  std::vector<float> Z_resolution_vec;      // for QA for MC, zvtx diff btw reso and MC
  std::vector<double> N_group_info;         // QA,  note : the information of groups remaining in the histogram after the strong background suppression
  std::vector<double> N_group_info_detail;  // good_vtx note : detail

  double final_zvtx{0};
  double tight_offset_width{0};  // tight zvertex QA & draw
  double tight_offset_peak{0};   // tight zvertex QA & draw
  double loose_offset_peak{0};   // z-vertex!
  double loose_offset_peakE{0};  // z-vertex error!
  bool good_zvtx_tag{0};
  double good_zvtx_tag_int{0};
  int good_comb_id{0};  // used for tracklet reco

  TGraphErrors* z_range_gr_draw{nullptr};  // draw processEvt
  TGraph* temp_event_xy{nullptr};          // draw processEvt
  TGraph* temp_event_rz{nullptr};          // draw processEvt
  TGraph* bkg{nullptr};                    // draw in tempbkg

  // note : in the event process
  std::vector<float> N_comb{};       // tracklet
  std::vector<float> N_comb_e{};     // tracklet
  std::vector<double> N_comb_phi{};  // tracklet
  std::vector<float> z_mid{};        // tracklet
  std::vector<float> z_range{};      // tracklet

  // function for analysis
  std::pair<double, double> Get_possible_zvtx(double rvtx, std::vector<double> p0, std::vector<double> p1);
  std::vector<double> find_Ngroup(TH1* hist_in);
  double get_radius(double x, double y);
  double calculateAngleBetweenVectors(double x1, double y1, double x2, double y2, double targetX, double targetY);
  double Get_extrapolation(double given_y, double p0x, double p0y, double p1x, double p1y);
  void line_breakdown(TH1* hist_in, std::pair<double, double> line_range);

  // tracklet reco
  double get_delta_phi(double angle_1, double angle_2);
  double get_track_phi(double inner_clu_phi_in, double delta_phi_in);

  // for Tree
  double LB_geo_mean(TH1* hist_in, std::pair<double, double> search_range, int event_i);

  // InitCanvas
  void Characterize_Pad(TPad* pad, float left = 0.15, float right = 0.1,
                        float top = 0.1, float bottom = 0.12,
                        bool set_logY = false, int setgrid_bool = 0);
};

#endif
