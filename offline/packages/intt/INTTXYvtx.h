#ifndef INTT_INTTXYVTX_H
#define INTT_INTTXYVTX_H

#include "InttVertexUtil.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TProfile.h>

#include <filesystem>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

class TH1;
class TH1F;
class TH2;

// note : this class mainly focus on two things
// note : 1. find a single vertex for the whole run
// note :     a. the functions prepared for this purpose are : ProcessEvt(), GetFinalVTXxy(), MacroVTXSquare()
// note : 2. find the vertex for each event
// note :     a. the functions prepared for this purpose are : ProcessEvt()

struct type_pos
{
  double x;
  double y;
};

class INTTXYvtx
{
 public:
  struct clu_info
  {
    int column;
    //--    int chip_id;
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

 public:
  INTTXYvtx(const std::string& runType,
            const std::string& outFolderDirectory,
            std::pair<double, double> beamOrigin,
            double phiDiffCut = 0.11,
            std::pair<double, double> DCACut = {-1, 1},
            int NCluCutl = 20,
            int NCluCut = 10000,
            double angleDiffNew_l = 0.0,
            double angleDiffNew_r = 3,
            double peekCut = 3.32405,
            bool printMessageOpt = true);

  virtual ~INTTXYvtx();

  void Init();

  void SetBeamOrigin(double beamx, double beamy)
  {
    beam_origin = std::make_pair(beamx, beamy);
    current_vtxX = beam_origin.first;
    current_vtxY = beam_origin.second;
  }

  virtual void SetSaveHisto(const bool save) { m_savehist = save; }
  virtual void EnableDrawHisto(const bool enable) { m_enable_drawhist = enable; }
  virtual void EnableQA(const bool enable) { m_enable_qa = enable; }
  virtual void PrintMessageOpt(const bool flag) { print_message_opt = flag; }

  //////////////////////////////////////////////////////
  // read cluster and make cluster pair in vector
  void ProcessEvt(int event_i,
                  const std::vector<clu_info>& temp_sPH_inner_nocolumn_vec,
                  const std::vector<clu_info>& temp_sPH_outer_nocolumn_vec,
                  const std::vector<std::vector<double>>& temp_sPH_nocolumn_vec,
                  const std::vector<std::vector<double>>& temp_sPH_nocolumn_rz_vec,
                  int NvtxMC,
                  double TrigZvtxMC,
                  bool PhiCheckTag,
                  Long64_t bco_full);

  //////////////////////////////////////////////////////
  // calculate XY vertex
  //  quadorant method
  std::vector<std::pair<double, double>> MacroVTXSquare(double length, int N_trial);
  // |- subMacroVTXxyCorrection(i,i1, draw_plot_opt);
  //     |- GetVTXxyCorrection_new(true_trial_index);
  //         |- subMacroPlotWorking
  //             |- calculateDistanceAndClosestPoint(
  //             |- calculateAngleBetweenVectors(
  //
  // |- DrawTGraphErrors(
  // |- Draw2TGraph(

  //  linefilled method
  std::vector<std::pair<double, double>> FillLine_FindVertex(
      std::pair<double, double> window_center,
      double segmentation = 0.005,
      double window_width = 5.0,
      int N_bins = 100);
  // |- calculateDistanceAndClosestPoint(
  // |- calculateAngleBetweenVectors(
  // |- TH2F_FakeClone(
  // |- TH2F_threshold_advanced_2(

  virtual void ClearEvt();
  virtual void PrintPlots();
  virtual void EndRun();

  //////////////////////////////////////////////////////
  // access to internal variables, necessary?
  unsigned long GetVecNele();
  std::pair<double, double> GetFinalVTXxy();
  std::pair<std::vector<TH2*>, std::vector<TH1F*>> GetHistFinal();

 private:
  std::string run_type;
  std::string out_folder_directory;
  std::pair<double, double> beam_origin;
//  double phi_diff_cut;
  std::pair<double, double> DCA_cut;
  int N_clu_cutl;
  int N_clu_cut;
  double angle_diff_new_l;
  double angle_diff_new_r;
  double peek;
  bool print_message_opt{false};

  bool m_savehist{false};
  bool m_enable_drawhist{false};
  bool m_enable_qa{false};

  bool m_initialized{false};
  std::string m_quad_pdfname{"New_Trial_square.pdf"};

  ////////////////////////////
  std::vector<double> subMacroVTXxyCorrection(int test_index, int trial_index, bool draw_plot_opt);
  std::vector<double> GetVTXxyCorrection_new(int trial_index);
  virtual void subMacroPlotWorking(bool phi_correction, double cos_fit_rangel, double cos_fit_ranger, double guas_fit_range);

  std::vector<double> calculateDistanceAndClosestPoint(double x1, double y1, double x2, double y2, double target_x, double target_y);
  double calculateAngleBetweenVectors(double x1, double y1, double x2, double y2, double targetX, double targetY);

  // void                    PrintPlotsVTXxy(std::string sub_out_folder_name);
  void PrintPlotsVTXxy();

  std::vector<std::pair<double, double>> Get4vtx(std::pair<double, double> origin, double length);
  void TH1F_FakeClone(TH1F* hist_in, TH1F* hist_out);
  void TH2F_FakeClone(TH2* hist_in, TH2* hist_out);

  void ClearHist(int print_option = 0);

  void DrawTGraphErrors(std::vector<double> x_vec, std::vector<double> y_vec,
                        std::vector<double> xE_vec, std::vector<double> yE_vec,
                        const std::string& output_directory, std::vector<std::string> plot_name);
  void Draw2TGraph(std::vector<double> x1_vec, std::vector<double> y1_vec,
                   std::vector<double> x2_vec, std::vector<double> y2_vec,
                   const std::string& output_directory, std::vector<std::string> plot_name);

  void TH2F_threshold(TH2* hist, double threshold);
  std::vector<double> SumTH2FColumnContent(TH2* hist_in);
  void Hist_1D_bkg_remove(TH1F* hist_in, double factor);
  void TH2F_threshold_advanced_2(TH2* hist, double threshold);

  // note : from the INTTXYvtxEvt.h // for FillLine_FindVertex
  void TH2FSampleLineFill(TH2* hist_in,
                          double segmentation,
                          std::pair<double, double>
                              inner_clu,
                          std::pair<double, double> outer_clu) const;

  std::vector<TH1*> m_v_hist{};

  /////////////////
  // QA histograms in process_evt // in m_v_hist
  TH2* N_cluster_correlation{nullptr};        // QA fill: ProcessEvt, draw: PrintPlot
  TH2* N_cluster_correlation_close{nullptr};  // QA fill: ProcessEvt, draw: PrintPlot
  TH2* inner_pos_xy{nullptr};                 // QA fill: ProcessEvt, draw: PrintPlot
  TH2* outer_pos_xy{nullptr};                 // QA fill: ProcessEvt, draw: PrintPlot
  TH2* inner_outer_pos_xy{nullptr};           // QA fill: ProcessEvt, draw: PrintPlot

  // Quadorant method // in m_v_hist
  TH2* DCA_distance_inner_phi{nullptr};  // fill: subMacroPlotWorking
  TH2* angle_diff_inner_phi{nullptr};    // fill: subMacroPlotWorking

  TH1F* angle_diff{nullptr};                 // QA fill: subMacroPlotWorking
  TH1F* angle_diff_new{nullptr};             // QA fill: subMacroPlotWorking
  TH1F* angle_diff_new_bkg_remove{nullptr};  // QA fill: subMacroPlotWorking
  TH1F* DCA_distance{nullptr};               // QA fill: subMacroPlotWorking

  TH2* DCA_distance_outer_phi{nullptr};  // QA fill: subMacroPlotWorking
  TH2* angle_diff_outer_phi{nullptr};    // QA fill: subMacroPlotWorking

  TH2* angle_correlation{nullptr};     // QA Fill: subMacroPlotWorking
  TH2* angle_diff_DCA_dist{nullptr};   // QA Fill: subMacroPlotWorking
  TH2* DCA_point{nullptr};             // QA fill: subMacroPlotWorking
  TH2* DCA_distance_inner_X{nullptr};  // QA fill: subMacroPlotWorking
  TH2* DCA_distance_inner_Y{nullptr};  // QA fill: subMacroPlotWorking
  TH2* DCA_distance_outer_X{nullptr};  // QA fill: subMacroPlotWorking
  TH2* DCA_distance_outer_Y{nullptr};  // QA fill: subMacroPlotWorking

  /// histograms & graphs created in subMacroPlotWorking
  TH2* DCA_distance_inner_phi_peak{nullptr};                  // fill: subMacroPlotWorking
  TProfile* DCA_distance_inner_phi_peak_profile{nullptr};      // fill: subMacroPlotWorking
  TGraph* DCA_distance_inner_phi_peak_profile_graph{nullptr};  // fill: subMacroPlotWorking

  TH2* angle_diff_inner_phi_peak{nullptr};                  // fill: subMacroPlotWorking
  TProfile* angle_diff_inner_phi_peak_profile{nullptr};      // fill: subMacroPlotWorking
  TGraph* angle_diff_inner_phi_peak_profile_graph{nullptr};  // fill: subMacroPlotWorking

  TH2* DCA_distance_outer_phi_peak{nullptr};                  // QA fill: subMacroPlotWorking
  TProfile* DCA_distance_outer_phi_peak_profile{nullptr};      // QA fill: subMacroPlotWorking
  TGraph* DCA_distance_outer_phi_peak_profile_graph{nullptr};  // QA fill: subMacroPlotWorking

  TH2* angle_diff_outer_phi_peak{nullptr};                  // QA fill: subMacroPlotWorking
  TProfile* angle_diff_outer_phi_peak_profile{nullptr};      // QA fill: subMacroPlotWorking
  TGraph* angle_diff_outer_phi_peak_profile_graph{nullptr};  // QA fill: subMacroPlotWorking

  // note : it's for the geometry correction // in m_v_hist
  TH2* angle_diff_inner_phi_peak_final{nullptr};    // fill: MacroVTXSquare
  TH2* DCA_distance_inner_phi_peak_final{nullptr};  // fill: MacroVTXSquare

  TH2* angle_diff_outer_phi_peak_final{nullptr};    // QA fill: MacroVTXSquare
  TH2* DCA_distance_outer_phi_peak_final{nullptr};  // QA fill: MacroVTXSquare
  TH1F* angle_diff_new_bkg_remove_final{nullptr};    // QA fill: MacroVTXSquare

  TF1* horizontal_fit_inner{nullptr};             // subMacroPlotWorking
  TF1* horizontal_fit_angle_diff_inner{nullptr};  // subMacroPlotWorking
  TF1* horizontal_fit_outer{nullptr};             // subMacroPlotWorking
  TF1* horizontal_fit_angle_diff_outer{nullptr};  // subMacroPlotWorking

  TF1* cos_fit{nullptr};   // subMacroPlotWorking, QA, not used for ana,
  TF1* gaus_fit{nullptr};  // subMacroPlotWorking, QA, not used for ana

  // LineFill method
  TH2* xy_hist{nullptr};        // fill: FillLine_FindVertex(
  TH2* xy_hist_bkgrm{nullptr};  // fill: FillLine_FindVertex(

  // note : to keep the cluster pair information
  // note : this is the vector for the whole run, not event by event
  std::vector<std::pair<type_pos, type_pos>> cluster_pair_vec{};

  double Clus_InnerPhi_Offset{0};
  double Clus_OuterPhi_Offset{0};
  double current_vtxX{0};
  double current_vtxY{0};

  int zvtx_cal_require = 15;

  std::string plot_text;
  long total_NClus{0};

  TCanvas* c1{nullptr};        // PrintPlotsVTXxy, PrintPlots(), FillLine_FindVertex, DrawTGraphErrors, Draw2TGraph
  TLatex* ltx{nullptr};        // PrintPlotsVTXxy, PrintPlots(), FillLine_FindVertex, DrawTGraphErrors, Draw2TGraph
  TLatex* draw_text{nullptr};  // PrintPlotsVTXxy,               FillLine_FindVertex, DrawTGraphErrors, Draw2TGraph
  /////////////////////

  void InitHist();
  void InitGraph();
  void InitRest();

  //--        void InitTreeOut();

  //--        // TFile * file_out;
  //--        // TTree * tree_out;
  //--        // double out_quadrant_corner_X;
  //--        // double out_quadrant_corner_Y;
  //--        // double out_quadrant_center_X;
  //--        // double out_quadrant_center_X;
  //--        // double out_line_filled_mean_X;
  //--        // double out_line_filled_mean_Y;
  //--        // double out_line_filled_stddev_X;
  //--        // double out_line_filled_stddev_Y;
};

#endif
