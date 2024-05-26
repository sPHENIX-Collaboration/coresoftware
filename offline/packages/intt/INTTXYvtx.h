#ifndef INTTXYvtx_h
#define INTTXYvtx_h

#include "InttVertexUtil.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <numeric>
#include <filesystem>


#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TFile.h>

using std::string;
using std::vector;
using std::pair;
using std::cout;
using std::endl;

// note : this class mainly focus on two things 
// note : 1. find a single vertex for the whole run
// note :     a. the functions prepared for this purpose are : ProcessEvt(), GetFinalVTXxy(), MacroVTXSquare()
// note : 2. find the vertex for each event
// note :     a. the functions prepared for this purpose are : ProcessEvt()

struct type_pos{
    double x;
    double y;
};

double cos_func(double *x, double *par)
{
    return -1 * par[0] * cos(par[1] * (x[0] - par[2])) + par[3];
}

class INTTXYvtx {

  public:
    struct clu_info {
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
        //std::vector<double> bco_diff_vec; // note : for the multi-hit cluster, more than one hit was included. so more than one bco_diff
    };

    public:
        INTTXYvtx(string               runType, 
                  string               outFolderDirectory,
                  pair<double, double> beamOrigin, 
                  double               phiDiffCut       = 0.11, 
                  pair<double, double> DCACut           = {-1,1}, 
                  int                  NCluCutl         = 20, 
                  int                  NCluCut          = 10000, 
                  double               angleDiffNew_l   = 0.0, 
                  double               angleDiffNew_r   = 3, 
                  double               peekCut          = 3.32405, 
                  bool                 printMessageOpt  = true);

        virtual ~INTTXYvtx();

        void Init();

        void SetBeamOrigin(double beamx, double beamy) 
                  { beam_origin  = std::make_pair(beamx, beamy); 
                    current_vtxX = beam_origin.first;
                    current_vtxY = beam_origin.second;
                  }


        virtual void SetSaveHisto(   const bool save)   { m_savehist        = save; }
        virtual void EnableDrawHisto(const bool enable) { m_enable_drawhist = enable; }
        virtual void EnableQA(       const bool enable) { m_enable_qa       = enable; }
        virtual void PrintMessageOpt(const bool flag)   { print_message_opt = flag; }
 
        //////////////////////////////////////////////////////
        // read cluster and make cluster pair in vector
        void ProcessEvt(int                    event_i, 
                        vector<clu_info>       temp_sPH_inner_nocolumn_vec, 
                        vector<clu_info>       temp_sPH_outer_nocolumn_vec, 
                        vector<vector<double>> temp_sPH_nocolumn_vec, 
                        vector<vector<double>> temp_sPH_nocolumn_rz_vec, 
                        int                    NvtxMC, 
                        double                 TrigZvtxMC, 
                        bool                   PhiCheckTag, 
                        Long64_t               bco_full);
        
        //////////////////////////////////////////////////////
        // calculate XY vertex
        //  quadorant method
        vector<pair<double,double>> MacroVTXSquare(double length, int N_trial);
                                      // |- subMacroVTXxyCorrection(i,i1, draw_plot_opt);
                                      //     |- GetVTXxyCorrection_new(true_trial_index);
                                      //         |- subMacroPlotWorking
                                      //             |- calculateDistanceAndClosestPoint(
                                      //             |- calculateAngleBetweenVectors(
                                      //
                                      // |- DrawTGraphErrors(
                                      // |- Draw2TGraph(

        //  linefilled method
        vector<pair<double,double>> FillLine_FindVertex(
                                        pair<double,double> window_center, 
                                        double              segmentation = 0.005, 
                                        double              window_width = 5.0, 
                                        int                 N_bins       = 100
                                        );
                                      // |- calculateDistanceAndClosestPoint(
                                      // |- calculateAngleBetweenVectors(
                                      // |- TH2F_FakeClone(
                                      // |- TH2F_threshold_advanced_2(

        virtual void ClearEvt();
        virtual void PrintPlots();
        virtual void EndRun();

        //////////////////////////////////////////////////////
        // access to internal variables, necessary?
        unsigned long                       GetVecNele();
        pair<double,double>                 GetFinalVTXxy();
        pair<vector<TH2F *>, vector<TH1F*>> GetHistFinal();

    protected:
        string              run_type;
        string              out_folder_directory;
        pair<double,double> beam_origin;
        double              phi_diff_cut;        
        pair<double,double> DCA_cut;
        int                 N_clu_cutl;
        int                 N_clu_cut;
        double              angle_diff_new_l;
        double              angle_diff_new_r;
        double              peek;
        bool                print_message_opt{false};

        bool                m_savehist{false};
        bool                m_enable_drawhist{false};
        bool                m_enable_qa{false};

        bool                m_initialized{false};
        string              m_quad_pdfname{"New_Trial_square.pdf"};

    protected:
        ////////////////////////////
        vector<double>          subMacroVTXxyCorrection(int test_index, int trial_index, bool draw_plot_opt);
        vector<double>          GetVTXxyCorrection_new(int trial_index);
        virtual void            subMacroPlotWorking(bool phi_correction, double cos_fit_rangel, double cos_fit_ranger, double guas_fit_range);

        std::vector<double>     calculateDistanceAndClosestPoint(double x1, double y1, double x2, double y2, double target_x, double target_y);
        double                  calculateAngleBetweenVectors(double x1, double y1, double x2, double y2, double targetX, double targetY);


        //void                    PrintPlotsVTXxy(string sub_out_folder_name);
        void                    PrintPlotsVTXxy();

        vector<pair<double,double>> Get4vtx(pair<double,double> origin, double length);
        void                        TH1F_FakeClone(TH1F* hist_in, TH1F* hist_out);
        void                        TH2F_FakeClone(TH2F* hist_in, TH2F* hist_out);

        void                        ClearHist(int print_option = 0);


        void                    DrawTGraphErrors(vector<double> x_vec, vector<double> y_vec, 
                                                 vector<double> xE_vec, vector<double> yE_vec, 
                                                 string output_directory, vector<string> plot_name);
        void                    Draw2TGraph(vector<double> x1_vec, vector<double> y1_vec, 
                                            vector<double> x2_vec, vector<double> y2_vec, 
                                            string output_directory, vector<string> plot_name);

        void                    TH2F_threshold(TH2F * hist, double threshold);
        vector<double>          SumTH2FColumnContent(TH2F * hist_in);
        void                    Hist_1D_bkg_remove(TH1F * hist_in, double factor);
        void                    TH2F_threshold_advanced_2(TH2F * hist, double threshold);

        // note : from the INTTXYvtxEvt.h // for FillLine_FindVertex
        void                    TH2FSampleLineFill(TH2F * hist_in, 
                                                   double segmentation, 
                                                   std::pair<double,double> 
                                                   inner_clu, std::pair<double,double> outer_clu);
        


    protected:
        std::vector<TH1*> m_v_hist{};

        /////////////////
        // QA histograms in process_evt // in m_v_hist
        TH2F*    N_cluster_correlation      {nullptr};//QA fill: ProcessEvt, draw: PrintPlot
        TH2F*    N_cluster_correlation_close{nullptr};//QA fill: ProcessEvt, draw: PrintPlot
        TH2F*    inner_pos_xy               {nullptr};//QA fill: ProcessEvt, draw: PrintPlot
        TH2F*    outer_pos_xy               {nullptr};//QA fill: ProcessEvt, draw: PrintPlot
        TH2F*    inner_outer_pos_xy         {nullptr};//QA fill: ProcessEvt, draw: PrintPlot

        // Quadorant method // in m_v_hist
        TH2F*     DCA_distance_inner_phi{nullptr};//fill: subMacroPlotWorking
        TH2F*     angle_diff_inner_phi  {nullptr};//fill: subMacroPlotWorking

        TH1F*     angle_diff               {nullptr};//QA fill: subMacroPlotWorking
        TH1F*     angle_diff_new           {nullptr};//QA fill: subMacroPlotWorking
        TH1F*     angle_diff_new_bkg_remove{nullptr};//QA fill: subMacroPlotWorking
        TH1F*     DCA_distance             {nullptr};//QA fill: subMacroPlotWorking

        TH2F*     DCA_distance_outer_phi{nullptr};//QA fill: subMacroPlotWorking
        TH2F*     angle_diff_outer_phi  {nullptr};//QA fill: subMacroPlotWorking

        TH2F*     angle_correlation     {nullptr};//QA Fill: subMacroPlotWorking
        TH2F*     angle_diff_DCA_dist   {nullptr};//QA Fill: subMacroPlotWorking
        TH2F*     DCA_point             {nullptr};//QA fill: subMacroPlotWorking
        TH2F*     DCA_distance_inner_X  {nullptr};//QA fill: subMacroPlotWorking
        TH2F*     DCA_distance_inner_Y  {nullptr};//QA fill: subMacroPlotWorking
        TH2F*     DCA_distance_outer_X  {nullptr};//QA fill: subMacroPlotWorking
        TH2F*     DCA_distance_outer_Y  {nullptr};//QA fill: subMacroPlotWorking

        /// histograms & graphs created in subMacroPlotWorking
        TH2F*     DCA_distance_inner_phi_peak              {nullptr};//fill: subMacroPlotWorking
        TProfile* DCA_distance_inner_phi_peak_profile      {nullptr};//fill: subMacroPlotWorking
        TGraph*   DCA_distance_inner_phi_peak_profile_graph{nullptr};//fill: subMacroPlotWorking

        TH2F*     angle_diff_inner_phi_peak                {nullptr};//fill: subMacroPlotWorking
        TProfile* angle_diff_inner_phi_peak_profile        {nullptr};//fill: subMacroPlotWorking
        TGraph*   angle_diff_inner_phi_peak_profile_graph  {nullptr};//fill: subMacroPlotWorking

        TH2F*     DCA_distance_outer_phi_peak              {nullptr};//QA fill: subMacroPlotWorking
        TProfile* DCA_distance_outer_phi_peak_profile      {nullptr};//QA fill: subMacroPlotWorking
        TGraph*   DCA_distance_outer_phi_peak_profile_graph{nullptr};//QA fill: subMacroPlotWorking

        TH2F*     angle_diff_outer_phi_peak                {nullptr};//QA fill: subMacroPlotWorking
        TProfile* angle_diff_outer_phi_peak_profile        {nullptr};//QA fill: subMacroPlotWorking
        TGraph*   angle_diff_outer_phi_peak_profile_graph  {nullptr};//QA fill: subMacroPlotWorking

        // note : it's for the geometry correction // in m_v_hist
        TH2F*    angle_diff_inner_phi_peak_final  {nullptr};//fill: MacroVTXSquare
        TH2F*    DCA_distance_inner_phi_peak_final{nullptr};//fill: MacroVTXSquare

        TH2F*    angle_diff_outer_phi_peak_final  {nullptr};//QA fill: MacroVTXSquare
        TH2F*    DCA_distance_outer_phi_peak_final{nullptr};//QA fill: MacroVTXSquare
        TH1F*    angle_diff_new_bkg_remove_final  {nullptr};//QA fill: MacroVTXSquare

        TF1*     horizontal_fit_inner           {nullptr};//subMacroPlotWorking
        TF1*     horizontal_fit_angle_diff_inner{nullptr};//subMacroPlotWorking
        TF1*     horizontal_fit_outer           {nullptr};//subMacroPlotWorking
        TF1*     horizontal_fit_angle_diff_outer{nullptr};//subMacroPlotWorking

        TF1*     cos_fit                        {nullptr};//subMacroPlotWorking, QA, not used for ana, 
        TF1*     gaus_fit                       {nullptr};//subMacroPlotWorking, QA, not used for ana
                

        // LineFill method
        TH2F*    xy_hist      {nullptr};//fill: FillLine_FindVertex(
        TH2F*    xy_hist_bkgrm{nullptr};//fill: FillLine_FindVertex(
        
        // note : to keep the cluster pair information
        // note : this is the vector for the whole run, not event by event
        vector<pair<type_pos,type_pos>> cluster_pair_vec{};

        double   Clus_InnerPhi_Offset{0};
        double   Clus_OuterPhi_Offset{0};
        double   current_vtxX{0};
        double   current_vtxY{0};

        int      zvtx_cal_require = 15;

        string   plot_text;
        long     total_NClus;

        TCanvas* c1       {nullptr}; //PrintPlotsVTXxy, PrintPlots(), FillLine_FindVertex, DrawTGraphErrors, Draw2TGraph
        TLatex*  ltx      {nullptr}; //PrintPlotsVTXxy, PrintPlots(), FillLine_FindVertex, DrawTGraphErrors, Draw2TGraph
        TLatex*  draw_text{nullptr}; //PrintPlotsVTXxy,               FillLine_FindVertex, DrawTGraphErrors, Draw2TGraph
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

INTTXYvtx::INTTXYvtx(string              runType, 
                     string              outFolderDirectory, 
                     pair<double,double> beamOrigin, 
                     double              phiDiffCut, 
                     pair<double,double> DCACut, 
                     int                 NCluCutl, 
                     int                 NCluCut, 
                     double              angleDiffNew_l, 
                     double              angleDiffNew_r, 
                     double              peekCut, 
                     bool                printMessageOpt)
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
    gErrorIgnoreLevel = kWarning; // note : To not print the "print plot info."

    //Init();
    plot_text = (run_type == "MC") ? "Simulation" : "Work-in-progress";

    cluster_pair_vec.clear();

    current_vtxX = beam_origin.first;
    current_vtxY = beam_origin.second;

}

INTTXYvtx::~INTTXYvtx()
{
  // InitHist
  for(auto& itr: m_v_hist){
    //--cout<<"del : "<<itr->GetTitle()<<" "<<std::hex<<(long)itr<<std::hex<<endl;
    delete itr; itr=nullptr;
  }
  //--for(auto& itr: m_v_hist){
  //--  cout<<"after del : "<<std::hex<<(long)itr<<std::hex<<endl;
  //--}

  // InitGraph
  delete angle_diff_inner_phi_peak_profile_graph;   angle_diff_inner_phi_peak_profile_graph  =nullptr;
  delete angle_diff_outer_phi_peak_profile_graph;   angle_diff_outer_phi_peak_profile_graph  =nullptr;
  delete DCA_distance_inner_phi_peak_profile_graph; DCA_distance_inner_phi_peak_profile_graph=nullptr;
  delete DCA_distance_outer_phi_peak_profile_graph; DCA_distance_outer_phi_peak_profile_graph=nullptr;

  // created in FillLine_FindVertex
  delete xy_hist;       xy_hist      =nullptr;       
  delete xy_hist_bkgrm; xy_hist_bkgrm=nullptr;
  	
  // InitRest
  delete cos_fit;                         cos_fit                        =nullptr; 
  delete gaus_fit;                        gaus_fit                       =nullptr;  
  delete horizontal_fit_inner;            horizontal_fit_inner           =nullptr;
  delete horizontal_fit_angle_diff_inner; horizontal_fit_angle_diff_inner=nullptr;
  delete horizontal_fit_outer;            horizontal_fit_outer           =nullptr;
  delete horizontal_fit_angle_diff_outer; horizontal_fit_angle_diff_outer=nullptr;

  delete c1;        c1       = nullptr;
  delete ltx;       ltx      = nullptr;
  delete draw_text; draw_text= nullptr;
}

void INTTXYvtx::Init()
{
    if(m_enable_drawhist||m_savehist){
      if ( !std::filesystem::exists(out_folder_directory.c_str()) ) {
        system(Form("mkdir %s",out_folder_directory.c_str()));
      }
    }

    InitHist();
    InitGraph();
    InitRest();

//--    InitTreeOut();

    m_initialized=true;
}

//--void INTTXYvtx::InitTreeOut()
//--{
//--    // file_out = new TFile(Form("%s/run_XY_tree.root",out_folder_directory.c_str()),"RECREATE");
//--//--    // file_out -> cd();
//--//--
//--//--    // tree_out = new TTree("tree", "tree avg VtxXY");
//--//--    // tree_out -> Branch("quadrant_corner_X",&out_quadrant_corner_X);
//--//--    // tree_out -> Branch("quadrant_corner_Y",&out_quadrant_corner_Y);
//--//--    // tree_out -> Branch("quadrant_center_X",&out_quadrant_center_X);
//--//--    // tree_out -> Branch("quadrant_center_X",&out_quadrant_center_X);
//--//--    // tree_out -> Branch("line_filled_mean_X",&out_line_filled_mean_X);
//--//--    // tree_out -> Branch("line_filled_mean_Y",&out_line_filled_mean_Y);
//--//--    // tree_out -> Branch("line_filled_stddev_X",&out_line_filled_stddev_X);
//--//--    // tree_out -> Branch("line_filled_stddev_Y",&out_line_filled_stddev_Y);
//--
//--    return;
//--}

void INTTXYvtx::InitGraph()
{
    angle_diff_inner_phi_peak_profile_graph = new TGraph();
    angle_diff_outer_phi_peak_profile_graph = new TGraph();
    DCA_distance_inner_phi_peak_profile_graph = new TGraph();
    DCA_distance_outer_phi_peak_profile_graph = new TGraph();
}

void INTTXYvtx::InitRest()
{
    horizontal_fit_inner = new TF1("horizontal_fit_inner","pol0",0,360);
    horizontal_fit_inner -> SetLineWidth(2);
    horizontal_fit_inner -> SetLineColor(2);

    horizontal_fit_angle_diff_inner = new TF1("horizontal_fit_angle_diff_inner","pol0",0,360);
    horizontal_fit_angle_diff_inner -> SetLineWidth(2);
    horizontal_fit_angle_diff_inner -> SetLineColor(2);

    horizontal_fit_outer = new TF1("horizontal_fit_outer","pol0",0,360);
    horizontal_fit_outer -> SetLineWidth(2);
    horizontal_fit_outer -> SetLineColor(2);

    horizontal_fit_angle_diff_outer = new TF1("horizontal_fit_angle_diff_outer","pol0",0,360);
    horizontal_fit_angle_diff_outer -> SetLineWidth(2);
    horizontal_fit_angle_diff_outer -> SetLineColor(2);

    if(m_enable_qa){
      cos_fit = new TF1("cos_fit",cos_func,0,360,4);
      cos_fit -> SetParNames("[A]", "[B]", "[C]", "[D]");
      cos_fit -> SetLineColor(2);

      gaus_fit = new TF1("gaus_fit", InttVertexUtil::gaus_func,0,360,4);
      gaus_fit -> SetLineColor(4);
      gaus_fit -> SetLineWidth(1);
      gaus_fit -> SetParNames("size", "mean", "width", "offset");
      gaus_fit -> SetNpx(1000);
    }

    if(m_enable_drawhist) 
    {
      c1 = new TCanvas("","",950,800);
      c1 -> cd();

      ltx = new TLatex();
      ltx->SetNDC();
      ltx->SetTextSize(0.045);
      ltx->SetTextAlign(31);

      draw_text = new TLatex();
      draw_text -> SetNDC();
      draw_text -> SetTextSize(0.03);
    }
}

void INTTXYvtx::InitHist()
{
    ///////////////////////////////////////
    // QA histograms in process_evt
    if(m_enable_qa){
      N_cluster_correlation = new TH2F("N_cluster_correlation","N_cluster_correlation",100,0,4000,100,0,4000);
      N_cluster_correlation -> SetStats(0);
      N_cluster_correlation -> GetXaxis() -> SetTitle("inner N_cluster");
      N_cluster_correlation -> GetYaxis() -> SetTitle("Outer N_cluster");
      N_cluster_correlation -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(N_cluster_correlation);

      N_cluster_correlation_close = new TH2F("N_cluster_correlation_close","N_cluster_correlation_close",100,0,500,100,0,500);
      N_cluster_correlation_close -> SetStats(0);
      N_cluster_correlation_close -> GetXaxis() -> SetTitle("inner N_cluster");
      N_cluster_correlation_close -> GetYaxis() -> SetTitle("Outer N_cluster");
      N_cluster_correlation_close -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(N_cluster_correlation_close);

      inner_pos_xy = new TH2F("inner_pos_xy","inner_pos_xy",360,-100,100,360,-100,100);
      inner_pos_xy -> SetStats(0);
      inner_pos_xy -> GetXaxis() -> SetTitle("X axis [mm]");
      inner_pos_xy -> GetYaxis() -> SetTitle("Y axis [mm]");
      inner_pos_xy -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(inner_pos_xy);

      outer_pos_xy = new TH2F("outer_pos_xy","outer_pos_xy",360,-150,150,360,-150,150);
      outer_pos_xy -> SetStats(0);
      outer_pos_xy -> GetXaxis() -> SetTitle("X axis [mm]");
      outer_pos_xy -> GetYaxis() -> SetTitle("Y axis [mm]");
      outer_pos_xy -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(outer_pos_xy);

      //inner_outer_pos_xy = new TH2F("inner_outer_pos_xy","inner_outer_pos_xy",360,-150,150,360,-150,150);
      inner_outer_pos_xy = new TH2F("inner_outer_pos_xy","inner_outer_pos_xy",360,-150,150,360,-150,150);
      inner_outer_pos_xy -> SetStats(0);
      inner_outer_pos_xy -> GetXaxis() -> SetTitle("X axis [mm]");
      inner_outer_pos_xy -> GetYaxis() -> SetTitle("Y axis [mm]");
      inner_outer_pos_xy -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(inner_outer_pos_xy);
    }

    ///////////////////////////////////////
    // histograms for quadrant method
    DCA_distance_inner_phi = new TH2F("DCA_distance_inner_phi","DCA_distance_inner_phi",100,0,360,100,-10,10);
    DCA_distance_inner_phi -> SetStats(0);
    DCA_distance_inner_phi -> GetXaxis() -> SetTitle("Inner phi [degree]");
    DCA_distance_inner_phi -> GetYaxis() -> SetTitle("DCA [mm]");
    DCA_distance_inner_phi -> GetXaxis() -> SetNdivisions(505);
    m_v_hist.push_back(DCA_distance_inner_phi);

    angle_diff_inner_phi = new TH2F("angle_diff_inner_phi","angle_diff_inner_phi",361,0,361,100,-1.5,1.5);
    angle_diff_inner_phi -> SetStats(0);
    angle_diff_inner_phi -> GetXaxis() -> SetTitle("Inner phi [degree]");
    angle_diff_inner_phi -> GetYaxis() -> SetTitle("Inner - Outer [degree]");
    angle_diff_inner_phi -> GetXaxis() -> SetNdivisions(505);
    m_v_hist.push_back(angle_diff_inner_phi);

    ///////////////////////////////////////
    // note : it's for the geometry correction
    angle_diff_inner_phi_peak_final = new TH2F("angle_diff_inner_phi_peak_final",
                                               "angle_diff_inner_phi_peak_final",361,0,361,100,-1.5,1.5);
    angle_diff_inner_phi_peak_final -> SetStats(0);
    angle_diff_inner_phi_peak_final -> GetXaxis() -> SetTitle("Inner phi [degree]");
    angle_diff_inner_phi_peak_final -> GetYaxis() -> SetTitle("Inner - Outer [degree]");
    angle_diff_inner_phi_peak_final -> GetXaxis() -> SetNdivisions(505);
    m_v_hist.push_back(angle_diff_inner_phi_peak_final);

    DCA_distance_inner_phi_peak_final = new TH2F("DCA_distance_inner_phi_peak_final",
                                                 "DCA_distance_inner_phi_peak_final",100,0,360,100,-10,10);
    DCA_distance_inner_phi_peak_final -> SetStats(0);
    DCA_distance_inner_phi_peak_final -> GetXaxis() -> SetTitle("Inner phi [degree]");
    DCA_distance_inner_phi_peak_final -> GetYaxis() -> SetTitle("DCA [mm]");
    DCA_distance_inner_phi_peak_final -> GetXaxis() -> SetNdivisions(505);
    m_v_hist.push_back(DCA_distance_inner_phi_peak_final);

    // QA
      angle_diff = new TH1F("angle_diff","angle_diff",100,0,5);
      angle_diff -> SetStats(0);
      angle_diff -> GetXaxis() -> SetTitle("|Inner - Outer| [degree]");
      angle_diff -> GetYaxis() -> SetTitle("Entry");
      angle_diff -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(angle_diff);

      angle_diff_new = new TH1F("angle_diff_new","angle_diff_new",100, angle_diff_new_l, angle_diff_new_r);
      angle_diff_new -> SetStats(0);
      angle_diff_new -> GetXaxis() -> SetTitle("|Inner - Outer| [degree]");
      angle_diff_new -> GetYaxis() -> SetTitle("Entry");
      angle_diff_new -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(angle_diff_new);

      DCA_distance = new TH1F("DCA_distance","DCA_distance",100,-10,10);
      DCA_distance -> SetStats(0);
      DCA_distance -> GetXaxis() -> SetTitle("DCA [mm]");
      DCA_distance -> GetYaxis() -> SetTitle("Entry");
      DCA_distance -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(DCA_distance);


    if(m_enable_qa){
      DCA_distance_outer_phi = new TH2F("DCA_distance_outer_phi","DCA_distance_outer_phi",100,0,360,100,-10,10);
      DCA_distance_outer_phi -> SetStats(0);
      DCA_distance_outer_phi -> GetXaxis() -> SetTitle("Outer phi [degree]");
      DCA_distance_outer_phi -> GetYaxis() -> SetTitle("DCA [mm]");
      DCA_distance_outer_phi -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(DCA_distance_outer_phi);

      angle_diff_outer_phi = new TH2F("angle_diff_outer_phi","angle_diff_outer_phi",361,0,361,100,-1.5,1.5);
      angle_diff_outer_phi -> SetStats(0);
      angle_diff_outer_phi -> GetXaxis() -> SetTitle("Outer phi [degree]");
      angle_diff_outer_phi -> GetYaxis() -> SetTitle("Inner - Outer [degree]");
      angle_diff_outer_phi -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(angle_diff_outer_phi);

      angle_correlation = new TH2F("angle_correlation","angle_correlation",361,0,361,361,0,361);
      angle_correlation -> SetStats(0);
      angle_correlation -> GetXaxis() -> SetTitle("Inner Phi [degree]");
      angle_correlation -> GetYaxis() -> SetTitle("Outer Phi [degree]");
      angle_correlation -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(angle_correlation);

      angle_diff_DCA_dist = new TH2F("angle_diff_DCA_dist","angle_diff_DCA_dist",100,-1.5,1.5,100,-3.5,3.5);
      angle_diff_DCA_dist -> SetStats(0);
      angle_diff_DCA_dist -> GetXaxis() -> SetTitle("Inner - Outer [degree]");
      angle_diff_DCA_dist -> GetYaxis() -> SetTitle("DCA distance [mm]");
      angle_diff_DCA_dist -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(angle_diff_DCA_dist);

      DCA_point = new TH2F("DCA_point","DCA_point",100,-10,10,100,-10,10);
      DCA_point -> SetStats(0);
      DCA_point -> GetXaxis() -> SetTitle("X pos [mm]");
      DCA_point -> GetYaxis() -> SetTitle("Y pos [mm]");
      DCA_point -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(DCA_point);

      DCA_distance_inner_X = new TH2F("DCA_distance_inner_X","DCA_distance_inner_X",100,-100,100,100,-10,10);
      DCA_distance_inner_X -> SetStats(0);
      DCA_distance_inner_X -> GetXaxis() -> SetTitle("Inner cluster X [mm]");
      DCA_distance_inner_X -> GetYaxis() -> SetTitle("DCA [mm]");
      DCA_distance_inner_X -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(DCA_distance_inner_X);

      DCA_distance_inner_Y = new TH2F("DCA_distance_inner_Y","DCA_distance_inner_Y",100,-100,100,100,-10,10);
      DCA_distance_inner_Y -> SetStats(0);
      DCA_distance_inner_Y -> GetXaxis() -> SetTitle("Inner cluster Y [mm]");
      DCA_distance_inner_Y -> GetYaxis() -> SetTitle("DCA [mm]");
      DCA_distance_inner_Y -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(DCA_distance_inner_Y);

      DCA_distance_outer_X = new TH2F("DCA_distance_outer_X","DCA_distance_outer_X",100,-130,130,100,-10,10);
      DCA_distance_outer_X -> SetStats(0);
      DCA_distance_outer_X -> GetXaxis() -> SetTitle("Outer cluster X [mm]");
      DCA_distance_outer_X -> GetYaxis() -> SetTitle("DCA [mm]");
      DCA_distance_outer_X -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(DCA_distance_outer_X);

      DCA_distance_outer_Y = new TH2F("DCA_distance_outer_Y","DCA_distance_outer_Y",100,-130,130,100,-10,10);
      DCA_distance_outer_Y -> SetStats(0);
      DCA_distance_outer_Y -> GetXaxis() -> SetTitle("Outer cluster Y [mm]");
      DCA_distance_outer_Y -> GetYaxis() -> SetTitle("DCA [mm]");
      DCA_distance_outer_Y -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(DCA_distance_outer_Y);

      ///////////////////////////////////////
      // note : it's for the geometry correction
      angle_diff_outer_phi_peak_final = new TH2F("angle_diff_outer_phi_peak_final",
                                                 "angle_diff_outer_phi_peak_final",361,0,361,100,-1.5,1.5);
      angle_diff_outer_phi_peak_final -> SetStats(0);
      angle_diff_outer_phi_peak_final -> GetXaxis() -> SetTitle("Outer phi [degree]");
      angle_diff_outer_phi_peak_final -> GetYaxis() -> SetTitle("Inner - Outer [degree]");
      angle_diff_outer_phi_peak_final -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(angle_diff_outer_phi_peak_final);

      DCA_distance_outer_phi_peak_final = new TH2F("DCA_distance_outer_phi_peak_final",
                                                   "DCA_distance_outer_phi_peak_final",100,0,360,100,-10,10);
      DCA_distance_outer_phi_peak_final -> SetStats(0);
      DCA_distance_outer_phi_peak_final -> GetXaxis() -> SetTitle("Outer phi [degree]");
      DCA_distance_outer_phi_peak_final -> GetYaxis() -> SetTitle("DCA [mm]");
      DCA_distance_outer_phi_peak_final -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(DCA_distance_outer_phi_peak_final);

      angle_diff_new_bkg_remove_final = new TH1F("angle_diff_new_bkg_remove_final",
                                                 "angle_diff_new_bkg_remove_final",100, angle_diff_new_l, angle_diff_new_r);
      angle_diff_new_bkg_remove_final -> SetStats(0);
      angle_diff_new_bkg_remove_final -> GetXaxis() -> SetTitle("|Inner - Outer| [degree]");
      angle_diff_new_bkg_remove_final -> GetYaxis() -> SetTitle("Entry");
      angle_diff_new_bkg_remove_final -> GetXaxis() -> SetNdivisions(505);
      m_v_hist.push_back(angle_diff_new_bkg_remove_final);
    }
}

// note : this function only prepare the pairs for the vertex XY calculation, it's like a general vertex for the whole run
void INTTXYvtx::ProcessEvt(
    int                    event_i,
    vector<clu_info>       temp_sPH_inner_nocolumn_vec,
    vector<clu_info>       temp_sPH_outer_nocolumn_vec, 
    vector<vector<double>> /*temp_sPH_nocolumn_vec*/,
    vector<vector<double>> /*temp_sPH_nocolumn_rz_vec*/,
    int                    NvtxMC,
    double                 /*TrigZvtxMC*/,
    bool                   PhiCheckTag,
    Long64_t               /*bco_full*/
)
{
    if(!m_initialized) {
       cout<<"INTTXYvtx is not initialized, abort in ProcessEvt"<<endl;
       exit(1);
    }

    if (print_message_opt && event_i%10000 == 0) {cout<<"In INTTXYvtx class, running event : "<<event_i<<endl;}

    total_NClus = temp_sPH_inner_nocolumn_vec.size() + temp_sPH_outer_nocolumn_vec.size();

    // note : the Move these two in the beginning of the function, in order to avoid those event-reject cuts
    if(m_enable_qa){
      N_cluster_correlation -> Fill(temp_sPH_inner_nocolumn_vec.size(),temp_sPH_outer_nocolumn_vec.size());
      N_cluster_correlation_close -> Fill(temp_sPH_inner_nocolumn_vec.size(),temp_sPH_outer_nocolumn_vec.size());
    }

    if (total_NClus < zvtx_cal_require) {
       return; 
       cout<<"return confirmation"<<endl;
    }
    
    if (run_type == "MC" && NvtxMC != 1) { 
       return; 
       cout<<"In INTTXYvtx class, event : "<<event_i
           <<" Nvtx : "<<NvtxMC<<" Nvtx more than one "<<endl;
    }
    if (PhiCheckTag == false) { 
       return; 
       cout<<"In INTTXYvtx class, event : "<<event_i
           <<" Nvtx : "<<NvtxMC<<" Not full phi has hits "<<endl;
    }
    
    if (   temp_sPH_inner_nocolumn_vec.size() < 10 
        || temp_sPH_outer_nocolumn_vec.size() < 10 
        || total_NClus > N_clu_cut 
        || total_NClus < N_clu_cutl)
    {
        return;
        printf("In INTTXYvtx class, event : %i, low clu continue, NClus : %lu \n", event_i, total_NClus); 
    }

    //-------------------------------
    // tracklet reconstruction accumulated multiple events
    for ( unsigned int inner_i = 0; inner_i < temp_sPH_inner_nocolumn_vec.size(); inner_i++ )
    {    
        for ( unsigned int outer_i = 0; outer_i < temp_sPH_outer_nocolumn_vec.size(); outer_i++ )
        {
            // note : try to ease the analysis and also make it quick.
            if (fabs(temp_sPH_inner_nocolumn_vec[inner_i].phi - temp_sPH_outer_nocolumn_vec[outer_i].phi) < 7) // todo : the pre phi cut is here, can be optimized
            {
                cluster_pair_vec.push_back({{temp_sPH_inner_nocolumn_vec[inner_i].x,
                                             temp_sPH_inner_nocolumn_vec[inner_i].y},
                                            {temp_sPH_outer_nocolumn_vec[outer_i].x,
                                             temp_sPH_outer_nocolumn_vec[outer_i].y}});
            }
            
        }
    }

    //-------------------------------
    // QA histogram
    if(m_enable_qa){
      for ( unsigned int inner_i = 0; inner_i < temp_sPH_inner_nocolumn_vec.size(); inner_i++ )
      {
          inner_pos_xy -> Fill(temp_sPH_inner_nocolumn_vec[inner_i].x,temp_sPH_inner_nocolumn_vec[inner_i].y);
          inner_outer_pos_xy -> Fill(temp_sPH_inner_nocolumn_vec[inner_i].x,temp_sPH_inner_nocolumn_vec[inner_i].y);
      }

      for ( unsigned int outer_i = 0; outer_i < temp_sPH_outer_nocolumn_vec.size(); outer_i++ )
      {
          outer_pos_xy -> Fill(temp_sPH_outer_nocolumn_vec[outer_i].x,temp_sPH_outer_nocolumn_vec[outer_i].y);
          inner_outer_pos_xy -> Fill(temp_sPH_outer_nocolumn_vec[outer_i].x,temp_sPH_outer_nocolumn_vec[outer_i].y);
      }
    }

    //--  cout<<"  "<<event_i<<" clusterpair:size : "<<cluster_pair_vec.size()<<endl;
}

void INTTXYvtx::ClearEvt()
{
    return;
}


unsigned long INTTXYvtx::GetVecNele()
{
    return cluster_pair_vec.size();
}

pair<double,double> INTTXYvtx::GetFinalVTXxy()
{
    return {current_vtxX, current_vtxY};
}

pair<vector<TH2F *>, vector<TH1F*>> INTTXYvtx::GetHistFinal()
{
    return {
        {DCA_distance_inner_phi_peak_final, 
         angle_diff_inner_phi_peak_final, 
         DCA_distance_outer_phi_peak_final, 
         angle_diff_outer_phi_peak_final, 
         xy_hist, 
         xy_hist_bkgrm},
        {angle_diff_new_bkg_remove_final}
    };
}

void INTTXYvtx::PrintPlots()
{
    if(!m_initialized) {
       cout<<"INTTXYvtx is not initialized, abort in PrintPlots"<<endl;
       exit(1);
    }

    if(m_enable_drawhist&&m_enable_qa)
    {
      string s_inttlabel = Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str());
      // note : -----------------------------------------------------------------------------------------
      inner_outer_pos_xy -> Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01,  s_inttlabel.c_str());
      //c1 -> Print(Form("%s/inner_outer_pos_xy.pdf",out_folder_directory.c_str()));
      c1 -> Print(Form("%s/xyvtx_qa.pdf(", out_folder_directory.c_str()));
      c1 -> Clear();
      
      // note : -----------------------------------------------------------------------------------------
      inner_pos_xy -> Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01,  s_inttlabel.c_str());
      //c1 -> Print(Form("%s/inner_pos_xy.pdf",out_folder_directory.c_str()));
      c1 -> Print(Form("%s/xyvtx_qa.pdf", out_folder_directory.c_str()));
      c1 -> Clear();

      // note : -----------------------------------------------------------------------------------------
      outer_pos_xy -> Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01,  s_inttlabel.c_str());
      //c1 -> Print(Form("%s/outer_pos_xy.pdf",out_folder_directory.c_str()));
      c1 -> Print(Form("%s/xyvtx_qa.pdf", out_folder_directory.c_str()));
      c1 -> Clear();

      // note : -----------------------------------------------------------------------------------------
      N_cluster_correlation -> Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01,  s_inttlabel.c_str());
      //c1 -> Print(Form("%s/N_cluster_correlation.pdf",out_folder_directory.c_str()));
      c1 -> Print(Form("%s/xyvtx_qa.pdf", out_folder_directory.c_str()));
      c1 -> Clear();

      // note : -----------------------------------------------------------------------------------------
      N_cluster_correlation_close -> Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01,  s_inttlabel.c_str());
      //c1 -> Print(Form("%s/N_cluster_correlation_close.pdf",out_folder_directory.c_str()));
      c1 -> Print(Form("%s/xyvtx_qa.pdf)", out_folder_directory.c_str()));
      c1 -> Clear();
    }
}


//--------------------------------------------------------------------------
// quadrant method
vector<pair<double,double>> INTTXYvtx::MacroVTXSquare(double length, int N_trial)
{
    if(!m_initialized) {
       cout<<"INTTXYvtx is not initialized, abort in MacroVTXSquare"<<endl;
       exit(1);
    }

    bool draw_plot_opt = m_enable_drawhist;

    const double                original_length = length;
    pair<double,double>         origin          = {0,0};
    vector<pair<double,double>> vtx_vec         = Get4vtx(origin,length); // vtx_vec.push_back(origin);

    int            small_index{0};
    vector<double> small_info_vec(18, -999);
    vector<double> grr_x{};
    vector<double> grr_E{};
    vector<double> grr_y{};

    vector<double> All_FitError_DCA_Y{};  
    vector<double> All_FitError_DCA_X{};  
    vector<double> All_FitError_angle_Y{};
    vector<double> All_FitError_angle_X{};

    vector<double> Winner_FitError_DCA_Y{};  
    vector<double> Winner_FitError_DCA_X{};  
    vector<double> Winner_FitError_angle_Y{};
    vector<double> Winner_FitError_angle_X{};


    if (print_message_opt == true) {
       cout<<"In INTTXYvtx::MacroVTXSquare, N pairs : "<<cluster_pair_vec.size()<<endl;

       cout<<N_trial<<" runs, smart. which gives you the resolution down to "
           <<length/pow(2,N_trial)<<" mm"<<endl;
    }

    if(cluster_pair_vec.size()==0){ // minimum tracklet cut. need to be tuned
       return {
           beam_origin,          // note : the best vertex 
           origin,               // note : the origin in that trial
           {0, 0}, // note : horizontal_fit_inner -> GetParError(0),  horizontal_fit_angle_diff_inner -> GetParError(0)
           {0, 0}, // note : horizontal_fit_inner -> GetParameter(0), horizontal_fit_angle_diff_inner -> GetParameter(0)
           {0, 0}, // note : horizontal_fit_outer -> GetParError(0),  horizontal_fit_angle_diff_outer -> GetParError(0)
           {0, 0}, // note : horizontal_fit_outer -> GetParameter(0), horizontal_fit_angle_diff_outer -> GetParameter(0)
           {0, 0}, // note : the mean and stddev of angle_diff
           {0, 0}, // note : the mean and stddev of DCA_distance
           {0, 0}, // note : the mean and stddev of angle_diff, but with the background removed
       };
    }

    //-------------------------------
    //info_vec contents
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

    if(m_enable_drawhist){
      c1->cd();
      c1->Range(0,0,1,1);
      ltx->DrawLatex(0.5, 0.5,  Form("QA plots for quadrant method"));
      c1 -> Print( Form("%s/%s(", out_folder_directory.c_str(), m_quad_pdfname.c_str()) );
      c1 -> Clear(); 
    }

    // current algorithm uses info_vec[3] and [5], others are for QA
    for (int i = 0; i < N_trial; i++)
    {
        if (print_message_opt == true) {
          cout<<"~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~"
              <<" step "<<i<<" "
              <<"~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~"<<endl;
        }
        for (unsigned int i1 = 0; i1 < vtx_vec.size(); i1++)
        {
            if (print_message_opt == true) {cout<<"tested vertex : "<<vtx_vec[i1].first<<" "<<vtx_vec[i1].second<<endl;}

            if(m_enable_drawhist){
              c1->cd();
              c1->Range(0,0,1,1);
              ltx->DrawLatex(0.5, 0.5,  Form("New_trial_square_%i_%i",
                                               i, //test_index, 
                                               i1)); //trial_index
              c1 -> Print(Form("%s/%s", out_folder_directory.c_str(), m_quad_pdfname.c_str()));
              c1 -> Clear(); 
            }

            current_vtxX = vtx_vec[i1].first;
            current_vtxY = vtx_vec[i1].second;

            vector<double> info_vec = subMacroVTXxyCorrection(i,i1, draw_plot_opt);

 
            if (print_message_opt == true) {
              cout<<"trial : "<<i
                  <<" vertex : "<<i1
                  <<" DCA fit error : "<<info_vec[3]
                  <<" angle diff fit error : "<<info_vec[5]<<endl;
            }

            All_FitError_DCA_Y.push_back(info_vec[3]);
            All_FitError_DCA_X.push_back(i);
            All_FitError_angle_Y.push_back(info_vec[5]);
            All_FitError_angle_X.push_back(i);

            if ( (i1 == 0) ||
                 (info_vec[3] < small_info_vec[3] && 
                  info_vec[5] < small_info_vec[5]
                 ) // note : the fit error of the pol0 fit
               )
            {
                small_info_vec = info_vec;
                small_index    = i1;
                
                TH2F_FakeClone(DCA_distance_inner_phi_peak, DCA_distance_inner_phi_peak_final);
                TH2F_FakeClone(angle_diff_inner_phi_peak,   angle_diff_inner_phi_peak_final);
                if(m_enable_qa){
                  TH2F_FakeClone(DCA_distance_outer_phi_peak, DCA_distance_outer_phi_peak_final);
                  TH2F_FakeClone(angle_diff_outer_phi_peak,   angle_diff_outer_phi_peak_final);
                  TH1F_FakeClone(angle_diff_new_bkg_remove,   angle_diff_new_bkg_remove_final);
                }
            }
            if (print_message_opt == true){cout<<" "<<endl;}

            ClearHist(1);
        }

        if (print_message_opt == true) {cout<<"the Quadrant "<<small_index<<" won the competition"<<endl;}
        
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
            origin = {(vtx_vec[small_index].first + origin.first)/2., 
                      (vtx_vec[small_index].second + origin.second)/2.};

            // cout<<"test : "<<origin.first<<" "<<origin.second<<" length: "<<length<<endl;
            // if (small_index == 4) {length /= 1.5;}
            // else {length /= 2.;}
            length /= 2.;
            vtx_vec = Get4vtx(origin,length); // vtx_vec.push_back(origin);
        }
    }

    if(m_enable_drawhist){
      c1->cd();
      c1 -> Print(Form("%s/%s)", out_folder_directory.c_str(), m_quad_pdfname.c_str()));
      c1 -> Clear(); 
    }

    if (draw_plot_opt == true) {
       DrawTGraphErrors(grr_x, grr_y, grr_E, grr_E, out_folder_directory, 
                        {Form("Square_scan_history_%.1fmm_%iTrials", original_length, N_trial) // title
                         , "nth scan"          // x_title
                         , "Winner index"      // y_title
                         , "APL"               // draw option
                         , "quadorant_qa.pdf(" // pdf name
                        });
       Draw2TGraph(All_FitError_angle_X, All_FitError_angle_Y, Winner_FitError_angle_X, Winner_FitError_angle_Y, out_folder_directory, 
                        {Form("Angle_diff_fit_error_%iTrials", N_trial) // title
                         , "n iteration"                   // x_title
                         , "#Delta#phi fit error [degree]" // y_title
                         , ""            // draw option    // draw option
                         , "quadorant_qa.pdf"// pdf name   // pdf name
                        });
       Draw2TGraph(All_FitError_DCA_X, All_FitError_DCA_Y, Winner_FitError_DCA_X, Winner_FitError_DCA_Y, out_folder_directory, 
                        {Form("DCA_fit_error_%iTrials", N_trial) // title
                         , "n iteration"        // x_title     
                         , "DCA fit error [mm]" // y_title
                         , ""                   // draw option
                         , "quadorant_qa.pdf)"  // pdf name   
                        });
    }

    return {
        vtx_vec[small_index], // note : the best vertex 
        origin,               // note : the origin in that trial
        {small_info_vec[3], small_info_vec[5]},  // note : horizontal_fit_inner -> GetParError(0),  horizontal_fit_angle_diff_inner -> GetParError(0)
        {small_info_vec[10],small_info_vec[11]}, // note : horizontal_fit_inner -> GetParameter(0), horizontal_fit_angle_diff_inner -> GetParameter(0)
        {small_info_vec[7], small_info_vec[9]},  // note : horizontal_fit_outer -> GetParError(0),  horizontal_fit_angle_diff_outer -> GetParError(0)
        {small_info_vec[12],small_info_vec[13]}, // note : horizontal_fit_outer -> GetParameter(0), horizontal_fit_angle_diff_outer -> GetParameter(0)
        {small_info_vec[14],small_info_vec[15]}, // note : the mean and stddev of angle_diff
        {small_info_vec[16],small_info_vec[17]}, // note : the mean and stddev of DCA_distance
        {small_info_vec[0], small_info_vec[1]},  // note : the mean and stddev of angle_diff, but with the background removed
    };
}

vector<double> INTTXYvtx::subMacroVTXxyCorrection(int test_index, int trial_index, bool draw_plot_opt)
{
    int true_trial_index = test_index * 4 + trial_index;
    vector<double> out_vec = GetVTXxyCorrection_new(true_trial_index);

    string sub_out_folder_name{};
    if (draw_plot_opt == true){
        sub_out_folder_name = Form("%s/New_trial_square_%i_%i",
                                    out_folder_directory.c_str(), 
                                    test_index, 
                                    trial_index);

        //--if (std::filesystem::exists(sub_out_folder_name.c_str()) == false) {
        //--  system(Form("mkdir %s",sub_out_folder_name.c_str()));
        //--}

        //PrintPlotsVTXxy(sub_out_folder_name);
        PrintPlotsVTXxy();
    }

    return out_vec;
}

// note : {circle radius, possible correction angle, the chi2/NDF of pol0 fit}
vector<double> INTTXYvtx::GetVTXxyCorrection_new(int trial_index)
{
    if (print_message_opt == true) {
        cout<<"Trial : "<<trial_index
            <<"---------------------------- ---------------------------- ----------------------------"<<endl;
        cout<<"Given vertex: "<<current_vtxX <<" "<<current_vtxY<<endl;
    }

    if(m_enable_qa){
      cos_fit -> SetParameters(4, 1.49143e-02,  185, 0.3); // todo : we may have to apply more constaints on the fitting
      cos_fit -> SetParLimits(2,0,360); // note : the peak location has to be positive

      // note : here is the test with a gaus fitting to find the peak
      gaus_fit -> SetParameters(-4.5, 197, 50, 0);
      gaus_fit -> SetParLimits(0,-100,0); // note : the gaus distribution points down
      // DCA_distance_inner_phi_peak_profile -> Fit(gaus_fit, "N","",100, 260);
      // cout<<"test, gaus fit range : "<<gaus_fit->GetParameter(1) - 25<<" "<<gaus_fit->GetParameter(1) + 25<<endl;
    }
    
    subMacroPlotWorking(1,100,260,25);

    return {
        angle_diff_new_bkg_remove -> GetMean(),
        angle_diff_new_bkg_remove -> GetStdDev(),           // note : angle diff stddev and error (1D histogram)
        horizontal_fit_inner -> GetChisquare() / double(horizontal_fit_inner -> GetNDF()), 
        horizontal_fit_inner            -> GetParError(0),               // note : inner DCA, pol0
        horizontal_fit_angle_diff_inner -> GetChisquare() / double(horizontal_fit_angle_diff_inner -> GetNDF()), 
        horizontal_fit_angle_diff_inner -> GetParError(0),    // note : inner angle diff, pol0
        horizontal_fit_outer -> GetChisquare() / double(horizontal_fit_outer -> GetNDF()), 
        horizontal_fit_outer -> GetParError(0),               // note : outer DCA, pol0
        horizontal_fit_angle_diff_outer -> GetChisquare() / double(horizontal_fit_angle_diff_outer -> GetNDF()), 
        horizontal_fit_angle_diff_outer -> GetParError(0),    // note : outer angle diff, pol0
        horizontal_fit_inner            -> GetParameter(0), 
        horizontal_fit_angle_diff_inner -> GetParameter(0), // note : 10, 11
        horizontal_fit_outer            -> GetParameter(0), 
        horizontal_fit_angle_diff_outer -> GetParameter(0), // note : 12, 13
        angle_diff                      -> GetMean()  , 
        angle_diff                      -> GetStdDev(), // note : 14, 15
        DCA_distance                    -> GetMean()  , 
        DCA_distance                    -> GetStdDev(), // note : 16, 17
    };
}

void INTTXYvtx::subMacroPlotWorking(
              bool   phi_correction, 
              double cos_fit_rangel, 
              double cos_fit_ranger, 
              double guas_fit_range)
{
    
    //   3 : horizontal_fit_inner            -> ParErr(0),   // note : inner DCA, pol0
    //   5 : horizontal_fit_angle_diff_inner -> ParErr(0),   // note : inner angle diff, pol0

    for (unsigned int i = 0; i < cluster_pair_vec.size(); i++)
    {
        vector<double> DCA_info_vec = calculateDistanceAndClosestPoint(
            cluster_pair_vec[i].first.x, cluster_pair_vec[i].first.y,
            cluster_pair_vec[i].second.x, cluster_pair_vec[i].second.y,
            current_vtxX, current_vtxY
        );

        double DCA_sign = calculateAngleBetweenVectors(
            cluster_pair_vec[i].second.x, cluster_pair_vec[i].second.y,
            cluster_pair_vec[i].first.x, cluster_pair_vec[i].first.y,
            current_vtxX, current_vtxY
        );

        if (phi_correction == true)
        {
            // cout<<"option selected "<<endl;
            Clus_InnerPhi_Offset = (cluster_pair_vec[i].first.y - current_vtxY < 0) 
                               ? atan2(cluster_pair_vec[i].first.y - current_vtxY, cluster_pair_vec[i].first.x - current_vtxX) * (180./TMath::Pi()) + 360 
                               : atan2(cluster_pair_vec[i].first.y - current_vtxY, cluster_pair_vec[i].first.x - current_vtxX) * (180./TMath::Pi());
            Clus_OuterPhi_Offset = (cluster_pair_vec[i].second.y - current_vtxY < 0) 
                               ? atan2(cluster_pair_vec[i].second.y - current_vtxY, cluster_pair_vec[i].second.x - current_vtxX) * (180./TMath::Pi()) + 360 
                               : atan2(cluster_pair_vec[i].second.y - current_vtxY, cluster_pair_vec[i].second.x - current_vtxX) * (180./TMath::Pi());
        }
        else // note : phi_correction == false 
        {
            Clus_InnerPhi_Offset = (cluster_pair_vec[i].first.y < 0) 
                               ? atan2(cluster_pair_vec[i].first.y, cluster_pair_vec[i].first.x) * (180./TMath::Pi()) + 360 
                               : atan2(cluster_pair_vec[i].first.y, cluster_pair_vec[i].first.x) * (180./TMath::Pi()); 
            Clus_OuterPhi_Offset = (cluster_pair_vec[i].second.y < 0) 
                               ? atan2(cluster_pair_vec[i].second.y, cluster_pair_vec[i].second.x) * (180./TMath::Pi()) + 360 
                               : atan2(cluster_pair_vec[i].second.y, cluster_pair_vec[i].second.x) * (180./TMath::Pi()); 
        }

        //----------------------
        // this is used for quadrant method
        DCA_distance_inner_phi-> Fill(Clus_InnerPhi_Offset, DCA_sign); 
        angle_diff_inner_phi  -> Fill(Clus_InnerPhi_Offset, Clus_InnerPhi_Offset - Clus_OuterPhi_Offset);

        //----------------------
        if(m_enable_qa){
          DCA_distance_outer_phi-> Fill(Clus_OuterPhi_Offset, DCA_sign);
          angle_diff_outer_phi  -> Fill(Clus_OuterPhi_Offset, Clus_InnerPhi_Offset - Clus_OuterPhi_Offset);
    
          angle_diff            -> Fill(abs(Clus_InnerPhi_Offset - Clus_OuterPhi_Offset));
          angle_diff_new        -> Fill(abs(Clus_InnerPhi_Offset - Clus_OuterPhi_Offset));
          DCA_distance          -> Fill(DCA_sign);

          // draw only
          angle_correlation     -> Fill(Clus_InnerPhi_Offset, Clus_OuterPhi_Offset);
          angle_diff_DCA_dist   -> Fill(Clus_InnerPhi_Offset - Clus_OuterPhi_Offset, DCA_sign);
          DCA_point             -> Fill(DCA_info_vec[1], DCA_info_vec[2]);

          DCA_distance_inner_X  -> Fill(cluster_pair_vec[i].first.x, DCA_sign);
          DCA_distance_inner_Y  -> Fill(cluster_pair_vec[i].first.y, DCA_sign);
          DCA_distance_outer_X  -> Fill(cluster_pair_vec[i].second.x, DCA_sign);
          DCA_distance_outer_Y  -> Fill(cluster_pair_vec[i].second.y, DCA_sign);
        }

    } // note : end of the loop for the cluster pair

    // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
    delete DCA_distance_inner_phi_peak;
    DCA_distance_inner_phi_peak = (TH2F*)DCA_distance_inner_phi -> Clone("DCA_distance_inner_phi_peak");
    TH2F_threshold(DCA_distance_inner_phi_peak, 0.5); // todo : the background cut can be modified, the ratio 0.5
    delete DCA_distance_inner_phi_peak_profile;
    DCA_distance_inner_phi_peak_profile = DCA_distance_inner_phi_peak->ProfileX("DCA_distance_inner_phi_peak_profile");
    
    double point_index = 0;
    vector<double> hist_column_content = SumTH2FColumnContent(DCA_distance_inner_phi_peak);
    for (int i = 0; i < DCA_distance_inner_phi_peak_profile->GetNbinsX(); i++){
        if (hist_column_content[i] < 5){continue;} // note : in order to remove some remaining background

        DCA_distance_inner_phi_peak_profile_graph -> SetPoint(point_index, 
                                                              DCA_distance_inner_phi_peak_profile->GetBinCenter(i+1), 
                                                              DCA_distance_inner_phi_peak_profile->GetBinContent(i+1));
        // cout<<"("<<DCA_distance_inner_phi_peak_profile->GetBinCenter(i+1)<<", "<< DCA_distance_inner_phi_peak_profile->GetBinContent(i+1)<<")"<<endl;
        point_index += 1;
    }

    //------------------------------------------------------------------
    // this is used to constrain the quadrant
    // info_vec[3];
    // todo : the fit range of the gaussian fit can be modified here 
    //
    horizontal_fit_inner -> SetParameter(0,0);
    DCA_distance_inner_phi_peak_profile_graph -> Fit(horizontal_fit_inner,"NQ","",0,360);
    //------------------------------------------------------------------

    if(m_enable_qa){
      DCA_distance_inner_phi_peak_profile_graph -> Fit(gaus_fit, "NQ","",
                                                       cos_fit->GetParameter(2) - guas_fit_range, // note : what we want and need is the peak position, 
                                                       cos_fit->GetParameter(2) + guas_fit_range);// so we fit the peak again 
      DCA_distance_inner_phi_peak_profile_graph -> Fit(cos_fit, "NQ","",cos_fit_rangel, cos_fit_ranger);
    }

    // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
    delete angle_diff_inner_phi_peak;
    angle_diff_inner_phi_peak = (TH2F*)angle_diff_inner_phi -> Clone("angle_diff_inner_phi_peak");
    TH2F_threshold_advanced_2(angle_diff_inner_phi_peak, 0.5); // todo : threshold ratio can be modified here
    hist_column_content = SumTH2FColumnContent(angle_diff_inner_phi_peak);
    angle_diff_inner_phi_peak_profile = angle_diff_inner_phi_peak->ProfileX("angle_diff_inner_phi_peak_profile");
    point_index = 0;
    for (int i = 0; i < angle_diff_inner_phi_peak_profile->GetNbinsX(); i++){
        if (hist_column_content[i] < 5){continue;} // note : in order to remove some remaining background
        
        angle_diff_inner_phi_peak_profile_graph -> SetPoint(point_index, 
                                                            angle_diff_inner_phi_peak_profile->GetBinCenter(i+1), 
                                                            angle_diff_inner_phi_peak_profile->GetBinContent(i+1));
        // cout<<"("<<angle_diff_inner_phi_peak_profile->GetBinCenter(i+1)<<", "<< angle_diff_inner_phi_peak_profile->GetBinContent(i+1)<<")"<<endl;
        point_index += 1;
    }

    //------------------------------------------------------------------
    // this is used to constrain the quadrant
    // info_vec[5];
    horizontal_fit_angle_diff_inner -> SetParameter(0,0);
    angle_diff_inner_phi_peak_profile_graph -> Fit(horizontal_fit_angle_diff_inner,"NQ","",0,360);
    //------------------------------------------------------------------

    angle_diff_inner_phi_peak_profile_graph -> Set(0);
    DCA_distance_outer_phi_peak_profile_graph -> Set(0);

    //------------------------------------------------------------------------
    // all others are not used for XY vertex calculation. for QA 
    //------------------------------------------------------------------------
     
      // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
      delete angle_diff_new_bkg_remove;
      angle_diff_new_bkg_remove = (TH1F*)angle_diff_new -> Clone("angle_diff_new_bkg_remove");
      angle_diff_new_bkg_remove -> SetLineColor(2);
      Hist_1D_bkg_remove(angle_diff_new_bkg_remove, 1.5);

    if(m_enable_qa){
      // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
      delete DCA_distance_outer_phi_peak;
      DCA_distance_outer_phi_peak = (TH2F*)DCA_distance_outer_phi -> Clone("DCA_distance_outer_phi_peak");
      TH2F_threshold(DCA_distance_outer_phi_peak, 0.5); // todo : the background cut can be modified, the ratio 0.5
      DCA_distance_outer_phi_peak_profile =  DCA_distance_outer_phi_peak->ProfileX("DCA_distance_outer_phi_peak_profile");
      point_index = 0;
      hist_column_content = SumTH2FColumnContent(DCA_distance_outer_phi_peak);
      for (int i = 0; i < DCA_distance_outer_phi_peak_profile->GetNbinsX(); i++){
          if (hist_column_content[i] < 5){continue;} // note : in order to remove some remaining background

          DCA_distance_outer_phi_peak_profile_graph -> SetPoint(point_index, 
                                                                DCA_distance_outer_phi_peak_profile->GetBinCenter(i+1), 
                                                                DCA_distance_outer_phi_peak_profile->GetBinContent(i+1));
          // cout<<"("<<DCA_distance_outer_phi_peak_profile->GetBinCenter(i+1)<<", "<< DCA_distance_outer_phi_peak_profile->GetBinContent(i+1)<<")"<<endl;
          point_index += 1;
      }
         
      horizontal_fit_outer -> SetParameter(0,0);
      // todo : the fit range of the gaussian fit can be modified here 
      DCA_distance_outer_phi_peak_profile_graph -> Fit(horizontal_fit_outer,"NQ","",0,360);

      // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
      delete angle_diff_outer_phi_peak;
      angle_diff_outer_phi_peak = (TH2F*)angle_diff_outer_phi -> Clone("angle_diff_outer_phi_peak");
      TH2F_threshold_advanced_2(angle_diff_outer_phi_peak, 0.5); // todo : threshold ratio can be modified here
      hist_column_content = SumTH2FColumnContent(angle_diff_outer_phi_peak);
      angle_diff_outer_phi_peak_profile =  angle_diff_outer_phi_peak->ProfileX("angle_diff_outer_phi_peak_profile");
      point_index = 0;
      for (int i = 0; i < angle_diff_outer_phi_peak_profile->GetNbinsX(); i++){
          if (hist_column_content[i] < 5){continue;} // note : in order to remove some remaining background
          
          angle_diff_outer_phi_peak_profile_graph -> SetPoint(point_index, 
                                                              angle_diff_outer_phi_peak_profile->GetBinCenter(i+1), 
                                                              angle_diff_outer_phi_peak_profile->GetBinContent(i+1));
          // cout<<"("<<angle_diff_outer_phi_peak_profile->GetBinCenter(i+1)<<", "<< angle_diff_outer_phi_peak_profile->GetBinContent(i+1)<<")"<<endl;
          point_index += 1;
      }

      horizontal_fit_angle_diff_outer -> SetParameter(0,0);
      angle_diff_outer_phi_peak_profile_graph -> Fit(horizontal_fit_angle_diff_outer,"NQ","",0,360);
   
      // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------

      angle_diff_outer_phi_peak_profile_graph -> Set(0);
      DCA_distance_inner_phi_peak_profile_graph -> Set(0);

    }
    
    if (m_enable_qa && print_message_opt == true) {
        cout<<"circle radius : "<<abs(gaus_fit->GetParameter(0) + gaus_fit->GetParameter(3))
            <<" possible correction angle : "<<gaus_fit->GetParameter(1)<<endl;
    }
}

//void INTTXYvtx::PrintPlotsVTXxy(string sub_out_folder_name)
void INTTXYvtx::PrintPlotsVTXxy()
{
    if(m_enable_drawhist)
    {
      string s_inttlabel = Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str());
      string s_pdfname = out_folder_directory + "/" + m_quad_pdfname;
      cout<<s_pdfname<<endl;
      
      // note : -----------------------------------------------------------------------------------------
      DCA_distance_inner_phi -> Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      //c1 -> Print(Form("%s/DCA_distance_inner_phi.pdf", sub_out_folder_name.c_str()));
      //c1 -> Print(Form("%s.pdf(", sub_out_folder_name.c_str()));
      c1 -> Print(s_pdfname.c_str());
      c1 -> Clear(); 

      // note : -----------------------------------------------------------------------------------------
      angle_diff_inner_phi -> Draw("colz0");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      //c1 -> Print(Form("%s/angle_diff_inner_phi.pdf", sub_out_folder_name.c_str()));
      //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
      c1 -> Print(s_pdfname.c_str());
      c1 -> Clear();

      // note : -----------------------------------------------------------------------------------------
      DCA_distance_inner_phi_peak -> SetStats(0);
      DCA_distance_inner_phi_peak -> GetXaxis() -> SetTitle("Inner phi [degree]");
      DCA_distance_inner_phi_peak -> GetYaxis() -> SetTitle("DCA [mm]");
      DCA_distance_inner_phi_peak -> Draw("colz0");
      DCA_distance_inner_phi_peak_profile -> Draw("same");
      // cos_fit -> Draw("l same");
      // gaus_fit -> Draw("l same");
      horizontal_fit_inner -> Draw("l same");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      draw_text -> DrawLatex(0.25, 0.84, Form("#color[2]{Assumed vertex: %.3f mm, %.3f mm}", current_vtxX, current_vtxY));
      draw_text -> DrawLatex(0.25, 0.80, Form("#color[2]{Pol0 fit chi2/NDF: %.3f, fit error: %.3f}", 
                                                          horizontal_fit_inner -> GetChisquare() / double(horizontal_fit_inner -> GetNDF()), 
                                                          horizontal_fit_inner->GetParError(0)));
      //c1 -> Print(Form("%s/DCA_distance_inner_phi_peak.pdf", sub_out_folder_name.c_str()));
      //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
      c1 -> Print(s_pdfname.c_str());
      c1 -> Clear(); 
 
      // note : -----------------------------------------------------------------------------------------
      angle_diff_inner_phi_peak -> SetStats(0);
      angle_diff_inner_phi_peak -> GetXaxis() -> SetTitle("Inner phi [degree]");
      angle_diff_inner_phi_peak -> GetYaxis() -> SetTitle("Inner - Outer [degree]");
      angle_diff_inner_phi_peak -> Draw("colz0");
      angle_diff_inner_phi_peak_profile -> Draw("same");
      horizontal_fit_angle_diff_inner -> Draw("l same");
      ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
      draw_text -> DrawLatex(0.25, 0.84, Form("#color[2]{Assumed vertex: %.3f mm, %.3f mm}", current_vtxX, current_vtxY));
      //c1 -> Print(Form("%s/angle_diff_inner_phi_peak.pdf", sub_out_folder_name.c_str()));
      //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
      c1 -> Print(s_pdfname.c_str());
      c1 -> Clear(); 


      //----------------------
      if(m_enable_qa){
        // note : -----------------------------------------------------------------------------------------
        DCA_distance_outer_phi -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        //c1 -> Print(Form("%s/DCA_distance_outer_phi.pdf", sub_out_folder_name.c_str()));
        //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
        c1 -> Print(s_pdfname.c_str());
        c1 -> Clear(); 
        
        // note : -----------------------------------------------------------------------------------------
        angle_diff_outer_phi -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        //c1 -> Print(Form("%s/angle_diff_outer_phi.pdf", sub_out_folder_name.c_str()));
        //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
        c1 -> Print(s_pdfname.c_str());
        c1 -> Clear();

        // note : -----------------------------------------------------------------------------------------
        angle_diff -> SetMinimum(0);
        angle_diff -> Draw("hist");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        //c1 -> Print(Form("%s/angle_diff.pdf", sub_out_folder_name.c_str()));
        //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
        c1 -> Print(s_pdfname.c_str());
        c1 -> Clear();
      
        // note : -----------------------------------------------------------------------------------------
        angle_diff_new -> SetMinimum(0);
        angle_diff_new -> Draw("hist");
        angle_diff_new_bkg_remove -> Draw("hist same");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        draw_text -> DrawLatex(0.4, 0.80, Form("#color[2]{Dist. StdDev: %.4f}", angle_diff_new_bkg_remove->GetStdDev()));
        //c1 -> Print(Form("%s/angle_diff_new.pdf", sub_out_folder_name.c_str()));
        //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
        c1 -> Print(s_pdfname.c_str());
        c1 -> Clear();
        
        // note : -----------------------------------------------------------------------------------------
        DCA_distance -> Draw("hist");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        //c1 -> Print(Form("%s/DCA_distance.pdf", sub_out_folder_name.c_str()));
        //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
        c1 -> Print(s_pdfname.c_str());
        c1 -> Clear(); 

        // note : -----------------------------------------------------------------------------------------
        angle_correlation -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        //c1 -> Print(Form("%s/angle_correlation.pdf", sub_out_folder_name.c_str()));
        //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
        c1 -> Print(s_pdfname.c_str());
        c1 -> Clear();

        // note : -----------------------------------------------------------------------------------------
        angle_diff_DCA_dist -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        //c1 -> Print(Form("%s/angle_diff_DCA_dist.pdf", sub_out_folder_name.c_str()));
        //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
        c1 -> Print(s_pdfname.c_str());
        c1 -> Clear();

        // note : -----------------------------------------------------------------------------------------
        DCA_point -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        //c1 -> Print(Form("%s/DCA_point.pdf", sub_out_folder_name.c_str()));
        //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
        c1 -> Print(s_pdfname.c_str());
        c1 -> Clear();

        // note : -----------------------------------------------------------------------------------------
        DCA_distance_inner_X -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        //c1 -> Print(Form("%s/DCA_distance_inner_X.pdf", sub_out_folder_name.c_str()));
        //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
        c1 -> Print(s_pdfname.c_str());
        c1 -> Clear();

        // note : -----------------------------------------------------------------------------------------
        DCA_distance_inner_Y -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        //c1 -> Print(Form("%s/DCA_distance_inner_Y.pdf", sub_out_folder_name.c_str()));
        //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
        c1 -> Print(s_pdfname.c_str());
        c1 -> Clear();

        // note : -----------------------------------------------------------------------------------------
        DCA_distance_outer_X -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        //c1 -> Print(Form("%s/DCA_distance_outer_X.pdf", sub_out_folder_name.c_str()));
        //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
        c1 -> Print(s_pdfname.c_str());
        c1 -> Clear();

        // note : -----------------------------------------------------------------------------------------
        DCA_distance_outer_Y -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        //c1 -> Print(Form("%s/DCA_distance_outer_Y.pdf", sub_out_folder_name.c_str()));
        c1 -> Print(s_pdfname.c_str());
        c1 -> Clear();

        // note : -----------------------------------------------------------------------------------------
        DCA_distance_outer_phi_peak -> SetStats(0);
        DCA_distance_outer_phi_peak -> GetXaxis() -> SetTitle("Outer phi [degree]");
        DCA_distance_outer_phi_peak -> GetYaxis() -> SetTitle("DCA [mm]");
        DCA_distance_outer_phi_peak -> Draw("colz0");
        DCA_distance_outer_phi_peak_profile -> Draw("same");
        horizontal_fit_outer -> Draw("l same");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        draw_text -> DrawLatex(0.25, 0.80, Form("#color[2]{Assumed vertex: %.3f mm, %.3f mm}", current_vtxX, current_vtxY));
        //c1 -> Print(Form("%s/DCA_distance_outer_phi_peak.pdf", sub_out_folder_name.c_str()));
        //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
        c1 -> Print(s_pdfname.c_str());
        c1 -> Clear(); 
      
        // note : -----------------------------------------------------------------------------------------
        angle_diff_outer_phi_peak -> SetStats(0);
        angle_diff_outer_phi_peak -> GetXaxis() -> SetTitle("Outer phi [degree]");
        angle_diff_outer_phi_peak -> GetYaxis() -> SetTitle("Inner - Outer [degree]");
        angle_diff_outer_phi_peak -> Draw("colz0");
        angle_diff_outer_phi_peak_profile -> Draw("same");
        horizontal_fit_angle_diff_outer -> Draw("l same");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("%s, peak : %f", s_inttlabel.c_str(), peek));
        draw_text -> DrawLatex(0.25, 0.84, Form("#color[2]{Assumed vertex: %.3f mm, %.3f mm}", current_vtxX, current_vtxY));
        //c1 -> Print(Form("%s/angle_diff_outer_phi_peak.pdf", sub_out_folder_name.c_str()));
        //c1 -> Print(Form("%s.pdf", sub_out_folder_name.c_str()));
        c1 -> Print(s_pdfname.c_str());
        c1 -> Clear(); 
      }

    }
}

void INTTXYvtx::ClearHist(int /*print_option*/)
{
    // clear histograms for quadrant method
    DCA_distance_inner_phi -> Reset("ICESM");
    angle_diff_inner_phi -> Reset("ICESM");

    if(m_enable_qa){
      DCA_distance_outer_phi -> Reset("ICESM");
      angle_diff_outer_phi -> Reset("ICESM");

      angle_correlation -> Reset("ICESM");
      angle_diff_DCA_dist -> Reset("ICESM");
      angle_diff -> Reset("ICESM");    
      DCA_point -> Reset("ICESM");
      DCA_distance -> Reset("ICESM");
      DCA_distance_inner_X -> Reset("ICESM");
      DCA_distance_inner_Y -> Reset("ICESM");
      DCA_distance_outer_X -> Reset("ICESM");
      DCA_distance_outer_Y -> Reset("ICESM");

      
      angle_diff_new -> Reset("ICESM");
      //
      // DCA_distance_inner_phi_peak -> Reset("ICESM");
      // DCA_distance_outer_phi_peak -> Reset("ICESM");
      // angle_diff_inner_phi_peak -> Reset("ICESM");
      // angle_diff_outer_phi_peak -> Reset("ICESM");
      // angle_diff_new_bkg_remove -> Reset("ICESM");
    }

    delete angle_diff_new_bkg_remove;   angle_diff_new_bkg_remove=nullptr;
    delete angle_diff_inner_phi_peak;   angle_diff_inner_phi_peak=nullptr;
    delete angle_diff_outer_phi_peak;   angle_diff_outer_phi_peak=nullptr;
    delete DCA_distance_inner_phi_peak; DCA_distance_inner_phi_peak=nullptr;
    delete DCA_distance_outer_phi_peak; DCA_distance_outer_phi_peak=nullptr;
}

void INTTXYvtx::EndRun()
{
    if(!m_initialized) {
       cout<<"INTTXYvtx is not initialized, abort in EndRun"<<endl;
       exit(1);
    }

  //------------------------
  // write histograms
    if(m_savehist)
    {
      if ( !std::filesystem::exists(out_folder_directory.c_str()) ) {
        system(Form("mkdir %s",out_folder_directory.c_str()));
      }

      TFile* file_out = new TFile(Form("%s/run_XY_histo.root",out_folder_directory.c_str()),"RECREATE");

      for(auto& itr: m_v_hist){
        itr->Write();
      }

      if(xy_hist      !=nullptr) xy_hist->Write();
      if(xy_hist_bkgrm!=nullptr) xy_hist_bkgrm->Write();
  	

      file_out -> Close();
      delete file_out;
    }

//ProcessEvt QA histos
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

void INTTXYvtx::TH2F_threshold(TH2F * hist, double threshold)
{
    double max_cut = hist -> GetMaximum() * threshold;

    for (int xi = 0; xi < hist -> GetNbinsX(); xi++){
        for(int yi = 0; yi < hist -> GetNbinsY(); yi++){
            if (hist -> GetBinContent(xi + 1, yi +1) < max_cut){ hist -> SetBinContent(xi + 1, yi +1, 0); }
        }
    }
}

void INTTXYvtx::TH2F_threshold_advanced_2(TH2F * hist, double threshold)
{
    // note : this function is to remove the background of the 2D histogram
    // note : but the threshold is given by average of the contents of the top "chosen_bin" bins and timing the threshold
    double max_cut = 0;
    int chosen_bin = 7;

    vector<float> all_bin_content_vec{}; //--all_bin_content_vec.clear();
    for (int xi = 0; xi < hist -> GetNbinsX(); xi++){
        for(int yi = 0; yi < hist -> GetNbinsY(); yi++){
            all_bin_content_vec.push_back(hist -> GetBinContent(xi + 1, yi +1));
        }
    }
    vector<unsigned long> ind(all_bin_content_vec.size(),0);
    TMath::Sort(all_bin_content_vec.size(), &all_bin_content_vec[0], &ind[0]);
    for (int i = 0; i < chosen_bin; i++) {
       max_cut += all_bin_content_vec[ind[i]]; 
       /*cout<<"test : "<<all_bin_content_vec[ind[i]]<<endl;*/
    }

    max_cut = (max_cut / double(chosen_bin)) * threshold;

    for (int xi = 0; xi < hist -> GetNbinsX(); xi++){
        for(int yi = 0; yi < hist -> GetNbinsY(); yi++){
            if (hist -> GetBinContent(xi + 1, yi +1) < max_cut){ hist -> SetBinContent(xi + 1, yi +1, 0); }
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

        // cout<<"slope : y="<<a<<"x+"<<b<<endl;
        
        // Calculate the closest distance from (target_x, target_y) to the line y = ax + b
        double closest_distance = std::abs(a * target_x - target_y + b) / std::sqrt(a * a + 1);

        // Calculate the coordinates of the closest point (Xc, Yc) on the line y = ax + b
        double Xc = (target_x + a * target_y - a * b) / (a * a + 1);
        double Yc = a * Xc + b;

        return { closest_distance, Xc, Yc };
    }
    else 
    {
        double closest_distance = std::abs(x1 - target_x);
        double Xc = x1;
        double Yc = target_y;

        return { closest_distance, Xc, Yc };
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
    
    // cout<<" crossProduct : "<<crossProduct<<endl;

    // Calculate the magnitudes of vector_1 and vector_2
    double magnitude1 = std::sqrt(vector1X * vector1X + vector1Y * vector1Y);
    double magnitude2 = std::sqrt(vector2X * vector2X + vector2Y * vector2Y);

    // Calculate the angle in radians using the inverse tangent of the cross product and dot product
    //--double dotProduct = vector1X * vector2X + vector1Y * vector2Y;

    //--double angleInRadians = std::atan2(std::abs(crossProduct), dotProduct);
    // Convert the angle from radians to degrees and return it
    //--double angleInDegrees = angleInRadians * 180.0 / M_PI;
    
    double angleInRadians_new = std::asin( crossProduct/(magnitude1*magnitude2) );
  
    //--double angleInDegrees_new = angleInRadians_new * 180.0 / M_PI;
    
    // cout<<"angle : "<<angleInDegrees_new<<endl;

    double DCA_value = sin(angleInRadians_new) * magnitude2; // DCA_value insread of DCA_distance

    return DCA_value;
}



void INTTXYvtx::Hist_1D_bkg_remove(TH1F * hist_in, double factor)
{   
    // todo : N bins considered to be used in the background quantification
    vector<double> Nbin_content_vec{};
    for (int i = hist_in -> GetNbinsX() - 5; i < hist_in -> GetNbinsX(); i++) { 
       Nbin_content_vec.push_back(hist_in -> GetBinContent(i+1));
    }

    double bkg_level = accumulate( Nbin_content_vec.begin(), Nbin_content_vec.end(), 0.0 ) / Nbin_content_vec.size();
    // cout<<"test, bkg cut : "<<bkg_level * factor<<endl;

    for (int i = 0; i < hist_in -> GetNbinsX(); i++){
        // note : the background rejection is here : bkg_level * 1.5 for the time being
        double bin_content = ( hist_in -> GetBinContent(i+1) <= bkg_level * factor) 
                           ? 0. 
                           : ( hist_in -> GetBinContent(i+1) - bkg_level * factor );

        hist_in -> SetBinContent(i+1, bin_content);
    }
}

void INTTXYvtx::DrawTGraphErrors(
             vector<double> x_vec, 
             vector<double> y_vec, 
             vector<double> xE_vec, 
             vector<double> yE_vec, 
             string         output_directory, 
             vector<string> plot_name)
{
    if(m_enable_drawhist) 
    {
        c1 -> cd();

        TGraphErrors * g = new TGraphErrors(x_vec.size(), &x_vec[0], &y_vec[0], &xE_vec[0], &yE_vec[0]);
        g -> SetMarkerStyle(20);
        g -> SetMarkerSize(1.5);
        g -> SetMarkerColor(1);
        g -> SetLineColor(1);
        g -> GetXaxis() -> SetTitle(plot_name[1].c_str());
        g -> GetYaxis() -> SetTitle(plot_name[2].c_str());
        g -> SetTitle(plot_name[0].c_str());
        if (plot_name.size() == 4){g -> Draw(plot_name[3].c_str());}
        else {g -> Draw("AP");}

        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        //c1 -> Print(Form("%s/%s.pdf", output_directory.c_str(), plot_name[0].c_str()));
        c1 -> Print(Form("%s/%s", output_directory.c_str(), plot_name[4].c_str()));
        c1 -> Clear();
        
        delete g;
    }
}

void INTTXYvtx::Draw2TGraph(
              vector<double> x1_vec, 
              vector<double> y1_vec, 
              vector<double> x2_vec, 
              vector<double> y2_vec, 
              string         output_directory, 
              vector<string> plot_name)
{
    if(m_enable_drawhist) 
    {
        c1 -> cd();
        c1 -> SetLogy(1);

        TGraph * g1 = new TGraph(x1_vec.size(), &x1_vec[0], &y1_vec[0]);
        g1 -> SetMarkerStyle(5);
        g1 -> SetMarkerSize(1);
        g1 -> SetMarkerColor(1);
        g1 -> SetLineColor(1);
        g1 -> GetXaxis() -> SetTitle(plot_name[1].c_str());
        g1 -> GetYaxis() -> SetTitle(plot_name[2].c_str());
        g1 -> GetXaxis() -> SetNdivisions(505);
        g1 -> GetXaxis() -> SetLimits(-1, x1_vec[x1_vec.size()-1] + 1);
        g1 -> SetTitle(plot_name[0].c_str());
        g1 -> Draw("AP");

        TGraph * g2 = new TGraph(x2_vec.size(), &x2_vec[0], &y2_vec[0]);
        g2 -> SetMarkerStyle(5);
        g2 -> SetMarkerSize(1);
        g2 -> SetMarkerColor(2);
        g2 -> SetLineColor(2);
        g2 -> Draw("PL same");

        TLegend * legend = new TLegend(0.4,0.75,0.65,0.9);
        // legend -> SetMargin(0);
        legend->SetTextSize(0.03);
        legend -> AddEntry(g1, "Tested vertex candidates", "p");
        legend -> AddEntry(g2, "Better performed candidates", "p");
        legend -> Draw("same");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        //c1 -> Print(Form("%s/%s.pdf", output_directory.c_str(), plot_name[0].c_str()));
        c1 -> Print(Form("%s/%s", output_directory.c_str(), plot_name[4].c_str()));
        c1 -> Clear();
        c1 -> SetLogy(0);

        delete g1;
        delete g2;
        delete legend;
    }
}

vector<double> INTTXYvtx::SumTH2FColumnContent(TH2F * hist_in)
{
    vector<double> sum_vec; sum_vec.clear();
    for (int i = 0; i < hist_in -> GetNbinsX(); i++){
        double sum = 0;
        for (int j = 0; j < hist_in -> GetNbinsY(); j++){
            sum += hist_in -> GetBinContent(i+1, j+1);
        }
        sum_vec.push_back(sum);
    }
    return sum_vec;
}


vector<pair<double,double>> INTTXYvtx::Get4vtx(pair<double,double> origin, double length)
{
    vector<pair<double,double>> unit_vtx = {{1,1},{-1,1},{-1,-1},{1,-1}};
    vector<pair<double,double>> vec_out{};//-- vec_out.clear();

    for (pair i1 : unit_vtx)
    {
        vec_out.push_back({i1.first * length + origin.first, i1.second * length + origin.second});
    }

    return vec_out;
}

void INTTXYvtx::TH2F_FakeClone(TH2F*hist_in, TH2F*hist_out)
{
    if (hist_in -> GetNbinsX() != hist_out -> GetNbinsX() || 
        hist_in -> GetNbinsY() != hist_out -> GetNbinsY())
    {
        cout<<"In INTTXYvtx::TH2F_FakeClone, the input and output histogram have different binning!"<<endl;
        return;
    }

    for (int i = 0; i < hist_in -> GetNbinsX(); i++){
        for (int j = 0; j < hist_in -> GetNbinsY(); j++){
            hist_out -> SetBinContent(i+1, j+1, hist_in -> GetBinContent(i+1, j+1));
        }
    }
}

void INTTXYvtx::TH1F_FakeClone(TH1F*hist_in, TH1F*hist_out)
{
    if (hist_in -> GetNbinsX() != hist_out -> GetNbinsX())
    {
        cout<<"In INTTXYvtx::TH1F_FakeClone, the input and output histogram have different binning!"<<endl;
        return;
    }

    for (int i = 0; i < hist_in -> GetNbinsX(); i++){
        hist_out -> SetBinContent(i+1, hist_in -> GetBinContent(i+1));
    }
}

void INTTXYvtx::TH2FSampleLineFill(
        TH2F*                    hist_in, 
        double                   segmentation, 
        std::pair<double,double> inner_clu, 
        std::pair<double,double> outer_clu)
{
    if(!m_initialized) {
       cout<<"INTTXYvtx is not initialized, abort in MacroVTXSquare"<<endl;
       exit(1);
    }

    double x_min = hist_in -> GetXaxis() -> GetXmin();
    double x_max = hist_in -> GetXaxis() -> GetXmax();
    double y_min = hist_in -> GetYaxis() -> GetXmin();
    double y_max = hist_in -> GetYaxis() -> GetXmax();

    double seg_x, seg_y;
    double angle;
    int n_seg = 0;

    while (true)
    {
        angle = atan2(inner_clu.second-outer_clu.second, inner_clu.first-outer_clu.first);
        seg_x = (n_seg * segmentation) * cos(angle) + outer_clu.first; // note : atan2(y,x), point.first is the radius
        seg_y = (n_seg * segmentation) * sin(angle) + outer_clu.second;
        
        if ( (seg_x > x_min && seg_x < x_max && seg_y > y_min && seg_y < y_max) != true ) {break;}

        hist_in -> Fill(seg_x, seg_y);
        n_seg += 1;
    }

    n_seg = 1;
    while (true)
    {
        angle = atan2(inner_clu.second-outer_clu.second, inner_clu.first-outer_clu.first);
        seg_x = (-1 * n_seg * segmentation) * cos(angle) + outer_clu.first; // note : atan2(y,x), point.first is the radius
        seg_y = (-1 * n_seg * segmentation) * sin(angle) + outer_clu.second;
        
        if ( (seg_x > x_min && seg_x < x_max && seg_y > y_min && seg_y < y_max) != true ) {break;}
        hist_in -> Fill(seg_x, seg_y);
        n_seg += 1;
    }
}

vector<pair<double,double>> 
INTTXYvtx::FillLine_FindVertex(
             pair<double,double> window_center, 
             double              segmentation, 
             double              window_width, 
             int                 N_bins
            )
{
    bool draw_plot = m_enable_drawhist;
  
    if(cluster_pair_vec.size()==0){ // minimum tracklet cut. should be tuned
      return {beam_origin,
              {0, 0}, 
              {0, 0}
             };   
    }

    delete xy_hist; //if xy_hist is nullptr, nothing happen;
    xy_hist = new TH2F("xy_hist","xy_hist", 
                          N_bins, -1 * window_width / 2. + window_center.first, 
                                       window_width / 2. + window_center.first, 
                          N_bins, -1 * window_width / 2. + window_center.second, 
                                       window_width / 2. + window_center.second);
    xy_hist -> SetStats(0);
    xy_hist -> GetXaxis() -> SetTitle("X axis [mm]");
    xy_hist -> GetYaxis() -> SetTitle("Y axis [mm]");
    xy_hist -> GetXaxis() -> SetNdivisions(505);

    delete xy_hist_bkgrm;
    xy_hist_bkgrm = new TH2F("xy_hist_bkgrm","xy_hist_bkgrm", 
                                      N_bins, -1 * window_width / 2. + window_center.first, 
                                                   window_width / 2. + window_center.first, 
                                      N_bins, -1 * window_width / 2. + window_center.second, 
                                                   window_width / 2. + window_center.second);
    xy_hist_bkgrm -> SetStats(0);
    xy_hist_bkgrm -> GetXaxis() -> SetTitle("X axis [mm]");
    xy_hist_bkgrm -> GetYaxis() -> SetTitle("Y axis [mm]");
    xy_hist_bkgrm -> GetXaxis() -> SetNdivisions(505);

    // cout<<"test test size and bin of the hist xy_hist : "<<xy_hist -> GetNbinsX()<<" "<<xy_hist -> GetNbinsY()<<endl;
    // cout<<"test test bin width of the hist xy_hist : "<<xy_hist -> GetXaxis() -> GetBinWidth(1)<<" "<<xy_hist -> GetYaxis() -> GetBinWidth(1)<<endl;
    // cout<<"draw_plot status : "<<draw_plot<<endl;

    
    for (unsigned int i = 0; i < cluster_pair_vec.size(); i++)
    {
        vector<double> DCA_info_vec = calculateDistanceAndClosestPoint(
            cluster_pair_vec[i].first.x,  cluster_pair_vec[i].first.y,
            cluster_pair_vec[i].second.x, cluster_pair_vec[i].second.y,
            window_center.first, window_center.second
        );

        double DCA_sign = calculateAngleBetweenVectors(
            cluster_pair_vec[i].second.x, cluster_pair_vec[i].second.y,
            cluster_pair_vec[i].first.x,  cluster_pair_vec[i].first.y,
            window_center.first, window_center.second
        );

        if (DCA_info_vec[0] != fabs(DCA_sign) && fabs( DCA_info_vec[0] - fabs(DCA_sign) ) > 0.1){
            cout<<"different DCA : "<<DCA_info_vec[0]<<" "<<DCA_sign<<" diff : "<<DCA_info_vec[0] - fabs(DCA_sign)<<endl;
        }

        Clus_InnerPhi_Offset = (cluster_pair_vec[i].first.y - window_center.second < 0) 
                                ? atan2(cluster_pair_vec[i].first.y - window_center.second, 
                                        cluster_pair_vec[i].first.x - window_center.first) * (180./TMath::Pi()) + 360 
                                : atan2(cluster_pair_vec[i].first.y - window_center.second, 
                                        cluster_pair_vec[i].first.x - window_center.first) * (180./TMath::Pi());
        Clus_OuterPhi_Offset = (cluster_pair_vec[i].second.y - window_center.second < 0) 
                                ? atan2(cluster_pair_vec[i].second.y - window_center.second, 
                                        cluster_pair_vec[i].second.x - window_center.first) * (180./TMath::Pi()) + 360 
                                : atan2(cluster_pair_vec[i].second.y - window_center.second, 
                                        cluster_pair_vec[i].second.x - window_center.first) * (180./TMath::Pi());

        if (fabs(Clus_InnerPhi_Offset - Clus_OuterPhi_Offset) < 5)
        {
            TH2FSampleLineFill(xy_hist, 
                               segmentation, 
                               {cluster_pair_vec[i].first.x, cluster_pair_vec[i].first.y}, 
                               {DCA_info_vec[1], DCA_info_vec[2]}
                              );
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

    double reco_vtx_x = xy_hist_bkgrm->GetMean(1);// note : the TH2F calculate the GetMean based on the bin center, no need to apply additional offset
    double reco_vtx_y = xy_hist_bkgrm->GetMean(2);// note : the TH2F calculate the GetMean based on the bin center, no need to apply additional offset

    // cout<<"test : in the line filled, the process is almost done"<<endl;

    if (draw_plot)
    {
        TGraph * reco_vertex_gr = new TGraph(); 
        reco_vertex_gr -> SetMarkerStyle(50);
        reco_vertex_gr -> SetMarkerColor(2);
        reco_vertex_gr -> SetMarkerSize(1);
        reco_vertex_gr -> SetPoint(reco_vertex_gr -> GetN(), reco_vtx_x, reco_vtx_y);


        string s_inttlabel = Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str());
        // note : -----------------------------------------------------------------------------------------
        xy_hist -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        //c1 -> Print(Form("%s/Run_xy_hist.pdf",out_folder_directory.c_str()));
        c1 -> Print(Form("%s/linefill_qa.pdf(",out_folder_directory.c_str()));
        c1 -> Clear();

        // note : -----------------------------------------------------------------------------------------
        xy_hist_bkgrm -> Draw("colz0");
        draw_text -> DrawLatex(0.21, 0.71+0.13, Form("Vertex of the Run: %.4f mm, %.4f mm", 
                                                     reco_vtx_x, reco_vtx_y));
        draw_text -> DrawLatex(0.21, 0.67+0.13, Form("Vertex error: %.4f mm, %.4f mm", 
                                                     xy_hist_bkgrm->GetMeanError(1), xy_hist_bkgrm->GetMeanError(2)));
        reco_vertex_gr -> Draw("p same");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, s_inttlabel.c_str());
        //c1 -> Print(Form("%s/Run_xy_hist_bkgrm.pdf",out_folder_directory.c_str()));
        c1 -> Print(Form("%s/linefill_qa.pdf)",out_folder_directory.c_str()));
        c1 -> Clear();

        // cout<<"test : hello, can you see me ?"<<endl;
    }

    return {{reco_vtx_x,reco_vtx_y},
            {xy_hist_bkgrm->GetMeanError(1), xy_hist_bkgrm->GetMeanError(2)}, 
            {xy_hist_bkgrm->GetStdDev(1),    xy_hist_bkgrm->GetStdDev(2)}
           };   
}

#endif
