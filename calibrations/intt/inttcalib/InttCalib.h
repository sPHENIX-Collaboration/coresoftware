#ifndef INTTCALIB_H
#define INTTCALIB_H

#include <intt/InttFeeMap.h>
#include <intt/InttMap.h>
#include <intt/InttSurveyMap.h>

#include <fun4all/SubsysReco.h>

#include <RtypesCore.h>

#include <array>
#include <map>
#include <string>
#include <utility>

class PHCompositeNode;
class TF1;
class TH1D;
class TH2D;

class InttCalib : public SubsysReco
{
 public:
  InttCalib(std::string const& = "InttCalib");
  ~InttCalib() override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  /// Called at the end of each run.
  int EndRun(int const) override;

  void SetHotMapCdbFile(std::string const& file) { m_hotmap_cdb_file = file; }
  void SetHotMapPngFile(std::string const& file) { m_hotmap_png_file = file; }

  void SetBcoMapCdbFile(std::string const& file) { m_bcomap_cdb_file = file; }
  void SetBcoMapPngFile(std::string const& file) { m_bcomap_png_file = file; }
  void SetStreamingMode(bool mode) { m_streaming = mode; }
  void SetppMode(bool mode) { m_ppmode = mode; }
  void SetBcoMaximumEvent(int mext) {m_evts_bco = mext; }
  void SetRunNumber(int runnum) { m_run_num = runnum; }
  void SetDoFeebyFee(bool in) {m_do_fee = in;}
  int SaveHitrates();
  int LoadHitrates();

  /// For debugging
  void Debug();
 private:

  int ConfigureHotMap_v3();
  int MakeHotMapCdb_v3();
  int MakeHotMapPng_v3();

  
  int ConfigureHotMap_v2();
  int MakeHotMapCdb_v2();
  int MakeHotMapPng_v2();

  int ConfigureHotMap();
  int MakeHotMapCdb();
  int MakeHotMapPng();

  int ConfigureBcoMap();
  int MakeBcoMapCdb();
  int MakeBcoMapPng();

  int ConfigureHist(TH1D*&, TF1*&, std::map<double, int> const&, std::string const&, std::string const&);
  int ConfigureHist_v2(TH1D*&, TF1*&, std::map<double, int> const&, std::string const&, std::string const&);
  int ConfigureHist_v3(TH1D*&, TF1*&, double, std::map<double, int> const&, std::string const&, std::string const&);
  int adjust_hitrate(InttMap::Offline_s const&, double&) const;
  int GetIndex(InttMap::RawData_s const&, InttMap::Offline_s const&) const;
  int GetFeeIndex(InttMap::RawData_s const&, InttMap::Offline_s const&) const;
  void SetColdSigmaCut(double in) {m_NUM_SIGMA_COLD = in;}
  void SetHotSigmaCut(double in) {m_NUM_SIGMA_HOT = in;}
  // For Fee by Fee
  int ConfigureHotMap_fee();
  int MakeHotMapCdb_fee();
  int MakeHotMapROOT_fee();
  // For Fee by Fee
  
  std::pair<double, double> CalculateStandardDeviation(const std::vector<int>& data);
  Color_t GetFeeColor(int) const;

  int m_evts{0};
  int m_run_num{0};
  double m_bco_stdDev{0};  
  double m_bco_mean{0};  
  double m_mean[8] = {0.};
  double m_sigma[8] = {0.};
  double m_mean_fee[112] = {0.};
  double m_sigma_fee[112] = {0.};
  int m_evts_bco = 50000;
  // int static const m_MAX_INDEX = 32;
  int static const m_MAX_INDEX = 8;
  int static const m_MAX_LADDER = 112;
  double static constexpr m_NUM_CHANNELS = 8 * 14 * 26 * 128;
  double  m_NUM_SIGMA_HOT = 5.0;
  double  m_NUM_SIGMA_COLD = 3.0;

  double m_min_hitrate{0.0};
  double m_min_fraction{0.0};
  double m_max_hitrate{0.0};
  double m_max_fraction{0.0};

  std::string m_hotmap_cdb_file;
  std::string m_hotmap_png_file;
  std::string m_bcomap_cdb_file;
  std::string m_bcomap_png_file;

  InttFeeMap m_feemap;
  InttSurveyMap m_survey;
  Eigen::Vector3d m_vertex{0.0, 0.0, 0.0};
  // int m_hitmap[8][14][26][128][129]
  std::array<std::array<std::array<std::array<std::array<double, 129>, 128>, 26>, 14>, 8> m_hitmap{};
  std::array<std::array<std::array<int, 26>, 14>, 8> m_hitmap_half{};
  
// TH1D* m_hist[8][14]
  std::array<TH1D*, m_MAX_INDEX> m_hist{};
  std::array<TH1D*, m_MAX_INDEX> m_hist_half{};
  std::array<TF1*, m_MAX_INDEX> m_fit{};
  std::array<double, m_MAX_INDEX> m_min{};
  std::array<double, m_MAX_INDEX> m_max{};
  std::array<double, m_MAX_LADDER> m_min_fee{};
  std::array<double, m_MAX_LADDER> m_max_fee{};
  std::array<double, m_MAX_INDEX> m_half_min{};
  std::array<double, m_MAX_INDEX> m_half_max{};

  std::array<TH1D*, m_MAX_LADDER> m_hist_fee{};
  std::array<TF1*, m_MAX_LADDER> m_fit_fee{};
  
  std::map<double, double> m_hitrates;
  std::map<double, double> m_invcdf;

  std::map<InttMap::RawData_s, int[128], InttMap::RawDataWildcardComparator> m_bcorates;
  std::map<InttMap::RawData_s, int, InttMap::RawDataWildcardComparator> m_bcopeaks;

  bool m_do_nothing = false;
  bool m_streaming = false;
  bool m_ppmode = true;
  bool m_do_make_bco = true;
  bool m_do_fee = false;
};

#endif  // INTTCALIB_H
