#ifndef INTTCALIB_H
#define INTTCALIB_H

#include <intt/InttFeeMapv1.h>
#include <intt/InttMap.h>
#include <intt/InttSurveyMapv1.h>

#include <fun4all/SubsysReco.h>

#include <RtypesCore.h>

#include <array>
#include <map>
#include <string>

class PHCompositeNode;
class TF1;
class TH1D;

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

  int SaveHitrates();
  int LoadHitrates();

  /// For debugging
  void Debug();

 private:
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
  int adjust_hitrate(InttMap::Offline_s const&, double&) const;
  int GetIndex(InttMap::RawData_s const&, InttMap::Offline_s const&) const;

  Color_t GetFeeColor(int) const;

  int m_evts{0};
  int m_run_num{0};
  int m_evts_bco = 30000;

  // int static const m_MAX_INDEX = 32;
  int static const m_MAX_INDEX = 8;
  double static constexpr m_NUM_CHANNELS = 8 * 14 * 26 * 128;
  double static constexpr m_NUM_SIGMA = 3.0;

  double m_min_hitrate{0.0};
  double m_min_fraction{0.0};
  double m_max_hitrate{0.0};
  double m_max_fraction{0.0};

  std::string m_hotmap_cdb_file;
  std::string m_hotmap_png_file;
  std::string m_bcomap_cdb_file;
  std::string m_bcomap_png_file;

  InttFeeMapv1 m_feemap;
  InttSurveyMapv1 m_survey;
  Eigen::Vector3d m_vertex{0.0, 0.0, 0.0};

  // int m_hitmap[8][14][26][128][129]
  std::array<std::array<std::array<std::array<std::array<double, 129>, 128>, 26>, 14>, 8> m_hitmap{};

  // TH1D* m_hist[8][14]
  std::array<TH1D*, m_MAX_INDEX> m_hist{};
  std::array<TF1*, m_MAX_INDEX> m_fit{};
  std::array<double, m_MAX_INDEX> m_min{};
  std::array<double, m_MAX_INDEX> m_max{};

  std::map<double, double> m_hitrates;
  std::map<double, double> m_invcdf;

  std::map<InttMap::RawData_s, int[128], InttMap::RawDataWildcardComparator> m_bcorates;
  std::map<InttMap::RawData_s, int, InttMap::RawDataWildcardComparator> m_bcopeaks;

  bool m_do_nothing = false;
  bool m_do_make_bco = true;
};

#endif  // INTTCALIB_H
