#ifndef INTTCALIB_H
#define INTTCALIB_H

#include <intt/InttMap.h>
#include <intt/InttFeeMapv1.h>
#include <intt/InttSurveyMapv1.h>

#include <fun4all/SubsysReco.h>

#include <RtypesCore.h>

#include <array>
#include <map>
#include <string>

class PHCompositeNode;

class InttCalib : public SubsysReco
{
 public:
  InttCalib(std::string const& = "InttCalib");
  ~InttCalib() override;

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode*) override;

  /// Called at the end of each run.
  int EndRun(int const) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode*) override;

  /// Reset
  int Reset(PHCompositeNode*) override;

  void Print(const std::string &what = "ALL") const override;

  void SetHotMapCdbFile(std::string const& file) {m_hotmap_cdb_file = file;}
  void SetHotMapPngFile(std::string const& file) {m_hotmap_png_file = file;}

  void SetBcoMapCdbFile(std::string const& file) {m_bcomap_cdb_file = file;}
  void SetBcoMapPngFile(std::string const& file) {m_bcomap_png_file = file;}

  /// For debugging
  void Debug() const;

 private:
  int ConfigureHotMap();
  int MakeHotMapCdb();
  int MakeHotMapPng();

  int ConfigureBcoMap();
  int MakeBcoMapCdb();
  int MakeBcoMapPng();

  int adjust_hitrate(InttMap::Offline_s const&, double&) const;
  Color_t GetFeeColor(int) const;

  int m_evts{0};
  int m_run_num{0};

  double static constexpr m_NUM_CHANNELS = 8 * 14 * 26 * 128;

  double m_min_hitrate{0.0};
  double m_min_fraction{0.0};
  double m_max_hitrate{0.0};
  double m_max_fraction{0.0};

  std::string m_hotmap_cdb_file;
  std::string m_hotmap_png_file;
  std::string m_bcomap_cdb_file;
  std::string m_bcomap_png_file;

  InttFeeMapv1    m_feemap;
  InttSurveyMapv1 m_survey;
  Eigen::Vector3d m_vertex{0.0, 0.0, 0.0};

  std::array<std::array<std::array<std::array<std::array<int, 129>, 128>, 26>, 14>, 8> m_hitmap;
  // m_hitmap[8][14][26][128][129]

  std::map<double, double> m_hitrates;
  std::map<double, double> m_invcdf;

  std::map<InttMap::RawData_s, int[128], InttMap::RawDataWildcardComparator> m_bcorates;
  std::map<InttMap::RawData_s, int,      InttMap::RawDataWildcardComparator> m_bcopeaks;

  bool m_do_nothing = false;
};

#endif // INTTCALIB_H
