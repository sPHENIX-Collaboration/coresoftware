#ifndef QA_TPCCALIBRATIONQA_H
#define QA_TPCCALIBRATIONQA_H

#include <string>

// Analyze existing QA histogram files, without rerunning event reconstruction.
// Each run() writes a new output directory. Use separate processes for concurrent jobs.
class TpcCalibrationQA
{
 public:
  void set_pt_range(double minimum, double maximum)
  {
    m_ptMin = minimum;
    m_ptMax = maximum;
  }
  void set_min_entries(int value) { m_minEntries = value; }
  void set_phi_group(int value) { m_phiGroup = value; }
  void set_eta_group(int value) { m_etaGroup = value; }
  void set_histogram_prefix(const std::string &value) { m_prefix = value; }
  void set_save_slices(bool value) { m_saveSlices = value; }
  void set_ratio_floor(double value) { m_ratioFloor = value; }
  void set_phi_qoverpt(bool value) { m_phiQpt = value; }
  void set_conditions_json(const std::string &value) { m_conditions = value; }
  void set_raw_time_shift(double value) { m_rawTimeShift = value; }
  void set_make_plots(bool value) { m_makePlots = value; }
  void set_phgarfield_table(const std::string &value) { m_phgReferences = value; }
  void set_fit_background(bool value) { m_fitBackground = value; }

  // Manifest rows: label<TAB>file-list[<TAB>drift-file-list<TAB>raw-time-shift].
  int run(const std::string &manifest, const std::string &output,
          const std::string &reference = "") const;

  static void cache(const char *listing, const char *output,
                    const char *prefix = "h_TpcSeedsQA_");
  static void collect(const char *manifest, const char *references, const char *output,
                      bool maps = false, bool plots = false, double ratioFloor = .02);
  static void attach_conditions(const char *table, const char *output);

 private:
  double m_ptMin{.5};
  double m_ptMax{1.};
  int m_minEntries{100};
  int m_phiGroup{3};
  int m_etaGroup{4};
  std::string m_prefix{"h_TpcSeedsQA_"};
  bool m_saveSlices{false};
  double m_ratioFloor{.02};
  bool m_phiQpt{false};
  std::string m_conditions{"[]"};
  double m_rawTimeShift{40.};
  bool m_makePlots{true};
  std::string m_phgReferences;
  bool m_fitBackground{false};
};
#endif
