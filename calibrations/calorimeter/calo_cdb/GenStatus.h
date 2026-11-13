#ifndef CALOCDB_GENSTATUS_H
#define CALOCDB_GENSTATUS_H

#include "geometry_constants.h"

// -- root includes --
#include <TProfile2D.h>

// -- c++ includes --
#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

class GenStatus
{
 public:
  explicit GenStatus(std::string calofittingPrefix = "HIST_DST_CALOFITTING");

  void process(const std::string &input, const std::string &output = "output");

 private:
  static void histToCaloCDBTree(const std::string &outputfile, const std::string &fieldName, int icalo, TProfile2D *hist);
  static void histToCaloCDBTree(const std::string &outputfile, const std::vector<std::pair<std::string, TProfile2D *>> &fields, int icalo);
  void analyze(const std::string &output);
  void generateHotMaps(const std::string &outputDir, const std::string &qaDir);

  // utils
  void setRunDataset(const std::string &input);
  int readHists(const std::string &input);

  std::string m_run;
  std::string m_dataset;
  std::string m_dataset_jetqa;
  std::string m_dataset_calofittingqa;
  std::string extractTag(const std::string& filename, const std::string& prefix);

  std::unique_ptr<TProfile2D> h_CaloValid_cemc_etaphi_badChi2;
  std::unique_ptr<TProfile2D> h_CaloValid_ihcal_etaphi_badChi2;
  std::unique_ptr<TProfile2D> h_CaloValid_ohcal_etaphi_badChi2;

  std::unique_ptr<TProfile2D> h_CaloValid_cemc_etaphi_time_raw;
  std::unique_ptr<TProfile2D> h_CaloValid_ihcal_etaphi_time_raw;
  std::unique_ptr<TProfile2D> h_CaloValid_ohcal_etaphi_time_raw;

  std::unique_ptr<TProfile2D> h_CaloFittingQA_cemc_etaphi_ZScrosscalib;
  std::unique_ptr<TProfile2D> h_CaloFittingQA_ihcal_etaphi_ZScrosscalib;
  std::unique_ptr<TProfile2D> h_CaloFittingQA_ohcal_etaphi_ZScrosscalib;

  // ZSCrossCalibT0: per-detector (cemc, ihcal, ohcal) x per-variant (minus5, main, plus5)
  static constexpr std::array<const char *, 3> m_zsT0_detectors{"cemc", "ihcal", "ohcal"};
  static constexpr std::array<const char *, 3> m_zsT0_variants{"minus5", "main", "plus5"};
  std::array<std::array<std::unique_ptr<TProfile2D>, 3>, 3> h_CaloFittingQA_etaphi_ZScrosscalibT0;

  // dsttype prefix of the calofitting QA histogram files (e.g. HIST_DST_CALOFITTING or HIST_CALOFITTINGQA)
  std::string m_calofittingPrefix;

  int cemc_bins_eta{CaloGeometry::CEMC_ETA_BINS};
  int cemc_bins_phi{CaloGeometry::CEMC_PHI_BINS};
  int hcal_bins_eta{CaloGeometry::HCAL_ETA_BINS};
  int hcal_bins_phi{CaloGeometry::HCAL_PHI_BINS};

  std::string m_CaloValid_list{"CaloValid.list"};
};

#endif
