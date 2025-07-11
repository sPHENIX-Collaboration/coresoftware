#ifndef CALOCDB_GENSTATUS_H
#define CALOCDB_GENSTATUS_H

#include "geometry_constants.h"

// -- root includes --
#include <TProfile2D.h>

// -- c++ includes --
#include <memory>
#include <string>

class GenStatus
{
 public:
  GenStatus() = default;

  void process(const std::string &input, const std::string &output = "output");

 private:
  static void histToCaloCDBTree(const std::string &outputfile, const std::string &fieldName, int icalo, TProfile2D *hist);
  void analyze(const std::string &output);

  // utils
  void setRunDataset(const std::string &input);
  int readHists(const std::string &input);

  std::string m_run;
  std::string m_dataset;

  std::unique_ptr<TProfile2D> h_CaloValid_cemc_etaphi_badChi2;
  std::unique_ptr<TProfile2D> h_CaloValid_ihcal_etaphi_badChi2;
  std::unique_ptr<TProfile2D> h_CaloValid_ohcal_etaphi_badChi2;

  std::unique_ptr<TProfile2D> h_CaloValid_cemc_etaphi_time_raw;
  std::unique_ptr<TProfile2D> h_CaloValid_ihcal_etaphi_time_raw;
  std::unique_ptr<TProfile2D> h_CaloValid_ohcal_etaphi_time_raw;

  int cemc_bins_eta{CaloGeometry::CEMC_ETA_BINS};
  int cemc_bins_phi{CaloGeometry::CEMC_PHI_BINS};
  int hcal_bins_eta{CaloGeometry::HCAL_ETA_BINS};
  int hcal_bins_phi{CaloGeometry::HCAL_PHI_BINS};
};

#endif
