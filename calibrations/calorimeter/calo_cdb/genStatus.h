#pragma once

// -- c++ includes --
#include <memory>
#include <string>

// -- root includes --
#include <TProfile2D.h>

class GenStatus
{
 public:
  GenStatus();

  void process(const std::string &input, const std::string &output = "output");

 private:

  static void histToCaloCDBTree(const std::string &outputfile, const std::string &fieldName, Int_t icalo, TProfile2D* hist);
  void analyze(const std::string &output);

  // utils
  void setRunDataset(const std::string &input);
  Int_t readHists(const std::string &input);

  std::string m_run;
  std::string m_dataset;

  std::unique_ptr<TProfile2D> h_CaloValid_cemc_etaphi_badChi2;
  std::unique_ptr<TProfile2D> h_CaloValid_ihcal_etaphi_badChi2;
  std::unique_ptr<TProfile2D> h_CaloValid_ohcal_etaphi_badChi2;

  std::unique_ptr<TProfile2D> h_CaloValid_cemc_etaphi_time_raw;
  std::unique_ptr<TProfile2D> h_CaloValid_ihcal_etaphi_time_raw;
  std::unique_ptr<TProfile2D> h_CaloValid_ohcal_etaphi_time_raw;

  Int_t cemc_bins_eta;
  Int_t cemc_bins_phi;
  Int_t hcal_bins_eta;
  Int_t hcal_bins_phi;
};
