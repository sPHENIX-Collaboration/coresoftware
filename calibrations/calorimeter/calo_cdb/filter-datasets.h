#pragma once

// -- c++ includes --
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <memory>

// -- ROOT includes --
#include <RtypesCore.h>

// -- sPHENIX includes --
#include <sphenixnpc/CDBUtils.h>

class FilterDatasets
{
 public:
  explicit FilterDatasets(Bool_t debug = false);

  void process(const std::string &input, const std::string &output = ".");

 private:

  void analyze(const std::string& input, const std::string &outputDir);
  void readRunInfo(const std::string &line);

  std::string getCalibration(const std::string &pl_type, uint64_t iov);
  Int_t setGlobalTag(const std::string &tagname);

  std::vector<std::pair<std::string, std::string>> m_runInfo;
  std::map<std::string, Int_t> m_ctr;

  std::vector<std::string> m_cdbName = {"CEMC_BadTowerMap"
                                      , "CEMC_meanTime", "HCALIN_meanTime", "HCALOUT_meanTime"
                                      , "CEMC_hotTowers_fracBadChi2", "HCALIN_hotTowers_fracBadChi2", "HCALOUT_hotTowers_fracBadChi2"};

  Bool_t m_debug;

  std::unique_ptr<CDBUtils> uti;
};
