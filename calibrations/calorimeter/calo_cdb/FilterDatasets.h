#ifndef CALOCDB_FILTERDATASETS_H
#define CALOCDB_FILTERDATASETS_H

// -- c++ includes --
#include <map>
#include <cstdint>
#include <string>
#include <vector>

class FilterDatasets
{
 public:
  explicit FilterDatasets(bool debug = false);

  void process(const std::string &input, const std::string &output = ".");

 private:
  void analyze(const std::string &input, const std::string &outputDir);
  void readRunInfo(const std::string &line);

  std::string getCalibration(const std::string &pl_type, uint64_t iov);

  std::vector<std::pair<std::string, std::string>> m_runInfo;
  std::map<std::string, int> m_ctr;

  std::vector<std::string> m_cdbName = {"CEMC_BadTowerMap"
                                      , "CEMC_meanTime", "HCALIN_meanTime", "HCALOUT_meanTime"
                                      , "CEMC_hotTowers_fracBadChi2", "HCALIN_hotTowers_fracBadChi2", "HCALOUT_hotTowers_fracBadChi2"
                                      , "CEMC_ZSCrossCalib", "HCALIN_ZSCrossCalib", "HCALOUT_ZSCrossCalib"};

  bool m_debug;
};

#endif
