// Tell emacs that this is a C++ source
// -*- C++ -*-
#ifndef CALOTOWERTIMECALIBRATION_H
#define CALOTOWERTIMECALIBRATION_H

#include <caloreco/CaloTowerDefs.h>

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class CDBTTree;
class PHCompositeNode;
class TH1;
class TH2;

/**
 * Timing-only sidecar calibration for calibrated TowerInfo nodes.
 *
 * Default input node:
 *   TOWERINFO_CALIB_<detector>
 *
 * Default output node:
 *   TOWERINFO_CALIB_TIMING_<detector>
 *
 * The input node is never modified. The output container is copied from the
 * input container, preserving calibrated energy and tower metadata/status.
 * Only the output tower time is replaced for valid, non-ZS towers:
 *
 *   raw_time = standard_time + official_mean_time
 *
 *   slew(E) = slew_p0 + slew_p1 * exp(slew_p2 * calibrated_energy)
 *
 *   corrected_time = raw_time
 *                  + phase1_shift
 *                  - sector_offset
 *                  - tower_offset
 *                  - slew(E)
 *
 * QA is intentionally reduced to:
 *   - standard and corrected tower-time distributions;
 *   - standard and corrected energy-vs-time distributions;
 *   - ZS-tower time, energy, and energy-vs-time distributions.
 *
 * ZS towers are copied into the timing sidecar without a custom timing
 * correction. Their QA therefore uses the time and calibrated energy from the
 * standard input tower directly.
 *
 * There are two tower-quality populations:
 *   (1) all accepted non-ZS towers, with no get_isGood() requirement;
 *   (2) strict TowerInfo::get_isGood().
 *
 * QA histograms are filled only for events with a finite reconstructed global
 * z vertex, but no |zvtx| magnitude cut is applied.
 */
class CaloTowerTimeCalibration : public SubsysReco
{
 public:
  explicit CaloTowerTimeCalibration(
      const std::string &name = "CaloTowerTimeCalibration");
  ~CaloTowerTimeCalibration() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void set_detector_type(CaloTowerDefs::DetectorSystem detector)
  {
    m_detectorType = detector;
  }

  void set_inputNodePrefix(const std::string &prefix)
  {
    m_inputNodePrefix = prefix;
  }

  void set_outputNodePrefix(const std::string &prefix)
  {
    m_outputNodePrefix = prefix;
  }

  void set_inputNodeName(const std::string &name)
  {
    m_inputNodeName = name;
  }

  void set_outputNodeName(const std::string &name)
  {
    m_outputNodeName = name;
  }

  void set_meanTimeCalibName(const std::string &name)
  {
    m_meanTimeCalibName = name;
  }

  void set_timeCorrectionCalibName(const std::string &name)
  {
    m_timeCorrectionCalibName = name;
  }

  void set_directURL_timeCorrection(const std::string &url)
  {
    m_directTimeCorrectionURL = url;
  }

  void set_localTimeCorrectionDirectory(const std::string &directory)
  {
    m_localTimeCorrectionDirectory = directory;
  }

  void set_directURL_meanTime(const std::string &url)
  {
    m_directMeanTimeURL = url;
  }

  void set_doAbortMissingCalibration(bool value = true)
  {
    m_abortMissingCalibration = value;
  }

  void set_doQA(bool value = true)
  {
    m_doQA = value;
  }

  // QA threshold only; calibration is still applied below this energy.
  void set_qaEnergyThreshold(float value)
  {
    m_qaEnergyThreshold = value;
  }

 private:
  struct TimingInfo
  {
    float meanTime{0.0F};
    float phase1Shift{0.0F};
    float sectorOffset{0.0F};
    float towerOffset{0.0F};
    float slewP0{0.0F};
    float slewP1{0.0F};
    float slewP2{0.0F};
    bool valid{false};
  };

  struct TowerQASet
  {
    TH1 *hStandardTime{nullptr};
    TH1 *hCorrectedTime{nullptr};
    TH2 *hStandardEnergyVsTime{nullptr};
    TH2 *hCorrectedEnergyVsTime{nullptr};
  };

  bool ResolveDetector();
  void ResolveNames();
  std::string BuildLocalTimeCorrectionURL() const;
  void CreateNodeTree(PHCompositeNode *topNode);
  void LoadCalibration(PHCompositeNode *topNode);
  void CreateQAHistograms();

  CaloTowerDefs::DetectorSystem m_detectorType{CaloTowerDefs::HCALOUT};
  std::string m_detector{"HCALOUT"};

  std::string m_inputNodePrefix{"TOWERINFO_CALIB_"};
  std::string m_outputNodePrefix{"TOWERINFO_CALIB_TIMING_"};
  std::string m_inputNodeName;
  std::string m_outputNodeName;

  std::string m_meanTimeCalibName;
  std::string m_timeCorrectionCalibName;

  std::string m_directMeanTimeURL;
  std::string m_directTimeCorrectionURL;

  std::string m_localTimeCorrectionDirectory;

  bool m_abortMissingCalibration{true};
  bool m_calibrationAvailable{false};

  CDBTTree *m_meanTimeCDB{nullptr};
  CDBTTree *m_timeCorrectionCDB{nullptr};
  std::vector<TimingInfo> m_timingInfo;

  bool m_doQA{false};
  bool m_qaHistogramsInitialized{false};
  float m_qaEnergyThreshold;

  // ZS-tower QA. ZS timing is not custom-corrected.
  TH1 *m_hZSTime{nullptr};
  TH1 *m_hZSEnergy{nullptr};
  TH2 *m_hZSEnergyVsTime{nullptr};

  // Complete tower QA without a get_isGood() requirement.
  TowerQASet m_qaAllTowers;

  // Strict TowerInfo::get_isGood() QA.
  TowerQASet m_qaGoodTowers;
};

#endif  // CALOTOWERTIMECALIBRATION_H
