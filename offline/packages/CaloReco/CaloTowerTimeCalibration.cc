#include "CaloTowerTimeCalibration.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <qautils/QAHistManagerDef.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
// Towers with invalid timing calibration are assigned NaN.
// Downstream analyses should require std::isfinite(tower->get_time())
// before using tower timing.
//
// ZS tower timing is not recalibrated; the standard calibrated
// tower time is copied unchanged.
namespace
{
  constexpr float InvalidSentinelLimit = -998.5F;

  bool IsUsableConstant(float value)
  {
    return std::isfinite(value) && value > InvalidSentinelLimit;
  }

  std::string FoundText(const std::string &url)
  {
    return url.empty() ? " [missing]" : " [found]";
  }

  template <class HistogramType>
  void StyleHistogram(HistogramType *histogram)
  {
    if (!histogram)
    {
      return;
    }

    // Store HIST as the default ROOT draw option for every QA histogram.
    // This keeps bin outlines/steps visible instead of point/error-bar drawing
    // when the histogram is opened without an explicit Draw() option.
    histogram->SetOption("HIST");

    // ROOT font 62 is Helvetica bold. Use a slightly smaller stored axis
    // style so long labels like "(sample)" fit cleanly when users open the
    // histograms directly from the ROOT file.
    constexpr int boldFont = 62;
    constexpr float axisTitleSize = 0.045F;
    constexpr float axisLabelSize = 0.036F;
    constexpr float xTitleOffset = 0.95F;
    constexpr float yTitleOffset = 0.90F;
    constexpr float zTitleOffset = 0.95F;

    TAxis *axes[] = {
        histogram->GetXaxis(),
        histogram->GetYaxis(),
        histogram->GetZaxis()};

    for (TAxis *axis : axes)
    {
      if (!axis)
      {
        continue;
      }

      axis->SetTitleFont(boldFont);
      axis->SetLabelFont(boldFont);
      axis->SetTitleSize(axisTitleSize);
      axis->SetLabelSize(axisLabelSize);
    }

    if (histogram->GetXaxis())
    {
      histogram->GetXaxis()->SetTitleOffset(xTitleOffset);
    }
    if (histogram->GetYaxis())
    {
      histogram->GetYaxis()->SetTitleOffset(yTitleOffset);
    }
    if (histogram->GetZaxis())
    {
      histogram->GetZaxis()->SetTitleOffset(zTitleOffset);
    }
  }
}  // namespace

CaloTowerTimeCalibration::CaloTowerTimeCalibration(const std::string &name)
  : SubsysReco(name)
  , m_qaEnergyThreshold(0.3F)
{
}

CaloTowerTimeCalibration::~CaloTowerTimeCalibration()
{
  delete m_meanTimeCDB;
  delete m_timeCorrectionCDB;
}

bool CaloTowerTimeCalibration::ResolveDetector()
{
  if (m_detectorType == CaloTowerDefs::CEMC)
  {
    m_detector = "CEMC";
    return true;
  }
  if (m_detectorType == CaloTowerDefs::HCALIN)
  {
    m_detector = "HCALIN";
    return true;
  }
  if (m_detectorType == CaloTowerDefs::HCALOUT)
  {
    m_detector = "HCALOUT";
    return true;
  }

  return false;
}

void CaloTowerTimeCalibration::ResolveNames()
{
  if (!m_inputNodePrefix.empty())
  {
    m_inputNodeName = m_inputNodePrefix + m_detector;
  }

  if (!m_outputNodePrefix.empty())
  {
    m_outputNodeName = m_outputNodePrefix + m_detector;
  }

  if (!m_meanTimeCalibName.empty())
  {
    m_meanTimeCalibName = m_detector + "_meanTime";
  }

  if (!m_timeCorrectionCalibName.empty())
  {
    // Central-CDB fallback used after an official domain is available.
    m_timeCorrectionCalibName = m_detector + "_timeCalib_MV_test";
  }
}

std::string CaloTowerTimeCalibration::BuildLocalTimeCorrectionURL() const
{
  if (m_localTimeCorrectionDirectory.empty())
  {
    return {};
  }

  recoConsts *recoConst = recoConsts::instance();
  const uint64_t runnumber =
      recoConst->get_uint64Flag("TIMESTAMP");

  if (runnumber == 0)
  {
    return {};
  }

  std::string directory = m_localTimeCorrectionDirectory;
  while (!directory.empty() && directory.back() == '/')
  {
    directory.pop_back();
  }

  return directory
         + "/"
         + m_detector
         + "_timeCalib_MV_test_run"
         + std::to_string(runnumber)
         + ".root";
}

int CaloTowerTimeCalibration::InitRun(PHCompositeNode *topNode)
{
  if (!ResolveDetector())
  {
    std::cerr << Name() << ": unsupported detector type" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  ResolveNames();

  delete m_meanTimeCDB;
  m_meanTimeCDB = nullptr;
  delete m_timeCorrectionCDB;
  m_timeCorrectionCDB = nullptr;
  m_timingInfo.clear();
  m_calibrationAvailable = false;

  const std::string meanTimeURL =
    !m_directMeanTimeURL.empty()
          ? m_directMeanTimeURL
          : CDBInterface::instance()->getUrl(m_meanTimeCalibName);

  std::string timeCorrectionURL;
  std::string attemptedLocalURL;

  if (!m_directTimeCorrectionURL.empty())
  {
    timeCorrectionURL = m_directTimeCorrectionURL;
    std::cout << Name() << "::" << m_detector
              << ": using direct custom timing payload "
              << timeCorrectionURL << std::endl;
  }
  else
  {
    attemptedLocalURL = BuildLocalTimeCorrectionURL();

    if (!attemptedLocalURL.empty()
        && !gSystem->AccessPathName(attemptedLocalURL.c_str()))
    {
      timeCorrectionURL = attemptedLocalURL;
      std::cout << Name() << "::" << m_detector
                << ": using local custom timing payload "
                << timeCorrectionURL << std::endl;
    }
    else
    {
      if (!attemptedLocalURL.empty())
      {
        std::cout << Name() << "::" << m_detector
                  << ": local custom timing payload not found: "
                  << attemptedLocalURL << std::endl;
      }

      // Future fallback after the payload has an official CDB domain.
      timeCorrectionURL =
          CDBInterface::instance()->getUrl(m_timeCorrectionCalibName);
    }
  }

  if (meanTimeURL.empty() || timeCorrectionURL.empty())
  {
    std::cerr << Name() << "::" << m_detector
              << ": missing required timing calibration" << std::endl
              << "  mean-time CDB " << m_meanTimeCalibName
              << FoundText(meanTimeURL) << std::endl;

    if (!attemptedLocalURL.empty())
    {
      std::cerr << "  attempted local payload "
                << attemptedLocalURL << " [missing]" << std::endl;
    }

    std::cerr << "  custom timing CDB " << m_timeCorrectionCalibName
              << FoundText(timeCorrectionURL) << std::endl;

    if (m_abortMissingCalibration)
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }

    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_meanTimeCDB = new CDBTTree(meanTimeURL);
  m_timeCorrectionCDB = new CDBTTree(timeCorrectionURL);

  try
  {
    CreateNodeTree(topNode);
    CreateQAHistograms();
    LoadCalibration(topNode);
  }
  catch (const std::exception &error)
  {
    std::cerr << Name() << "::" << m_detector
              << ": " << error.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_calibrationAvailable = true;

  std::cout << Name() << "::" << m_detector
            << " input=" << m_inputNodeName
            << " output=" << m_outputNodeName
            << " custom_payload=" << timeCorrectionURL
            << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTowerTimeCalibration::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iterator(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(
      iterator.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    throw std::runtime_error("DST node is missing");
  }

  TowerInfoContainer *inputTowers =
      findNode::getClass<TowerInfoContainer>(topNode, m_inputNodeName);
  if (!inputTowers)
  {
    throw std::runtime_error(
        "missing input node " + m_inputNodeName
        + "; register this subsystem after standard CaloTowerCalib");
  }

  PHCompositeNode *detectorNode = dynamic_cast<PHCompositeNode *>(
      iterator.findFirst("PHCompositeNode", m_detector));
  if (!detectorNode)
  {
    detectorNode = new PHCompositeNode(m_detector);
    dstNode->addNode(detectorNode);
  }

  TowerInfoContainer *outputTowers =
      findNode::getClass<TowerInfoContainer>(topNode, m_outputNodeName);

  if (!outputTowers)
  {
    outputTowers =
        dynamic_cast<TowerInfoContainer *>(inputTowers->CloneMe());
    if (!outputTowers)
    {
      throw std::runtime_error("failed to clone input tower container");
    }

    auto *outputNode = new PHIODataNode<PHObject>(
        outputTowers,
        m_outputNodeName,
        "PHObject");

    detectorNode->addNode(outputNode);
  }
}

void CaloTowerTimeCalibration::LoadCalibration(PHCompositeNode *topNode)
{
  TowerInfoContainer *inputTowers =
      findNode::getClass<TowerInfoContainer>(topNode, m_inputNodeName);
  if (!inputTowers)
  {
    throw std::runtime_error("missing input node " + m_inputNodeName);
  }

  if (!m_meanTimeCDB || !m_timeCorrectionCDB)
  {
    throw std::runtime_error("timing CDB objects were not created");
  }

  const unsigned int numberOfTowers = inputTowers->size();
  m_timingInfo.resize(numberOfTowers);

  const int payloadTowerCount =
      m_timeCorrectionCDB->GetSingleIntValue("ntowers");
  const int formulaVersion =
      m_timeCorrectionCDB->GetSingleIntValue("formula_version");
  const int payloadRun =
      m_timeCorrectionCDB->GetSingleIntValue("runnumber");

  recoConsts *recoConst = recoConsts::instance();
  const int requestedRun = static_cast<int>(
      recoConst->get_uint64Flag("TIMESTAMP"));

  if (payloadTowerCount > 0
      && static_cast<unsigned int>(payloadTowerCount) != numberOfTowers)
  {
    throw std::runtime_error(
        "payload ntowers does not match " + m_inputNodeName);
  }

  if (payloadRun > 0 && requestedRun > 0 && payloadRun != requestedRun)
  {
    throw std::runtime_error(
        "payload run " + std::to_string(payloadRun)
        + " does not match TIMESTAMP " + std::to_string(requestedRun));
  }

  unsigned int invalidTowers = 0;

  for (unsigned int channel = 0;
       channel < numberOfTowers;
       ++channel)
  {
    const unsigned int encodedKey = inputTowers->encode_key(channel);
    TimingInfo &timing = m_timingInfo[channel];

    // Standard CaloTowerCalib used this value as:
    // standard_time = raw_time - mean_time.
    timing.meanTime =
        m_meanTimeCDB->GetFloatValue(encodedKey, "time");

    // MakeAllTimingCDB_MV_test writes these values by flat channel.
    timing.phase1Shift =
        m_timeCorrectionCDB->GetFloatValue(channel, "phase1_shift");
    timing.sectorOffset =
        m_timeCorrectionCDB->GetFloatValue(channel, "sector_offset");
    timing.towerOffset =
        m_timeCorrectionCDB->GetFloatValue(channel, "tower_offset");
    timing.slewP0 =
        m_timeCorrectionCDB->GetFloatValue(channel, "slew_p0");
    timing.slewP1 =
        m_timeCorrectionCDB->GetFloatValue(channel, "slew_p1");
    timing.slewP2 =
        m_timeCorrectionCDB->GetFloatValue(channel, "slew_p2");

    bool valid =
        IsUsableConstant(timing.meanTime)
        && IsUsableConstant(timing.phase1Shift)
        && IsUsableConstant(timing.sectorOffset)
        && IsUsableConstant(timing.towerOffset)
        && IsUsableConstant(timing.slewP0)
        && IsUsableConstant(timing.slewP1)
        && IsUsableConstant(timing.slewP2);

    if (formulaVersion >= 3)
    {
      valid = valid
              && m_timeCorrectionCDB->GetIntValue(
                     channel,
                     "timing_valid") != 0;
    }

    timing.valid = valid;

    if (!valid)
    {
      ++invalidTowers;
    }

  }

  std::cout << Name() << "::" << m_detector
            << " loaded " << numberOfTowers
            << " channels; invalid=" << invalidTowers
            << ", formula_version=" << formulaVersion
            << std::endl;
}

int CaloTowerTimeCalibration::process_event(PHCompositeNode *topNode)
{
  if (!m_calibrationAvailable)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  TowerInfoContainer *standardTowers =
      findNode::getClass<TowerInfoContainer>(topNode, m_inputNodeName);
  TowerInfoContainer *timingTowers =
      findNode::getClass<TowerInfoContainer>(topNode, m_outputNodeName);

  if (!standardTowers || !timingTowers)
  {
    std::cerr << Name() << "::" << m_detector
              << ": required tower node is missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (standardTowers->size() != timingTowers->size()
      || standardTowers->size() != m_timingInfo.size())
  {
    std::cerr << Name() << "::" << m_detector
              << ": tower-container size mismatch" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  const unsigned int numberOfTowers = standardTowers->size();

  // QA is filled only for events with a finite reconstructed global z vertex.
  // There is no |zvtx| magnitude cut and no zvtx-vs-time histogram.
  bool haveValidZVertex = false;

  if (m_doQA)
  {
    GlobalVertexMap *vertexMap =
        findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

    if (vertexMap && !vertexMap->empty())
    {
      GlobalVertex *vertex = vertexMap->begin()->second;
      if (vertex)
      {
        haveValidZVertex = std::isfinite(vertex->get_z());
      }
    }
  }

  for (unsigned int channel = 0;
       channel < numberOfTowers;
       ++channel)
  {
    TowerInfo *standardTower =
        standardTowers->get_tower_at_channel(channel);
    TowerInfo *timingTower =
        timingTowers->get_tower_at_channel(channel);

    if (!standardTower || !timingTower)
    {
      continue;
    }

    // Preserve calibrated energy and all metadata/status in the sidecar node.
    timingTower->copy_tower(standardTower);

    // ZS timing is not custom-corrected. Fill dedicated ZS QA directly from
    // the standard calibrated tower, then leave the copied sidecar unchanged.
    // ZS QA is restricted to towers that pass get_isGood(). Do not apply
    // m_qaEnergyThreshold here: all energies of good ZS towers are retained.
    if (standardTower->get_isZS())
    {
      if (m_doQA && haveValidZVertex && standardTower->get_isGood())
      {
        const float zsTime = standardTower->get_time();
        const float zsEnergy = standardTower->get_energy();

        if (std::isfinite(zsTime) && m_hZSTime)
        {
          m_hZSTime->Fill(zsTime);
        }

        if (std::isfinite(zsEnergy) && m_hZSEnergy)
        {
          m_hZSEnergy->Fill(zsEnergy);
        }

        if (std::isfinite(zsTime)
            && std::isfinite(zsEnergy)
            && m_hZSEnergyVsTime)
        {
          // Time is the x axis, matching the rest of the timing QA.
          m_hZSEnergyVsTime->Fill(zsTime, zsEnergy);
        }
      }

      continue;
    }

    const TimingInfo &timing = m_timingInfo[channel];
    if (!timing.valid)
    {
      timingTower->set_time(std::numeric_limits<float>::quiet_NaN());
      continue;
    }

    const float standardTime = standardTower->get_time();
    const float calibratedEnergy = standardTower->get_energy();

    if (!std::isfinite(standardTime)
        || !std::isfinite(calibratedEnergy))
    {
      timingTower->set_time(std::numeric_limits<float>::quiet_NaN());
      continue;
    }

    // Standard CaloTowerCalib stores raw_time - official_mean_time.
    const float reconstructedRawTime =
        standardTime + timing.meanTime;

    const float slewCorrection =
        timing.slewP0
        + timing.slewP1
              * std::exp(timing.slewP2 * calibratedEnergy);

    const float correctedTime =
        reconstructedRawTime
        + timing.phase1Shift
        - timing.sectorOffset
        - timing.towerOffset
        - slewCorrection;

    if (!std::isfinite(reconstructedRawTime)
        || !std::isfinite(slewCorrection)
        || !std::isfinite(correctedTime))
    {
      timingTower->set_time(std::numeric_limits<float>::quiet_NaN());
      continue;
    }

    // Calibration is written before any QA quality definition is applied.
    timingTower->set_time(correctedTime);

    if (!m_doQA
        || !haveValidZVertex
        || calibratedEnergy < m_qaEnergyThreshold)
    {
      continue;
    }

    const auto fillTowerQA =
        [&](TowerQASet &qa)
        {
          qa.hStandardTime->Fill(standardTime);
          qa.hCorrectedTime->Fill(correctedTime);

          // Time is always the x axis.
          qa.hStandardEnergyVsTime->Fill(
              standardTime,
              calibratedEnergy);
          qa.hCorrectedEnergyVsTime->Fill(
              correctedTime,
              calibratedEnergy);

        };

    // No get_isGood() requirement.
    fillTowerQA(m_qaAllTowers);

    if (standardTower->get_isGood())
    {
      fillTowerQA(m_qaGoodTowers);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTowerTimeCalibration::CreateQAHistograms()
{
  if (!m_doQA || m_qaHistogramsInitialized)
  {
    return;
  }

  Fun4AllHistoManager *histogramManager =
      QAHistManagerDef::getHistoManager();
  if (!histogramManager)
  {
    throw std::runtime_error(
        Name() + "::" + m_detector
        + ": QAHistManagerDef returned a null histogram manager");
  }

  const std::string prefix =
      "h_CaloTowerTimeCalibration_" + m_detector;
  const std::string title =
      "CaloTowerTimeCalibration " + m_detector;

  constexpr int timeBins = 2000;
  constexpr float timeMin = -10.0F;
  constexpr float timeMax = 10.0F;
  constexpr int energyBins = 250;
  constexpr float energyMin = 0.0F;
  constexpr float energyMax = 25.0F;

  // ZS towers are concentrated at low energy. Use 10 MeV (0.01 GeV) bins
  // over 0--1 GeV for the dedicated ZS energy QA.
  constexpr int zsEnergyBins = 100;
  constexpr float zsEnergyMin = 0.0F;
  constexpr float zsEnergyMax = 1.0F;

  // Dedicated good-ZS QA. The timing sidecar intentionally does not modify ZS
  // tower times, so these are filled from standard towers passing get_isGood().
  m_hZSTime = new TH1F(
      (prefix + "_ZS_time").c_str(),
      (title + ", good ZS towers;ZS tower time (sample);towers").c_str(),
      timeBins,
      timeMin,
      timeMax);

  m_hZSEnergy = new TH1F(
      (prefix + "_ZS_energy").c_str(),
      (title + ", good ZS towers;ZS tower calibrated energy (GeV);towers").c_str(),
      zsEnergyBins,
      zsEnergyMin,
      zsEnergyMax);

  m_hZSEnergyVsTime = new TH2F(
      (prefix + "_ZS_energy_vs_time").c_str(),
      (title
       + ", good ZS towers;ZS tower time (sample);ZS tower calibrated energy (GeV)")
          .c_str(),
      timeBins,
      timeMin,
      timeMax,
      zsEnergyBins,
      zsEnergyMin,
      zsEnergyMax);

  StyleHistogram(m_hZSTime);
  StyleHistogram(m_hZSEnergy);
  StyleHistogram(m_hZSEnergyVsTime);

  histogramManager->registerHisto(m_hZSTime);
  histogramManager->registerHisto(m_hZSEnergy);
  histogramManager->registerHisto(m_hZSEnergyVsTime);

  const auto bookTowerQASet =
      [&](TowerQASet &qa,
          const std::string &nameSuffix,
          const std::string &titleSuffix)
      {
        const std::string qaPrefix = prefix + nameSuffix;
        const std::string qaTitle = title + titleSuffix;

        qa.hStandardTime = new TH1F(
            (qaPrefix + "_standard_time").c_str(),
            (qaTitle + ";standard time (sample);towers").c_str(),
            timeBins,
            timeMin,
            timeMax);

        qa.hCorrectedTime = new TH1F(
            (qaPrefix + "_corrected_time").c_str(),
            (qaTitle + ";custom corrected time (sample);towers").c_str(),
            timeBins,
            timeMin,
            timeMax);

        qa.hStandardEnergyVsTime = new TH2F(
            (qaPrefix + "_standard_energy_vs_time").c_str(),
            (qaTitle + ";standard time (sample);calibrated energy (GeV)").c_str(),
            timeBins,
            timeMin,
            timeMax,
            energyBins,
            energyMin,
            energyMax);

        qa.hCorrectedEnergyVsTime = new TH2F(
            (qaPrefix + "_corrected_energy_vs_time").c_str(),
            (qaTitle + ";custom corrected time (sample);calibrated energy (GeV)").c_str(),
            timeBins,
            timeMin,
            timeMax,
            energyBins,
            energyMin,
            energyMax);

        StyleHistogram(qa.hStandardTime);
        StyleHistogram(qa.hCorrectedTime);
        StyleHistogram(qa.hStandardEnergyVsTime);
        StyleHistogram(qa.hCorrectedEnergyVsTime);

        histogramManager->registerHisto(qa.hStandardTime);
        histogramManager->registerHisto(qa.hCorrectedTime);
        histogramManager->registerHisto(qa.hStandardEnergyVsTime);
        histogramManager->registerHisto(qa.hCorrectedEnergyVsTime);

        // No |zvtx| magnitude-cut histograms are booked. These inclusive
        // histograms are filled only when the event has a finite global z vertex.
      };

  bookTowerQASet(
      m_qaAllTowers,
      "",
      ", all accepted towers (no get_isGood cut)");

  bookTowerQASet(
      m_qaGoodTowers,
      "_isGood",
      ", strict TowerInfo::get_isGood() towers");

  // Intentionally no Sumw2 calls; all histograms are filled with unit weight.
  m_qaHistogramsInitialized = true;

  std::cout << Name() << "::" << m_detector
            << " booked all-tower and strict-get_isGood energy/time QA sets"
            << "; plus get_isGood ZS time, energy, and energy-vs-time QA"
            << "; QA requires a finite global z vertex but applies no |zvtx| cut"
            << std::endl;
}
