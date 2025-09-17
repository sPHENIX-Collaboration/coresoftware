#ifndef CALOPACKETSKIMMER_H
#define CALOPACKETSKIMMER_H

#include <caloreco/CaloTowerDefs.h>
#include <fun4all/SubsysReco.h>
#include <algorithm>
#include <string>
#include <vector>
#include <TH1D.h>

class PHCompositeNode;

class CaloPacketSkimmer : public SubsysReco
{
 public:
  explicit CaloPacketSkimmer(const std::string& name = "CaloPacketSkimmer");
  ~CaloPacketSkimmer() override = default;

  int process_event(PHCompositeNode* /*topNode*/) override;
  int InitRun(PHCompositeNode* topNode) override;

  int EndRun(const int /*runnumber*/) override;

  void set_all_detectors_off() { m_CaloDetectors.clear(); }

  void set_all_detectors_on()
  {
    m_CaloDetectors = {CaloTowerDefs::CEMC, CaloTowerDefs::HCALIN,
                       CaloTowerDefs::HCALOUT, CaloTowerDefs::SEPD,
                       CaloTowerDefs::ZDC, CaloTowerDefs::MBD};
  }

  void set_detector_on(CaloTowerDefs::DetectorSystem det)
  {
    if (std::find(m_CaloDetectors.begin(), m_CaloDetectors.end(), det) == m_CaloDetectors.end())
    {
      m_CaloDetectors.push_back(det);
    }
  }

  void set_offlineflag(const bool f = true)
  {
    m_UseOfflinePacketFlag = f;
  }

 private:
  std::vector<CaloTowerDefs::DetectorSystem> m_CaloDetectors{
      CaloTowerDefs::CEMC, CaloTowerDefs::HCALIN, CaloTowerDefs::HCALOUT,
      CaloTowerDefs::SEPD, CaloTowerDefs::ZDC, CaloTowerDefs::MBD};
  bool m_UseOfflinePacketFlag{true};
  bool m_PacketNodesFlag{false};

  int processDetector(PHCompositeNode* topNode, CaloTowerDefs::DetectorSystem det);

  std::string getDetectorName(CaloTowerDefs::DetectorSystem det)
  {
    switch (det)
    {
    case CaloTowerDefs::CEMC:
      return "CEMC";
    case CaloTowerDefs::HCALIN:
      return "HCALIN";
    case CaloTowerDefs::HCALOUT:
      return "HCALOUT";
    case CaloTowerDefs::SEPD:
      return "SEPD";
    case CaloTowerDefs::ZDC:
      return "ZDC";
    case CaloTowerDefs::MBD:
      return "MBD";
    default:
      return "UNKNOWN";
    }
  }

  std::pair<int, int> getPacketRange(CaloTowerDefs::DetectorSystem det)
  {
    switch (det)
    {
    case CaloTowerDefs::CEMC:
      return {6001, 6128};
    case CaloTowerDefs::HCALIN:
      return {7001, 7008};
    case CaloTowerDefs::HCALOUT:
      return {8001, 8008};
    case CaloTowerDefs::SEPD:
      return {9001, 9006};
    case CaloTowerDefs::ZDC:
      return {12001, 12001};
    case CaloTowerDefs::MBD:
      return {1001, 1002};
    default:
      return {0, 0};  // Invalid range
    }
  }

  TH1D* h_aborted_events{nullptr};
  TH1D* h_kept_events{nullptr};
  TH1D* h_missing_packets{nullptr};
  TH1D* h_empty_packets{nullptr};
};

#endif
