// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EPD_EPDRECO_H__
#define EPD_EPDRECO_H__

//===========================================================
/// \author Ejiro Umaka
//===========================================================

#include <fun4all/SubsysReco.h>

#include <array>
#include <string> // for string

class CDBTTree;
class PHCompositeNode;

class EpdReco : public SubsysReco {
public:
  EpdReco(const std::string &name = "EpdReco");
  ~EpdReco() override = default;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

private:
  void CreateNodes(PHCompositeNode *topNode);
  float GetTilePhi(int thisphi);
  float GetTilePhi0(int thisphi0);
  float GetTileR(int thisr);
  float GetTileZ(int thisz);
  void FillTilePhiArray();
  void FillTilePhi0Array();

  CDBTTree *cdbttree{nullptr};
  bool m_overrideCalibName{false};
  bool m_overrideFieldName{false};

  std::string m_fieldname;
  std::string m_calibName;
  std::string m_Detector {"SEPD"};
  std::string m_TowerInfoNodeName_calib {"TOWERINFO_CALIB_SEPD"};

  std::array<float,24> tilephi{};
  std::array<float,12> tilephi0{};
};

#endif // EPD_EPDRECO_H
