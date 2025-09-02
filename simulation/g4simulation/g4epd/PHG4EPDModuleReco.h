// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4EPD_PHG4EPDMODULERECO_H
#define G4EPD_PHG4EPDMODULERECO_H

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <array>
#include <limits>
#include <string>

class PHCompositeNode;

class PHG4EPDModuleReco : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4EPDModuleReco(const std::string &name = "EpdModuleReco");

  ~PHG4EPDModuleReco() override = default;

  //! module initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  //! Reset after every event
  int ResetEvent(PHCompositeNode * /*topNode*/) override;

  void SetDefaultParameters() override;

  void Detector(const std::string &detector);

  const std::string &get_epd_sim_tower_node_prefix() const
  {
    return m_EPDSimTowerNodePrefix;
  }

  void
  set_epd_sim_tower_node_prefix(const std::string &epdsimTowerNodePrefix)
  {
    m_EPDSimTowerNodePrefix = epdsimTowerNodePrefix;
  }

  const std::string &get_epd_calib_tower_node_prefix() const
  {
    return m_EPDCalibTowerNodePrefix;
  }

  void
  set_epd_calib_tower_node_prefix(const std::string &epdcalibTowerNodePrefix)
  {
    m_EPDCalibTowerNodePrefix = epdcalibTowerNodePrefix;
  }

  void set_timing_window(const double tmi, const double tma);

 private:
  void FillTilePhiArray();
  void FillTilePhi0Array();
  static int Getrmap(int rindex);
  static int Getphimap(int phiindex);
  float GetTilePhi(int thisphi);
  float GetTilePhi0(int thisphi0);
  static float GetTileR(int thisr);
  static float GetTileZ(int thisz);
  void CreateNodes(PHCompositeNode *topNode);

  double m_EpdMpv{std::numeric_limits<double>::quiet_NaN()};
  double tmin{std::numeric_limits<double>::quiet_NaN()};
  double tmax{std::numeric_limits<double>::quiet_NaN()};
  double m_DeltaT{std::numeric_limits<double>::quiet_NaN()};

  std::array<float, 24> m_tilephi{};
  std::array<float, 12> m_tilephi0{};
  std::array<std::array<std::array<double, 31>, 12>, 2> m_EpdTile_e = {};
  std::array<std::array<std::array<double, 31>, 12>, 2> m_EpdTile_Calib_e = {};

  std::string m_Detector;
  std::string m_Hitnodename;
  std::string m_EPDSimTowerNodePrefix{"SIM"};
  std::string m_EPDCalibTowerNodePrefix{"CALIB"};

  std::string m_TowerInfoNodeName;
  std::string m_TowerInfoNodeName_calib;
};

#endif
