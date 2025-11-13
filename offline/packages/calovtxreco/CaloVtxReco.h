#ifndef CALOVTXRECO_CALOVTXRECO_H
#define CALOVTXRECO_CALOVTXRECO_H

#include <calobase/RawTowerDefs.h>

#include <fun4all/SubsysReco.h>

class CaloVertexMap;
class PHCompositeNode;
class RawTowerGeomContainer;
class TowerInfoContainer;

class CaloVtxReco : public SubsysReco
{
 public:
  CaloVtxReco(const std::string &name = "CaloVtxReco", const std::string &jetnodename = "zzjets06", const bool use_z_energy_dep = false);

  virtual ~CaloVtxReco() = default;

  int createNodes(PHCompositeNode *topNode);

  float new_eta(int channel, TowerInfoContainer *towers, RawTowerGeomContainer *geom, RawTowerDefs::CalorimeterId caloID, float testz);

  int InitRun(PHCompositeNode *topNode) override;

  int calo_tower_algorithm(PHCompositeNode *topNode) const;

  int process_event(PHCompositeNode *topNode) override;

  float get_jet_threshold() { return m_jet_threshold; }

  void set_jet_threshold(float new_thresh) { m_jet_threshold = new_thresh; }

  float get_calib_factor() { return m_calib_factor; }

  void set_calib_factor(float new_calib) { m_calib_factor = new_calib; }

  float get_energy_cut() { return m_energy_cut; }

  void set_energy_cut(float new_energy) { m_energy_cut = new_energy; }

 private:
  CaloVertexMap *m_calovtxmap{nullptr};
  float m_jet_threshold{15};
  float m_zvtx{std::numeric_limits<float>::quiet_NaN()};
  float m_calib_factor{1.406};
  float m_energy_cut{0.1};
  float m_radius_EM{std::numeric_limits<float>::quiet_NaN()};
  float m_radius_OH{std::numeric_limits<float>::quiet_NaN()};
  bool m_use_z_energy_dep;
  std::string m_jetnodename;
};

#endif  // CALOVTXRECO_CALOVTXRECO_H
