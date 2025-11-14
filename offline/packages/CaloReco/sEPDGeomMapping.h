#ifndef CALORECO_SEPDGEOMMAPPING_H
#define CALORECO_SEPDGEOMMAPPING_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>

class PHCompositeNode;
class RawTowerGeomContainer;

class sEPDGeomMapping : public SubsysReco
{
 public:
  sEPDGeomMapping(const std::string &name = "sEPDGeomMapping");
  ~sEPDGeomMapping() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *) override { return 0; }  // Geometry only needs InitRun

 private:
  void CreateGeometry(PHCompositeNode *topNode);

  // SEPD geometry helper functions (from EpdReco)
  float GetTilePhi(int phibin) const;
  float GetTilePhi0(int phibin) const;
  float GetTileR(int rbin) const;
  float GetTileZ(int arm) const;
  int GetRmap(int tileindex) const;
  int GetPhimap(int tileindex) const;
  void FillTilePhiArray();
  void FillTilePhi0Array();

  // Geometry arrays for SEPD tiles
  std::array<float, 24> tilephi{};
  std::array<float, 12> tilephi0{};

  RawTowerGeomContainer *m_RawTowerGeomContainer = nullptr;
};

#endif  // CALORECO_SEPDGEOMMAPPING_H