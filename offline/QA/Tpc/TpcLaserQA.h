#ifndef QA_TPC_TPCLASERQA_H
#define QA_TPC_TPCLASERQA_H

#include <fun4all/SubsysReco.h>

#include <tpc/LaserClusterHelper.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <string>

class PHCompositeNode;
class TH1;
class TH2;

class TpcLaserQA : public SubsysReco
{
 public:
  explicit TpcLaserQA(const std::string& name = "TpcLaserQA");

  ~TpcLaserQA() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;

  void set_useZ(bool use) { m_useZ = use; }

 private:
  void createHistos();
  std::string getHistoPrefix() const;

  double rBinEdges[5] = {0.0, 0.256, 0.504, 0.752, 1.0};

  TH1* m_nLaserEvents{nullptr};
  TH2* m_TPCWheel[2]{nullptr};

  TH1* m_nLaserClusters[2]{nullptr};
  TH2* m_saturation[2]{nullptr};

  TH1* m_sample_R1[2][12]{{nullptr}};
  TH1* m_sample_R2[2][12]{{nullptr}};
  TH1* m_sample_R3[2][12]{{nullptr}};

  ActsGeometry *m_tGeometry{nullptr};
  PHG4TpcGeomContainer *m_geom_container{nullptr};

  LaserClusterHelper m_laserClusterHelper;
  bool m_useZ{false};
};

#endif
