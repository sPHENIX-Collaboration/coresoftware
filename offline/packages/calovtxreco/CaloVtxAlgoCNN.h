#ifndef CALOVTXALGOCNN_H
#define CALOVTXALGOCNN_H

#include "CaloVtxAlgo.h"

#include <array>
#include <memory>
#include <string>

class PHCompositeNode;

// Calo-image CNN (ONNX) as a CaloVtxAlgo. Input is the 3x24x64 tower-energy image; preprocessing and calibration are inside the graph.
class CaloVtxAlgoCNN : public CaloVtxAlgo
{
 public:
  static constexpr int kNLayer = 3;
  static constexpr int kNEtaImg = 24;
  static constexpr int kNPhiImg = 64;
  static constexpr int kNEtaEmcFine = 96;
  static constexpr int kNPhiEmcFine = 256;
  enum Layer
  {
    kEMC = 0,
    kIHC = 1,
    kOHC = 2
  };

  explicit CaloVtxAlgoCNN() = default;
  ~CaloVtxAlgoCNN() override; 

  int Init(PHCompositeNode *topNode) override;
  int CalculateVertex(PHCompositeNode *topNode, float &zvtx) override;
  std::string Name() const override { return "CNN"; }
  VertexDefs::CALOALGO Algo() const override { return VertexDefs::CALOALGO::CNN; }

  void setModelFile(const std::string &path) { m_modelFile = path; }
  void setTowerEMin(Layer layer, float e) { m_towerEMin.at(layer) = e; }  // [GeV]
  void setMinTotalEnergy(float e) { m_minTotalEnergy = e; }               // [GeV]

 private:
  struct OnnxSession;  // pImpl, defined in the .cc

  int fillTowerImage(PHCompositeNode *topNode);
  int fillHcalLayer(PHCompositeNode *topNode, int layer);
  int fillEmcRetower(PHCompositeNode *topNode);
  int buildEmcRetowerMap(PHCompositeNode *topNode);
  bool predict(float &z);

  std::string m_modelFile{"vertex_cnn.onnx"};
  std::array<float, kNLayer> m_towerEMin{{0.068, 0.005, 0.035}};
  float m_minTotalEnergy{0.};

  std::unique_ptr<OnnxSession> m_onnx;

  std::array<std::array<std::array<float, kNPhiImg>, kNEtaImg>, kNLayer> m_image{};
  std::array<std::array<double, kNPhiEmcFine>, kNEtaEmcFine> m_rawEmcFine{};
  bool m_retowerMapReady{false};
  int m_retowerPhiOffset{-1};
  std::array<int, kNEtaImg> m_retowerLowerEta{};
  std::array<int, kNEtaImg> m_retowerUpperEta{};
  std::array<double, kNEtaImg> m_retowerLowerFrac{};
  std::array<double, kNEtaImg> m_retowerUpperFrac{};

  bool m_warnedMissingTowers{false};
  bool m_warnedRetowerMap{false};
  bool m_warnedPredict{false};
};

#endif
