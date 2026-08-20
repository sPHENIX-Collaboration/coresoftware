#ifndef CALOVTXALGOVIT_H
#define CALOVTXALGOVIT_H

#include "CaloVtxAlgo.h"

#include <array>
#include <memory>
#include <string>
#include <vector>

class PHCompositeNode;

// Calo vision-transformer (CaloViTv2_gate, ONNX) as a CaloVtxAlgo.
// Inputs are three 2-channel tower images at native tower granularity:
//   emcal (1,2,96,256), ihcal (1,2,24,64), ohcal (1,2,24,64)
// channel 0 = raw tower energy (negative energies floored at 0, NO log1p --
// the graph applies it), channel 1 = tower time. Output is pred_z [cm].
// The EMCal is used at fine 96x256 granularity, so no retowering is needed.
class CaloVtxAlgoVit : public CaloVtxAlgo
{
 public:
  static constexpr int kNCalo = 3;
  static constexpr int kNChan = 2;  // 0 = energy, 1 = time
  enum Calo
  {
    kEMC = 0,
    kIHC = 1,
    kOHC = 2
  };
  static constexpr std::array<int, kNCalo> kNEta{{96, 24, 24}};
  static constexpr std::array<int, kNCalo> kNPhi{{256, 64, 64}};

  // ctor/dtor defined in the .cc, where OnnxSession is a complete type
  // (required for the std::unique_ptr pImpl member, notably under cling)
  CaloVtxAlgoVit();
  ~CaloVtxAlgoVit() override;

  int Init(PHCompositeNode *topNode) override;
  int CalculateVertex(PHCompositeNode *topNode, float &zvtx) override;
  std::string Name() const override { return "ViT"; }
  VertexDefs::CALOALGO Algo() const override { return VertexDefs::CALOALGO::VIT; }

  void setModelFile(const std::string &path) { m_modelFile = path; }
  void setTowerNode(Calo calo, const std::string &node) { m_towerNode.at(calo) = node; }
  // tensor names in the exported graph (defaults match the CaloViTv2_gate export)
  void setInputName(Calo calo, const std::string &name) { m_inputName.at(calo) = name; }
  void setOutputName(const std::string &name) { m_outputName = name; }
  void setMinTotalEnergy(float e) { m_minTotalEnergy = e; }  // [GeV]
  void setUseGoodTowersOnly(bool b) { m_useGoodTowersOnly = b; }

 private:
  struct OnnxSession;  // pImpl, defined in the .cc

  int fillInputs(PHCompositeNode *topNode);
  int fillCalo(PHCompositeNode *topNode, int calo);
  bool predict(float &z);

  std::string m_modelFile{"calovit.onnx"};
  std::array<std::string, kNCalo> m_towerNode{{"TOWERINFO_CALIB_CEMC", "TOWERINFO_CALIB_HCALIN", "TOWERINFO_CALIB_HCALOUT"}};
  std::array<std::string, kNCalo> m_inputName{{"emcal", "ihcal", "ohcal"}};
  std::string m_outputName{"pred_z"};
  float m_minTotalEnergy{0.};
  bool m_useGoodTowersOnly{true};

  std::unique_ptr<OnnxSession> m_onnx;

  // flat (1,2,eta,phi) input tensors, one per calorimeter
  std::array<std::vector<float>, kNCalo> m_input;
  double m_etot{0.};

  bool m_warnedMissingTowers{false};
  bool m_warnedPredict{false};
};

#endif
