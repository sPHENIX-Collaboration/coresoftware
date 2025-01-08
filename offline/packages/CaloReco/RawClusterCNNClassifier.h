#ifndef RAWCLUSTERCNNCLASSIFIER_H
#define RAWCLUSTERCNNCLASSIFIER_H

#include <fun4all/SubsysReco.h>

#include <phool/onnxlib.h>

class PHCompositeNode;
class RawClusterContainer;

class RawClusterCNNClassifier : public SubsysReco
{
 public:
  RawClusterCNNClassifier(const std::string &name = "RawClusterCNNClassifier");

  ~RawClusterCNNClassifier() override;

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  void Print(const std::string &what = "ALL") const override;
  
  void set_modelPath(const std::string &modelPath) { m_modelPath = modelPath; }

  void set_inputNodeName(const std::string &inputNodeName) { m_inputNodeName = inputNodeName; }

  void set_outputNodeName(const std::string &outputNodeName) { m_outputNodeName = outputNodeName; }

  void set_towerNodeName(const std::string &towerNodeName) { m_towerNodeName = towerNodeName; }

  void set_min_cluster_e(const float min_cluster_e) { m_min_cluster_e = min_cluster_e; }

 private:
  Ort::Session *onnxmodule{nullptr};
  const int inputDimx{5};
  const int inputDimy{5};
  const int inputDimz{1};
  const int outputDim{1};

  const float minET{3};

  std::string m_modelPath{"/sphenix/u/shuhang98/core_patch/coresoftware/offline/packages/CaloReco/functional_model.onnx"};

  std::string m_inputNodeName{"CLUSTERINFO_CEMC"};
  std::string m_outputNodeName{"CLUSTERINFO_CEMC_CNN"};
  std::string m_towerNodeName{"TOWERINFO_CALIB_CEMC"};

  bool inplace{false};

  RawClusterContainer *_clusters{nullptr};

  float m_min_cluster_e{3};

  void CreateNodes(PHCompositeNode* topNode);


};

#endif