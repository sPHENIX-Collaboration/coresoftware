#ifndef RAWCLUSTERCNNCLASSIFIER_H
#define RAWCLUSTERCNNCLASSIFIER_H

#include <fun4all/SubsysReco.h>

#include <phool/onnxlib.h>

class PHCompositeNode;

class RawClusterCNNClassifier : public SubsysReco
{
 public:
  RawClusterCNNClassifier(const std::string &name = "RawClusterCNNClassifier");

  ~RawClusterCNNClassifier() override;

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  void Print(const std::string &what = "ALL") const override;

 private:
  Ort::Session *onnxmodule{nullptr};
  const int inputDimx = 5;
  const int inputDimy = 5;
  const int inputDimz = 1;
  const int outputDim = 1;

  const float minET = 3;
};

#endif