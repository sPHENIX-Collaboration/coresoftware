// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef JETBASE_JETCALIB_H
#define JETBASE_JETCALIB_H

#include <calobase/TowerInfoContainer.h>
#include <fun4all/SubsysReco.h>
#include <iostream>
#include <limits>
#include <string>

class CDBTF;
class PHCompositeNode;
class TF1;

class JetCalib : public SubsysReco
{
 public:
  explicit JetCalib(const std::string &name = "JetCalib");
  ~JetCalib() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int CreateNodeTree(PHCompositeNode *topNode);

  // Getters.
  void set_InputNode(const std::string &inputNode) {m_inputNode = inputNode;}
  void set_OutputNode(const std::string &outputNode) {m_outputNode = outputNode;}
  
 private:
  // Functions.
  std::string fetchCalibDir(const char* calibType);
  float getTotalCorrFactor(TF1* JetCalibFunc, float jetPt);

  // Input.
  std::string m_inputNode{"AntiKt_Tower_r04"};
  std::string m_outputNode{"AntiKt_Tower_r04_Calib"};

  // Variables.
  CDBTF *m_JetCalibFile{nullptr};
  TF1 *m_JetCalibFunc{nullptr};
};

#endif  // JETBASE_JETCALIB_H