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
  void set_InputNode(const std::string &inputNode) { m_inputNode = inputNode; }
  void set_OutputNode(const std::string &outputNode) { m_outputNode = outputNode; }
  void set_JetRadius(float radius) { jet_radius = radius; }
  void set_ZvrtxNode(const std::string &zvrtxNode) { m_zvrtxNode = zvrtxNode; }
  void set_ApplyZvrtxDependentCalib(bool apply) { ApplyZvrtxDependentCalib = apply; }
  void set_ApplyEtaDependentCalib(bool apply) { ApplyEtaDependentCalib = apply; }

 private:
  // Functions.
  static std::string fetchCalibDir(const char *calibType);
  float doCalibration(const std::vector<std::vector<TF1 *>> &JetCalibFunc, float jetPt, float zvrtx, float eta) const;

  // Input.
  std::string m_inputNode{"AntiKt_Tower_r04"};
  std::string m_outputNode{"AntiKt_Tower_r04_Calib"};
  std::string m_zvrtxNode{"GlobalVertexMap"};
  float jet_radius{0.4};                 // Jet radius.
  bool ApplyZvrtxDependentCalib{false};  // Apply Z-vertex dependent calibration.
  bool ApplyEtaDependentCalib{false};    // Apply eta dependent calibration.

  // Variables.
  CDBTF *m_JetCalibFile{nullptr};
  std::vector<std::vector<TF1 *>> m_JetCalibFunc{};
};

#endif  // JETBASE_JETCALIB_H
