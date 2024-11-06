// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef JETBASE_JETCALIB_H
#define JETBASE_JETCALIB_H

#include <calobase/TowerInfoContainer.h>  // for TowerInfoContainer, TowerIn...

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

  void set_DoEtaJes(int doJes) {m_doEtaJes = doJes;}
  void set_doInsitu(int doInsitu) {m_doInsitu = doInsitu;}
  
 private:
  
  int getEtaBin(float eta);
  std::string fetchCalibDir(const char* calibType);
  float getTotalCorrFactor(TF1* etaJesTF, TF1 *rTrkTF, TF1 *gammeJetTF, float jetPt, float zvtx);

  //ints and floats
  int m_isEtaDependent{0};
  int m_doEtaJes{0};
  int m_doInsitu{0};
  static const int m_nEtaBins{4};
  int m_radius{4};
  
  //kinematic selections
  float m_etaStart{-0.7};
  float m_etaEnd{0.7};
  float m_zvtx{std::numeric_limits<float>::quiet_NaN()};
  std::vector<float> m_etaBins;
  
  //calibration designators
  std::string  m_jetType{nullptr};
  std::string m_jetInstrument{nullptr};
  int m_doBackgroundSub{1};
  int m_calibyear{24};

  //CDBTFs
  CDBTF *m_etaJesFile{nullptr};
  CDBTF *m_rTrkFile{nullptr};
  CDBTF *m_gammaJetFile{nullptr};

  TF1 *m_etaJesFunc[m_nEtaBins] {nullptr};
  TF1 *m_gammaJetFunc {nullptr};
  TF1 *m_rTrkFunc {nullptr};
 
  
  //int m_runNumber;
};

#endif  // JETBASE_JETCALIB_H
