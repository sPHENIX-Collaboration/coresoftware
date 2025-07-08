// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_TRACKING_TPCSILICONQA_H
#define QA_TRACKING_TPCSILICONQA_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>
#include <vector>

class PHCompositeNode;
class TH1;
class TH2;

class TpcSiliconQA : public SubsysReco
{
 public:
  TpcSiliconQA(const std::string &name = "TpcSiliconQA");

  ~TpcSiliconQA() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;

  void setXCut(float cutVal) { m_xcut = cutVal; }
  void setYCut(float cutVal) { m_ycut = cutVal; }
  void setEtaCut(float cutVal) { m_etacut = cutVal; }
  void setPhiCut(float cutVal) { m_phicut = cutVal; }

 private:
  void createHistos();

  float m_crossing = std::numeric_limits<float>::quiet_NaN();
  float m_silseedx = std::numeric_limits<float>::quiet_NaN();
  float m_silseedy = std::numeric_limits<float>::quiet_NaN();
  float m_silseedz = std::numeric_limits<float>::quiet_NaN();
  float m_silseedphi = std::numeric_limits<float>::quiet_NaN();
  float m_silseedeta = std::numeric_limits<float>::quiet_NaN();
  float m_tpcseedx = std::numeric_limits<float>::quiet_NaN();
  float m_tpcseedy = std::numeric_limits<float>::quiet_NaN();
  float m_tpcseedz = std::numeric_limits<float>::quiet_NaN();
  float m_tpcseedphi = std::numeric_limits<float>::quiet_NaN();
  float m_tpcseedeta = std::numeric_limits<float>::quiet_NaN();

  float m_xcut = 1.0;
  float m_ycut = 1.0;
  float m_etacut = 0.05;
  float m_phicut = 0.25;

  std::string getHistoPrefix() const;
  int m_event = 0;

  TH1 *h_crossing = nullptr;
  TH1 *h_phiDiff[8] = {nullptr};
  TH1 *h_etaDiff[8] = {nullptr};
  TH1 *h_xDiff[8] = {nullptr};
  TH1 *h_yDiff[8] = {nullptr};
  TH1 *h_zDiff[8] = {nullptr};
};

#endif  // QA_TRACKING_TPCSILICONQA_H
