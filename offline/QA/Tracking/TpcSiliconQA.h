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

 private:
  void createHistos();

  std::string m_trackMapName = "SvtxTrackMap"; 
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

  std::string getHistoPrefix() const;
  int m_event = 0;

  TH1 *h_crossing = nullptr;
  TH1 *h_trackMatch = nullptr;
  TH1 *h_phiDiff = nullptr;
  TH1 *h_etaDiff = nullptr;
  TH1 *h_xDiff = nullptr;
  TH1 *h_yDiff = nullptr;
  TH1 *h_zDiff = nullptr; 
};

#endif  // QA_TRACKING_TPCSILICONQA_H
