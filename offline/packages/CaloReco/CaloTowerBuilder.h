// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTOWERBUILDER_H
#define CALOTOWERBUILDER_H

#include <fun4all/SubsysReco.h>

#include <climits>
#include <string>

class CaloWaveformProcessing;
class PHCompositeNode;
class TowerInfoContainer;

class CaloTowerBuilder : public SubsysReco
{
 public:
  explicit CaloTowerBuilder(const std::string &name = "CaloTowerBuilder");
  ~CaloTowerBuilder() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void CreateNodeTree(PHCompositeNode *topNode);

  enum DetectorSystem
  {
    CEMC = 0,
    HCALIN = 1,
    HCALOUT = 2,
    EPD = 3
  };

  void set_detector_type(CaloTowerBuilder::DetectorSystem dettype)
  {
    m_dettype = dettype;
    return;
  }

  void set_nsamples(int _nsamples)
  {
    m_nsamples = _nsamples;
    return;
  }
  void set_dataflag(bool flag)
  {
    m_isdata = flag;
    return;
  }

 private:
  CaloWaveformProcessing *WaveformProcessing = nullptr;
  CaloTowerBuilder::DetectorSystem m_dettype = CaloTowerBuilder::CEMC;
  TowerInfoContainer *m_CaloInfoContainer = nullptr;  //! Calo info
  std::string m_detector = "CEMC";
  int m_packet_low = INT_MIN;
  int m_packet_high = INT_MIN;
  int m_nsamples = 16;
  int m_nchannels = 192;
  int m_nzerosuppsamples = 2;
  bool m_isdata = true;
};

#endif  // CALOTOWERBUILDER_H
