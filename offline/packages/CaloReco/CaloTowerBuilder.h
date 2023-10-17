// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTOWERBUILDER_H
#define CALOTOWERBUILDER_H

#include <fun4all/SubsysReco.h>
#include "CaloWaveformProcessing.h"

#include <climits>
#include <string>

class CaloWaveformProcessing;
class PHCompositeNode;
class TowerInfoContainer;
class TowerInfoContainerv3;


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
    EPD = 3,
    MBD = 4,
    ZDC = 5
  };

  enum BuilderType
  {
    kPRDFTowerv1 = 0,
    kPRDFWaveform = 1,
    kWaveformTowerv2 = 2
  };

  void set_detector_type(CaloTowerBuilder::DetectorSystem dettype)
  {
    m_dettype = dettype;
    return;
  }

  void set_builder_type(CaloTowerBuilder::BuilderType buildertype)
  {
    m_buildertype = buildertype;
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

  void set_processing_type(CaloWaveformProcessing::process processingtype)
  {
    _processingtype = processingtype;
  }

  void set_softwarezerosuppression(bool usezerosuppression, int softwarezerosuppression)
  {
    _nsoftwarezerosuppression = softwarezerosuppression;
    _bdosoftwarezerosuppression = usezerosuppression;
  }

 private:
  CaloWaveformProcessing *WaveformProcessing;
  CaloTowerBuilder::DetectorSystem m_dettype;
  CaloTowerBuilder::BuilderType m_buildertype;
  TowerInfoContainer *m_CaloInfoContainer;  //! Calo info
  TowerInfoContainerv3 *m_CaloWaveformContainer; 
  std::string m_detector;
  int m_packet_low;
  int m_packet_high;
  int m_nsamples;
  int m_nchannels;
  int m_nzerosuppsamples;
  bool m_isdata;
  int _nsoftwarezerosuppression;
  bool _bdosoftwarezerosuppression;
  CaloWaveformProcessing::process _processingtype;
};

#endif  // CALOTOWERBUILDER_H
