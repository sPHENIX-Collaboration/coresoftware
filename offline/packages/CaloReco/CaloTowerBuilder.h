// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTOWERBUILDER_H
#define CALOTOWERBUILDER_H

#include "CaloWaveformProcessing.h"
#include "CaloTowerDefs.h"

#include <fun4all/SubsysReco.h>

#include <limits>
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

  void set_detector_type(CaloTowerDefs::DetectorSystem dettype)
  {
    m_dettype = dettype;
    return;
  }

  void set_builder_type(CaloTowerDefs::BuilderType buildertype)
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
  CaloWaveformProcessing *WaveformProcessing {nullptr};
  TowerInfoContainer *m_CaloInfoContainer {nullptr};  //! Calo info
  bool m_isdata {true};
  bool _bdosoftwarezerosuppression {false};
  int m_packet_low {std::numeric_limits<int>::min()};
  int m_packet_high {std::numeric_limits<int>::min()};
  int m_nsamples {16};
  int m_nchannels {192};
  int m_nzerosuppsamples {2};
  int _nsoftwarezerosuppression {40};
  CaloTowerDefs::DetectorSystem m_dettype {CaloTowerDefs::CEMC};
  CaloTowerDefs::BuilderType m_buildertype {CaloTowerDefs::kPRDFTowerv1};
  CaloWaveformProcessing::process _processingtype {CaloWaveformProcessing::NONE};
  std::string m_detector {"CEMC"};
};

#endif  // CALOTOWERBUILDER_H
