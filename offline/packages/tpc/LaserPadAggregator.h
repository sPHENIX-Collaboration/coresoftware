#ifndef TPC_LASERPADAGGREGATOR_H
#define TPC_LASERPADAGGREGATOR_H

#include <fun4all/SubsysReco.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>

#include <map>
#include <string>
#include <vector>

class EventHeader;
class LaserEventInfo;
class LaserAggregatedPadContainer;
class LaserAggregatedPad;
class PHCompositeNode;
class TrkrHitSet;
class TrkrHitSetContainer;

class LaserPadAggregator : public SubsysReco
{
 public:
  explicit LaserPadAggregator(const std::string &name = "LaserPadAggregator");

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_adc_threshold(double val) { m_adc_threshold = val; }
  void set_padContainerNodeName(std::string &nodeName) { m_padContainerNodeName = nodeName; }

 private:
  int m_event {-1};

  EventHeader *eventHeader{nullptr};

  LaserEventInfo *m_laserEventInfo {nullptr};

  TrkrHitSetContainer *m_hits {nullptr};
  LaserAggregatedPadContainer *m_padlist {nullptr};
  double m_adc_threshold {74.4};

  unsigned int nLaserEvents {0};

  std::string m_padContainerNodeName {"LASER_AGGREGATED_PAD_CONTAINER"};

  std::map<std::pair<uint8_t, uint16_t>, double> m_adcMap[2][12][3];
  std::map<std::pair<uint8_t, uint16_t>, double> m_hitMap[2][12][3];

};

#endif
