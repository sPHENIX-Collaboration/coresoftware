#ifndef TPC_LASERAGGREGATOR_H
#define TPC_LASERAGGREGATOR_H

#include <fun4all/SubsysReco.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>


// BOOST for combi seeding
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

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
//class PHG4TpcGeom;
//class PHG4TpcGeomContainer;

class LaserPadAggregator : public SubsysReco
{
 public:
  LaserPadAggregator(const std::string &name = "LaserPadAggregator");
  ~LaserPadAggregator() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  //void calc_cluster_parameter(std::vector<pointKeyLaser> &clusHits, std::multimap<unsigned int, std::pair<std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey>, std::array<int, 3>>> &adcMap, bool isLamination);
  //void remove_hits(std::vector<pointKeyLaser> &clusHits, boost::geometry::index::rtree<pointKeyLaser, boost::geometry::index::quadratic<16>> &rtree, std::multimap<unsigned int, std::pair<std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey>, std::array<int, 3>>> &adcMap);

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

  // TPC shaping offset correction parameter
  // From Tony Frawley July 5, 2022
  // double m_sampa_tbias {39.6};  // ns

};

#endif
