#ifndef TPC_LASERCLUSTERIZER_H
#define TPC_LASERCLUSTERIZER_H

#include <fun4all/SubsysReco.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
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
class LaserClusterContainer;
class LaserCluster;
class PHCompositeNode;
class TrkrHitSet;
class TrkrHitSetContainer;
class PHG4TpcCylinderGeom;
class PHG4TpcCylinderGeomContainer;

class LaserClusterizer : public SubsysReco
{
 public:
  LaserClusterizer(const std::string &name = "LaserClusterizer");
  ~LaserClusterizer() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  //void calc_cluster_parameter(std::vector<pointKeyLaser> &clusHits, std::multimap<unsigned int, std::pair<std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey>, std::array<int, 3>>> &adcMap, bool isLamination);
  //void remove_hits(std::vector<pointKeyLaser> &clusHits, boost::geometry::index::rtree<pointKeyLaser, boost::geometry::index::quadratic<16>> &rtree, std::multimap<unsigned int, std::pair<std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey>, std::array<int, 3>>> &adcMap);

  void set_adc_threshold(float val) { m_adc_threshold = val; }
  void set_min_clus_size(float val) { min_clus_size = val; }
  void set_min_adc_sum(float val) { min_adc_sum = val; }
  void set_max_time_samples(int val) { m_time_samples_max = val; }
  void set_lamination(bool val) { m_lamination = val; }
  void set_do_sequential(bool val) { m_do_sequential = val; }
  void set_do_fitting(bool val) { m_do_fitting = val; }

 private:
  int m_event {-1};
  int m_time_samples_max {360};

  EventHeader *eventHeader{nullptr};

  LaserEventInfo *m_laserEventInfo {nullptr};

  TrkrHitSetContainer *m_hits {nullptr};
  LaserClusterContainer *m_clusterlist {nullptr};
  ActsGeometry *m_tGeometry {nullptr};
  PHG4TpcCylinderGeomContainer *m_geom_container {nullptr};
  double min_clus_size {1};
  double min_adc_sum {10};
  double m_adc_threshold {74.4};

  bool m_lamination {false};

  bool m_do_sequential {false};

  bool m_do_fitting {true};
  
  double m_tdriftmax {0};
  double AdcClockPeriod {53.0};  // ns
  double NZBinsSide {249};

  // TPC shaping offset correction parameter
  // From Tony Frawley July 5, 2022
  // double m_sampa_tbias {39.6};  // ns

};

#endif
