#ifndef TPC_LASERCLUSTERIZER_H
#define TPC_LASERCLUSTERIZER_H

#include <fun4all/SubsysReco.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>

#include <phool/PHTimer.h>

#include <TFile.h>
#include <TH1I.h>
#include <TTree.h>

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
class RawHitSet;
class RawHitSetContainer;
class PHG4TpcCylinderGeom;
class PHG4TpcCylinderGeomContainer;

class LaserClusterizer : public SubsysReco
{
 public:
  typedef boost::geometry::model::point<float, 3, boost::geometry::cs::cartesian> point;
  typedef boost::geometry::model::box<point> box;
  typedef std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey> specHitKey;
  typedef std::pair<point, specHitKey> pointKeyLaser;

  LaserClusterizer(const std::string &name = "LaserClusterizer");
  ~LaserClusterizer() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  // void calc_cluster_parameter(std::vector<pointKeyLaser> &clusHits, std::multimap<unsigned int,std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>> &adcMap);
  void calc_cluster_parameter(std::vector<pointKeyLaser> &clusHits, std::multimap<unsigned int, std::pair<std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey>, std::array<int, 3>>> &adcMap, bool isLamination);
  // void remove_hits(std::vector<pointKeyLaser> &clusHits,  boost::geometry::index::rtree<pointKeyLaser, boost::geometry::index::quadratic<16> > &rtree, std::multimap <unsigned int, std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>> &adcMap, std::multimap <unsigned int, float*> &adcCoords);
  void remove_hits(std::vector<pointKeyLaser> &clusHits, boost::geometry::index::rtree<pointKeyLaser, boost::geometry::index::quadratic<16>> &rtree, std::multimap<unsigned int, std::pair<std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey>, std::array<int, 3>>> &adcMap);

  void set_debug(bool debug) { m_debug = debug; }
  void set_debug_name(const std::string &name) { m_debugFileName = name; }

  void set_adc_threshold(float val) { m_adc_threshold = val; }
  void set_min_clus_size(float val) { min_clus_size = val; }
  void set_min_adc_sum(float val) { min_adc_sum = val; }
  void set_max_time_samples(int val) { m_time_samples_max = val; }

 private:
  int m_event {-1};
  int m_time_samples_max {360};

  EventHeader *eventHeader{nullptr};

  LaserEventInfo *m_laserEventInfo {nullptr};

  TrkrHitSetContainer *m_hits {nullptr};
  RawHitSetContainer *m_rawhits {nullptr};
  LaserClusterContainer *m_clusterlist {nullptr};
  LaserClusterContainer *m_clusterlistLaminations {nullptr};
  ActsGeometry *m_tGeometry {nullptr};
  PHG4TpcCylinderGeomContainer *m_geom_container {nullptr};
  double min_clus_size {1};
  double min_adc_sum {10};
  double m_adc_threshold {74.4};

  double m_tdriftmax {0};
  double AdcClockPeriod {53.0};  // ns
  double NZBinsSide {249};

  bool do_read_raw {false};

  // TPC shaping offset correction parameter
  // From Tony Frawley July 5, 2022
  // double m_sampa_tbias {39.6};  // ns

  bool m_debug {false};
  std::string m_debugFileName {"LaserClusterizer_debug.root"};
  TFile *m_debugFile {nullptr};
  TTree *m_clusterTree {nullptr};
  // TTree *m_hitTree {nullptr};
  TH1I *m_itHist_0 {nullptr};
  TH1I *m_itHist_1 {nullptr};

  int m_nClus {0};
  double time_search {0};
  double time_clus {0};
  double time_erase {0};
  double time_all {0};

  LaserCluster *m_currentCluster {nullptr};
  LaserClusterContainer *m_eventClusters {nullptr};
  std::vector<float> m_currentHit;
  std::vector<float> m_currentHit_hardware;

  std::unique_ptr<PHTimer> t_all;
  std::unique_ptr<PHTimer> t_search;
  std::unique_ptr<PHTimer> t_clus;
  std::unique_ptr<PHTimer> t_erase;
};

#endif
