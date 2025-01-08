#ifndef TPC_TPC3DCLUSTERIZER_H
#define TPC_TPC3DCLUSTERIZER_H

#include <fun4all/SubsysReco.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>

#include <phool/PHTimer.h>

#include <TFile.h>
#include <TH1I.h>
#include <TTree.h>
#include <TNtuple.h>

// BOOST for combi seeding
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <map>
#include <string>
#include <vector>

class LaserClusterContainerv1;
class LaserClusterv1;
class PHCompositeNode;
class TrkrHitSet;
class TrkrHitSetContainer;
class PHG4TpcCylinderGeom;
class PHG4TpcCylinderGeomContainer;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<float, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey> specHitKey;
typedef std::pair<point, specHitKey> pointKeyLaser;

class Tpc3DClusterizer : public SubsysReco
{
 public:
  Tpc3DClusterizer(const std::string &name = "Tpc3DClusterizer");
  ~Tpc3DClusterizer() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  // void calc_cluster_parameter(std::vector<pointKeyLaser> &clusHits, std::multimap<unsigned int,std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>> &adcMap);
  void calc_cluster_parameter(std::vector<pointKeyLaser> &clusHits, std::multimap<unsigned int, std::pair<std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey>, std::array<int, 3>>> &adcMap);
  // void remove_hits(std::vector<pointKeyLaser> &clusHits,  bgi::rtree<pointKeyLaser, bgi::quadratic<16> > &rtree, std::multimap <unsigned int, std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>> &adcMap, std::multimap <unsigned int, float*> &adcCoords);
  void remove_hits(std::vector<pointKeyLaser> &clusHits, bgi::rtree<pointKeyLaser, bgi::quadratic<16>> &rtree, std::multimap<unsigned int, std::pair<std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey>, std::array<int, 3>>> &adcMap);

  void set_debug(bool debug) { m_debug = debug; }
  void set_debug_name(const std::string &name) { m_debugFileName = name; }
  void set_output(bool output) { m_output = output; }
  void set_output_name(const std::string &name) { m_outputFileName = name; }

  void set_pedestal(float val) { pedestal = val; }
  void set_min_clus_size(float val) { min_clus_size = val; }
  void set_min_adc_sum(float val) { min_adc_sum = val; }

 private:
  int m_event = -1;
  int m_seed = -1;

  TrkrHitSetContainer *m_hits = nullptr;
  LaserClusterContainerv1 *m_clusterlist = nullptr;
  ActsGeometry *m_tGeometry = nullptr;
  PHG4TpcCylinderGeomContainer *m_geom_container = nullptr;
  double pedestal = 74.4;
  double min_clus_size = 1;
  double min_adc_sum = 10;
  //  double m_pedestal = 74.4;

  double m_tdriftmax = 0;
  double AdcClockPeriod = 53.0;  // ns
  double NZBinsSide = 360-76;//249;


  // TPC shaping offset correction parameter
  // From Tony Frawley July 5, 2022
  // double m_sampa_tbias = 39.6;  // ns

  bool m_debug = false;
  bool m_output = false;
  std::string m_debugFileName = "Tpc3DClusterizer_debug.root";
  std::string m_outputFileName = "Tpc3DClusterizer_output.root";
  TFile *m_debugFile = nullptr;
  TFile *m_outputFile = nullptr;
  TTree *m_clusterTree = nullptr;
  TNtuple *m_clusterNT = nullptr;
  // TTree *m_hitTree = nullptr;
  TH1I *m_itHist_0 = nullptr;
  TH1I *m_itHist_1 = nullptr;

  TH1D *m_tHist_0 = nullptr;
  TH1D *m_tHist_1 = nullptr;

  int m_nClus = 0;
  double time_search = 0;
  double time_clus = 0;
  double time_erase = 0;
  double time_all = 0;

  LaserClusterv1 *m_currentCluster = nullptr;
  std::vector<LaserClusterv1 *> m_eventClusters = {nullptr};
  std::vector<float> m_currentHit;
  std::vector<float> m_currentHit_hardware;

  std::unique_ptr<PHTimer> t_all;
  std::unique_ptr<PHTimer> t_search;
  std::unique_ptr<PHTimer> t_clus;
  std::unique_ptr<PHTimer> t_erase;
};

#endif
