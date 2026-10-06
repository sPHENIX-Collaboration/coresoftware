#ifndef TPC_LASERAGGREGATEDCLUSTERIZER_H
#define TPC_LASERAGGREGATEDCLUSTERIZER_H

#include <fun4all/SubsysReco.h>

#include <trackbase/ActsGeometry.h>
#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <pthread.h>

#include <map>
#include <string>
#include <vector>

class ActsGeometry;
class LaserAggregatedPad;
class LaserAggregatedPadContainer;
class LaserClusterContainer;
class LaserCluster;
class PHCompositeNode;


class LaserAggregatedClusterizer : public SubsysReco
{
 public:
  explicit LaserAggregatedClusterizer(const std::string &name = "LaserAggregatedClusterizer");

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void set_padContainerNodeName(std::string &nodeName) { m_padContainerNodeName = nodeName; }
  void set_clusterNodeName(std::string &nodeName) { m_clusterNodeName = nodeName; }
  void set_nHitPerLaserEventMin(double nHitPerLaserEventMin) { m_nHitPerLaserEventMin = nHitPerLaserEventMin; }
  void set_QAFile(std::string &name) { m_QAName = name; }

  void AddSegmentFile(const std::string &f) { m_segfiles.push_back(f); }
  void SetSegmentList(const std::string &listfile);   // parse a text file of paths

 private:  
  pthread_mutex_t m_threadlock;
  LaserClusterContainer *m_clusterlist {nullptr};
  double m_nHitPerLaserEventMin {1.0};

  std::string m_padContainerNodeName {"LASER_AGGREGATED_PAD_CONTAINER"};
  std::string m_clusterNodeName {"LASER_AGGREGATED_CLUSTER"};

  std::vector<std::string> m_segfiles;
  LaserAggregatedPadContainer *m_aggregatedPads {nullptr};

  bool m_done {false};

  std::string m_QAName{""};


  // TPC shaping offset correction parameter
  // From Tony Frawley July 5, 2022
  // double m_sampa_tbias {39.6};  // ns

};

#endif
