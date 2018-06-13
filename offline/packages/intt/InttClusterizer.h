/**
 * @file intt/InttClusterizer.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Clusterizer for the INTT
 */
#ifndef INTT_INTTCLUSTERIZER_H
#define INTT_INTTCLUSTERIZER_H

#include <fun4all/SubsysReco.h>
#include <utility>

class TrkrHitSetContainer;
class TrkrClusterContainer;

/**
 * @brief Clusterizer for the INTT
 */
class InttClusterizer : public SubsysReco {

public:

  typedef std::pair<unsigned int, unsigned int> pixel;

  InttClusterizer(const std::string &name = "InttClusterizer");
  virtual ~InttClusterizer(){}
  
  //! module initialization
  int Init(PHCompositeNode *topNode){return 0;}
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
  //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode){return 0;}
  
  //! option to turn off z-dimension clustering
  void SetZClustering(const bool make_z_clustering) {
    m_makeZClustering = make_z_clustering;
  }
  bool GetZClustering() const {
      return m_makeZClustering;
  }

private:

  bool are_adjacent(const pixel lhs, const pixel rhs );

  void ClusterIntt(PHCompositeNode *topNode);

  void PrintClusters(PHCompositeNode *topNode);
  
  // node tree storage pointers
  TrkrHitSetContainer* m_hits;
  TrkrClusterContainer* m_clusterlist;

  // settings
  bool m_makeZClustering;    // z_clustering_option
};

#endif // INTT_INTTCLUSTERIZER_H
