/**
 * @file mvtx/MvtxClusterizer.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Clusterizer for the MVTX
 */
#ifndef MVTX_MVTXCLUSTERIZER_H
#define MVTX_MVTXCLUSTERIZER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>

#include <string>                // for string
#include <utility>
class PHCompositeNode;
class TrkrHit;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;

/**
 * @brief Clusterizer for the MVTX
 */
class MvtxClusterizer : public SubsysReco
{
 public:
  typedef std::pair<unsigned int, unsigned int> pixel;

  MvtxClusterizer(const std::string &name = "MvtxClusterizer");
  ~MvtxClusterizer() override {}

  //! module initialization
  int Init(PHCompositeNode *topNode) override { return 0; }

  //! run initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  //! end of process
  int End(PHCompositeNode *topNode) override { return 0; }

  //! option to turn off z-dimension clustering
  void SetZClustering(const bool make_z_clustering)
  {
    m_makeZClustering = make_z_clustering;
  }
  bool GetZClustering() const
  {
    return m_makeZClustering;
  }

 private:
  //bool are_adjacent(const pixel lhs, const pixel rhs);
  bool are_adjacent(const std::pair<TrkrDefs::hitkey, TrkrHit*> &lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*> &rhs);

  void ClusterMvtx(PHCompositeNode *topNode);
  void PrintClusters(PHCompositeNode *topNode);

  // node tree storage pointers
  TrkrHitSetContainer *m_hits;
  TrkrClusterContainer *m_clusterlist; 

  TrkrClusterHitAssoc *m_clusterhitassoc;



  // settings
  bool m_makeZClustering;  // z_clustering_option
};

#endif  // MVTX_MVTXCLUSTERIZER_H
