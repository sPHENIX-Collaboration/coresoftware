/**
 * @file mvtx/MvtxClusterPruner.h
 * @author Hugo Pereira Da Costa <hugo.pereira-da-costa@lanl.gov>
 * @date May 2025
 * @brief Pruner to handle multiple copies of clusters for the MVTX
 */
#ifndef MVTX_MVTXCLUSTERPRUNER_H
#define MVTX_MVTXCLUSTERPRUNER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>  // for string
#include <utility>

// forward declarations
class PHCompositeNode;


/**
 * @brief cluster pruner for the MVTX
 */
class MvtxClusterPruner : public SubsysReco
{
  public:

  //! constructor
  MvtxClusterPruner(const std::string &name = "MvtxClusterPruner");

  //! module initialization
  int Init(PHCompositeNode * /*topNode*/) override { return 0; }

  //! run initialization
  int InitRun(PHCompositeNode * /*topNode*/) override;

  //! event processing
  int process_event(PHCompositeNode * /*topNode*/) override;

  //! end of process
  int End(PHCompositeNode * /*topNode*/) override;

  //! strict matching
  void set_use_strict_matching( bool value )
  { m_use_strict_matching = value; }

  private:

  //! strict matching
  bool m_use_strict_matching = false;

  //!@name counters
  //@{
  uint64_t m_cluster_counter_total = 0;
  uint64_t m_cluster_counter_deleted = 0;
  //@}
};

#endif  // MVTX_MVTXHITPRUNER_H
