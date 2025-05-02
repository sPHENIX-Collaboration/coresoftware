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
  typedef std::pair<unsigned int, unsigned int> pixel;

  MvtxClusterPruner(const std::string &name = "MvtxClusterPruner");
  ~MvtxClusterPruner() override = default;

  //! module initialization
  int Init(PHCompositeNode * /*topNode*/) override { return 0; }

  //! run initialization
  int InitRun(PHCompositeNode * /*topNode*/) override;

  //! event processing
  int process_event(PHCompositeNode * /*topNode*/) override;

  //! end of process
  int End(PHCompositeNode * /*topNode*/) override;

};

#endif  // MVTX_MVTXHITPRUNER_H
