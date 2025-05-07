/**
 * @file mvtx/MvtxHitPruner.h
 * @author Tony Frawley
 * @date July 2022
 * @brief Pruner to handle multiple copies of hits for the MVTX
 */
#ifndef MVTX_MVTXHITPRUNER_H
#define MVTX_MVTXHITPRUNER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>  // for string
#include <utility>
class PHCompositeNode;
class TrkrHit;
class TrkrHitSetContainer;

/**
 * @brief Hit pruner for the MVTX
 */
class MvtxHitPruner : public SubsysReco
{
 public:

  //! constructor
  MvtxHitPruner(const std::string &name = "MvtxHitPruner");

  //! destructor
  ~MvtxHitPruner() override = default;

  //! module initialization
  int Init(PHCompositeNode * /*topNode*/) override { return 0; }

  //! run initialization
  int InitRun(PHCompositeNode * /*topNode*/) override;

  //! event processing
  int process_event(PHCompositeNode * /*topNode*/) override;

  //! end of process
  int End(PHCompositeNode * /*topNode*/) override { return 0; }

 private:
  // node tree storage pointers
  TrkrHitSetContainer *m_hits = nullptr;
};

#endif  // MVTX_MVTXHITPRUNER_H
