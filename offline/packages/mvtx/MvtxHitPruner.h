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

#include <string>                // for string
#include <utility>
class PHCompositeNode;
class TrkrHit;
class TrkrHitSetContainer;

/**
 * @brief Clusterizer for the MVTX
 */
class MvtxHitPruner : public SubsysReco
{
 public:
  typedef std::pair<unsigned int, unsigned int> pixel;

  MvtxHitPruner(const std::string &name = "MvtxHitPruner");
  ~MvtxHitPruner() override {}

  //! module initialization
  int Init(PHCompositeNode */*topNode*/) override { return 0; }

  //! run initialization
  int InitRun(PHCompositeNode * /*topNode*/) override;

  //! event processing
  int process_event(PHCompositeNode * /*topNode*/) override;

  //! end of process
  int End(PHCompositeNode */*topNode*/) override { return 0; }

 private:

  // node tree storage pointers
  TrkrHitSetContainer *m_hits;

  // settings

};

#endif  // MVTX_MVTXHITPRUNER_H
