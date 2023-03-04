// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPC_HITS_H
#define TPC_HITS_H

#include "TPCMap.h"
#include "TPC_RawHit.h"

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TrkrHitSetContainer;
class TrkrHitSet;
class TrkrHit;

class tpc_hits : public SubsysReco
{
 public:
  tpc_hits(const std::string &name = "tpc_hits");

  ~tpc_hits() override;

  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
  */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
  */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

 protected:
  static const int layercount = 16;
  static const int layeroffset = 7 + 16;

  TrkrHitSetContainer *m_hits = nullptr;
  TrkrHitSet *m_hitset[layercount] = {};
  TrkrHit *m_hit = nullptr;

  // RawHitSetContainer *m_rawhits __attribute__ ((unused)) = nullptr;

  TPCMap M;

  int starting_BCO;
  int rollover_value;
  int current_BCOBIN;
};

#endif  // TPC_HITS_H
