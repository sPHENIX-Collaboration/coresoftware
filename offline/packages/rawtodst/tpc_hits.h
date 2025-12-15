// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPC_HITS_H
#define TPC_HITS_H

#include "TPCMap.h"
#include "TPC_RawHit.h"

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>

class PHCompositeNode;
class Fun4AllHistoManager;
class TH2;

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

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

 private:
  static const int layercount {16};
  static const int layeroffset {7 + 16};
  Fun4AllHistoManager *hm {nullptr};
  TH2 *_h_hit_XY {nullptr};
  TrkrHitSet *m_hitset[layercount] = {nullptr};

  int starting_BCO {-1};
  int rollover_value {0};
  int current_BCOBIN {0};

  TPCMap M;

  std::string _filename {"./outputfile.root"};
};

#endif  // TPC_HITS_H
