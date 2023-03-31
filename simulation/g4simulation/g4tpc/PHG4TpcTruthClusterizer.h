#ifndef G4TPC_TRUTHCLUSTERIZER__H
#define G4TPC_TRUTHCLUSTERIZER__H

#include <g4tracking/TruthClusterizerBase.h>
  // This is the helper class for PHG4TpcElectronDrift. It will cluster the
  // PHG4Hits which come from each individual embedded ("Truth") track. It uses
  // pointers to the PHG4TpcDigitizer and TpcClusterizer modules , which also
  // do the clustering for all the PHG4Hits together (from all truth tracks
  // together + all other tracks + noise).
  //
  // Because it makes calls to these other modules, t needs to be initialized
  // from the Fun4All driving macro *after* adding the above listed modules.

class PHG4TpcDigitizer;
class TpcClusterizer;
class PHCompositeNode;

class  PHG4TpcTruthClusterizer : public TruthClusterizerBase {
  private:
  PHCompositeNode*   m_topNode     { nullptr };
  PHG4TpcDigitizer*  m_digitiser      ;
  TpcClusterizer*    m_clusterizer    ;

  public:
  PHG4TpcTruthClusterizer ();

  void init_run(PHCompositeNode*& _topNode, int _verbosity);
  int clusterize_hits(TrkrClusterContainer*);
  void check_g4hit(PHG4Hit*);
  void end_of_event();
  ~PHG4TpcTruthClusterizer(){};

  // all data below is implemented for the clusterize_hits() function
  // It does the work that the following modules do for SvtxTrack's:
  // (1) g4tpc/PHG4TpcDigitizer -- use prefix _D__
  // (2) tpc/TpcClusterizer     -- use prefix _C__ // note that this can probably be *greatly* simplified
};

#endif
