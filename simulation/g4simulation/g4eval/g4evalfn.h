#ifndef G4EVALFN__H
#define G4EVALFN__H

#include "TrkrClusLoc.h"

class TrkrClusterIsMatcher;
class EmbRecoMatchContainer;
class SvtxTrackMap;

namespace g4evalfn {

  enum DET { MVTX=0, INTT=1, TPC=2, TPOT=3 }; // 

  int trklayer_det(TrkrDefs::hitsetkey); // 0:MVTX 1:INTt 2:TPC 3:TPOT and beyond
  int trklayer_det(TrkrDefs::cluskey);   
  int trklayer_det(int layer);          

  TrkrClusLoc clusloc_PHG4(TrkrClusterIsMatcher*, TrkrDefs::cluskey);
  TrkrClusLoc clusloc_SVTX(TrkrClusterIsMatcher*, TrkrDefs::cluskey);

  inline float abs_dphi (float aphi, float bphi) {
    float phi_delta = fabs(aphi-bphi);
    while (phi_delta > M_PI) phi_delta = fabs(phi_delta-2*M_PI);
    return phi_delta;
  }

  std::vector<int> unmatchedSvtxTrkIds(EmbRecoMatchContainer*, SvtxTrackMap*);

  float calc_match_statistic(TrkrClusterIsMatcher* ismatcher, TrkrDefs::cluskey key_A, TrkrDefs::cluskey key_B);
}

#endif
