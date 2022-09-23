/**
 * @file g4tpc/TrkrTruthTrackv1.cc
 * @author D. Stewart
 * @date September 2022
 * @brief Version 1 of TrkrTruthTrack
 */
#include "TrkrTruthTrackv1.h"

#include <cmath>
#include <TLorentzVector.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>

void TrkrTruthTrackv1::identify(std::ostream &os) const 
{
  os << " TrkrTruthTrack: " << std::endl
     << "   trackid(" << trackid << ")  [X0,Y0,Z0]("<<X0<<","<<Y0<<","<<Z0<<")  "
     << " [pseudorapidity,pt,phi]("<<pseudoRapidity<<","<<pt<<","<<phi<<")" << std::endl;
  os << " Clusters HitSetKey(layer) : " << std::endl << "  ";
  int cnt = 0;
  for (auto cluster : clusters) {
    if (cnt == 8) {
      cnt = 0;
      os << std::endl << "  ";
    }
    if (cnt > 0) os << ", ";
    uint32_t i_hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluster);
    int layer = TrkrDefs::getLayer(cluster);
    os << " " << i_hitsetkey<<"("<<layer<< ")";
    ++cnt;
  }
  if (cnt != 0) os << std::endl;
}

TrkrTruthTrackv1::TrkrTruthTrackv1() 
  : X0 {0.}
  , Y0 {0.}
  , Z0 {0.}
  , pseudoRapidity {0.}
  , pt {0.}
  , phi {0.}
  , clusters {}
{ trackid = 0; }

TrkrTruthTrackv1::TrkrTruthTrackv1(unsigned int _trackid, PHG4Particle* p, PHG4VtxPoint* vtx) 
  : clusters {}
{
  trackid = _trackid;

  X0 = vtx->get_x();
  Y0 = vtx->get_y();
  Z0 = vtx->get_z();

  TLorentzVector v1;
  v1.SetPxPyPzE(p->get_px(), p->get_py(), p->get_pz(), p->get_e());
  phi = v1.Phi();
  pseudoRapidity = v1.PseudoRapidity();
  pt = v1.Pt();
}

void TrkrTruthTrackv1::addCluster(TrkrDefs::cluskey key) {
  clusters.push_back(key);
}
