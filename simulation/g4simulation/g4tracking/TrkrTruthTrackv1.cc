/**
 * @file g4tpc/TrkrTruthTrackv1.cc
 * @author D. Stewart
 * @date September 2022
 * @brief Version 1 of TrkrTruthTrack
 */
#include "TrkrTruthTrackv1.h"

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>

#include <TLorentzVector.h>

#include <cmath>

void TrkrTruthTrackv1::identify(std::ostream& os) const
{
  os << " TrkrTruthTrack: " << std::endl
     << "   trackid(" << trackid << ")  [X0,Y0,Z0](" << X0 << "," << Y0 << "," << Z0 << ")  "
     << " [pseudorapidity,pt,phi](" << pseudoRapidity << "," << pt << "," << phi << ")" << std::endl;
  os << " Clusters HitSetKey(layer) : " << std::endl
     << "  ";
  int cnt = 0;
  for (auto cluster : clusters)
  {
    if (cnt == 8)
    {
      cnt = 0;
      os << std::endl
         << "  ";
    }
    if (cnt > 0)
    {
      os << ", ";
    }
    uint32_t i_hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluster);
    int layer = TrkrDefs::getLayer(cluster);
    os << " " << i_hitsetkey << "(" << layer << ")";
    ++cnt;
  }
  if (cnt != 0)
  {
    os << std::endl;
  }
}

TrkrTruthTrackv1::TrkrTruthTrackv1()
  : trackid{UINT_MAX}
  , X0{NAN}
  , Y0{NAN}
  , Z0{NAN}
  , pseudoRapidity{NAN}
  , pt{NAN}
  , phi{NAN}
  , clusters{}
{
}

TrkrTruthTrackv1::TrkrTruthTrackv1(unsigned int _trackid, PHG4Particle* p, PHG4VtxPoint* vtx)
  : trackid{_trackid}
  , clusters{}
{
  X0 = vtx->get_x();
  Y0 = vtx->get_y();
  Z0 = vtx->get_z();

  TLorentzVector v1;
  v1.SetPxPyPzE(p->get_px(), p->get_py(), p->get_pz(), p->get_e());
  phi = v1.Phi();
  pseudoRapidity = v1.PseudoRapidity();
  pt = v1.Pt();
}

void TrkrTruthTrackv1::addCluster(TrkrDefs::cluskey key)
{
  clusters.push_back(key);
}


bool TrkrTruthTrackv1::has_hitsetkey(TrkrDefs::hitsetkey key) const {
  return std::binary_search(clusters.begin(), clusters.end(), key, CompHitSetKey() );
}

bool TrkrTruthTrackv1::has_hitsetkey(TrkrDefs::cluskey key) const {
  return std::binary_search(clusters.begin(), clusters.end(), key, CompHitSetKey() );
}

std::pair<bool, TrkrDefs::cluskey> TrkrTruthTrackv1::get_cluskey(TrkrDefs::hitsetkey hitsetkey) const { 
  auto lb = std::lower_bound(clusters.begin(), clusters.end(), hitsetkey, CompHitSetKey());
  if (lb == clusters.end() || TrkrDefs::getHitSetKeyFromClusKey(*lb) != hitsetkey) return { false, 0. };
  else                      return { true,  *lb };
}
