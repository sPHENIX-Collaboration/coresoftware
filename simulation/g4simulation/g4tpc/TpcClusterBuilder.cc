#include "TpcClusterBuilder.h"
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TpcDefs.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>

class TpcClusterBuilder;

TpcClusterBuilder& TpcClusterBuilder::operator+=(const TpcClusterBuilder& rhs) {
  // layer, side, and layerGeom won't change between different additions, but they
  // but need to be set by the first addition
  layerGeom = rhs.layerGeom;
  layer     = rhs.layer;
  side      = rhs.side;


  if (rhs.phi_bin_lo < phi_bin_lo) phi_bin_lo = rhs.phi_bin_lo;
  if (rhs.phi_bin_hi > phi_bin_hi) phi_bin_hi = rhs.phi_bin_hi;

  if (rhs.time_bin_lo < time_bin_lo) time_bin_lo = rhs.time_bin_lo;
  if (rhs.time_bin_hi > time_bin_hi) time_bin_hi = rhs.time_bin_hi;

  neff_electrons += rhs.neff_electrons;
  phi_integral   += rhs.phi_integral;
  time_integral  += rhs.time_integral;
  return *this;
}

void TpcClusterBuilder::reset() {
  //hitsetkey = 0;
  layerGeom      = nullptr; 
  layer          = 0;
  side           = 0;
  neff_electrons = 0;
  phi_integral   = 0.;
  time_integral  = 0.;
  phi_bin_lo     = 0;
  phi_bin_hi     = 0;
  time_bin_lo    = 0;
  time_bin_hi    = 0;
}

TpcClusterBuilder::PairCluskeyCluster TpcClusterBuilder::build(MapHitsetkeyUInt& cluster_cnt) const {
/* std::pair<TrkrDefs::cluskey, TrkrClusterv4*> TpcClusterBuilder::build_cluster( */
    /* std::map<TrkrDefs::hitsetkey, unsigned int>& cluster_cnt) const { */
  auto phi_mean { phi_integral / neff_electrons };
  TrkrClusterv4* cluster = new TrkrClusterv4();
  cluster->setPosition ( 0, phi_mean );
  cluster->setPosition ( 1, time_integral / neff_electrons);
  cluster->setPhiSize  ( phi_bin_hi  - phi_bin_lo  +1);
  cluster->setZSize    ( time_bin_hi - time_bin_lo +1);
  cluster->setAdc      ( neff_electrons );

  // generate the hitset key at the mean phi location:
  /* std::cout << " a0 " << std::endl; */
  const int phi_pad_number           = layerGeom->get_phibin(phi_mean);
  /* std::cout << " a1 " << std::endl; */
  const auto phibins                 = layerGeom->get_phibins();
  /* std::cout << " a2 " << std::endl; */
  const unsigned int pads_per_sector = phibins / 12;
  const unsigned int sector          = phi_pad_number / pads_per_sector;
  TrkrDefs::hitsetkey hitsetkey = TpcDefs::genHitSetKey(layer, sector, side);

  // generate the clusterkey 
  unsigned int which_cluster = 0;
  auto iter = cluster_cnt.find( hitsetkey );
  if (iter == cluster_cnt.end()) {
    cluster_cnt[hitsetkey] = which_cluster;
  } else {
    iter->second += 1;
    which_cluster = iter->second;
  }

  // return the clusterkey which is needed to index the track
  auto cluskey = TrkrDefs::genClusKey( hitsetkey, which_cluster );
  return std::make_pair( cluskey, cluster );
}

bool TpcClusterBuilder::has_data() const {
  return neff_electrons != 0;
}

