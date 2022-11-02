#include "TpcClusterBuilder.h"
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TpcDefs.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <algorithm>

class TpcClusterBuilder;
using std::cout, std::endl;

TpcClusterBuilder& TpcClusterBuilder::operator+=(const TpcClusterBuilder& rhs) {
  // layer, side, and layerGeom won't change between different additions, but they
  // but need to be set by the first addition
  if (rhs.neff_electrons==0) return *this; 

  layerGeom = rhs.layerGeom;
  layer     = rhs.layer;
  side      = rhs.side;

  if (rhs.hasPhiBins) {
    if (!hasPhiBins) {
      phi_bin_lo = rhs.phi_bin_lo;
      phi_bin_hi = rhs.phi_bin_hi;
    } else {
      int rhs_phi_bin_lo = rhs.phi_bin_lo;
      int rhs_phi_bin_hi = rhs.phi_bin_hi;

      // if this is the case, wrap the lower values to be higher
      if ( (phi_bin_lo - rhs_phi_bin_lo) > (nphibins/2)) {
        rhs_phi_bin_lo += nphibins;
        rhs_phi_bin_hi += nphibins;
      } else if ( (rhs_phi_bin_lo - phi_bin_lo) > (nphibins/2)) {
        phi_bin_lo += nphibins;
        phi_bin_hi += nphibins;
      }

      if (rhs_phi_bin_lo < phi_bin_lo) phi_bin_lo = rhs_phi_bin_lo;
      if (rhs_phi_bin_hi > phi_bin_hi) phi_bin_hi = rhs_phi_bin_hi;
    }
    hasPhiBins = true;
  } 


  if (rhs.hasTimeBins) {
    if (!hasTimeBins) {
      time_bin_lo = rhs.time_bin_lo;
      time_bin_hi = rhs.time_bin_hi;
    } else {
      if (rhs.time_bin_lo < time_bin_lo) time_bin_lo = rhs.time_bin_lo;
      if (rhs.time_bin_hi > time_bin_hi) time_bin_hi = rhs.time_bin_hi;
    }
    hasTimeBins = true;
  }

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
  phi_bin_lo     = INT_MAX;
  phi_bin_hi     = INT_MIN;
  time_bin_lo    = INT_MAX;
  time_bin_hi    = INT_MIN;
  nphibins       = 0;
  hasPhiBins     = false;
  hasTimeBins    = false;
}

TpcClusterBuilder::PairCluskeyCluster TpcClusterBuilder::build(MapHitsetkeyUInt& cluster_cnt) const {
/* std::pair<TrkrDefs::cluskey, TrkrClusterv4*> TpcClusterBuilder::build_cluster( */
    /* std::map<TrkrDefs::hitsetkey, unsigned int>& cluster_cnt) const { */
  auto phi_mean { phi_integral / neff_electrons };
  TrkrClusterv4* cluster = new TrkrClusterv4();
  cluster->setPosition ( 0, phi_mean );
  cluster->setPosition ( 1, time_integral / neff_electrons);
  cluster->setPhiSize  ( static_cast<char>(phi_bin_hi-phi_bin_lo+1) );
  if (phi_bin_hi < phi_bin_lo) std::cout << " WARNING 101: phi_bin_hi<phi_bin_lo: " 
     << phi_bin_hi << "<" << phi_bin_lo << std::endl;
  cluster->setZSize    ( static_cast<char>(time_bin_hi-time_bin_lo+1));
  cluster->setAdc      ( neff_electrons );

  const int phi_pad_number           = layerGeom->get_phibin(phi_mean);
  const auto phibins                 = layerGeom->get_phibins();
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

void TpcClusterBuilder::fillPhiBins(const std::vector<int>& bins) {
  if (bins.size() == 0) return;
  phi_bin_lo = bins[0];
  phi_bin_hi = bins[bins.size()-1];
  // Check if the phi_bin_hi "wrapped around" to be below phi_bin_lo; 
  // If it did wrap it "back above" phi_bin_lo
  if (phi_bin_hi < phi_bin_lo) phi_bin_hi += nphibins;
  hasPhiBins = true;
}

void TpcClusterBuilder::fillTimeBins(const std::vector<int>& bins) {
  if (bins.size() == 0) return;
  time_bin_lo = bins[0];
  time_bin_hi = bins[bins.size()-1];
  hasTimeBins = true;
}
