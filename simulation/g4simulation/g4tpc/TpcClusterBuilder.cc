#include "TpcClusterBuilder.h"
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TpcDefs.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <algorithm>
#include <Acts/Surfaces/Surface.hpp>

class TpcClusterBuilder;
using std::cout, std::endl;

TpcClusterBuilder& TpcClusterBuilder::operator+=(const TpcClusterBuilder& rhs) {
  // layer, side, and layerGeom won't change between different additions, but they
  // but need to be set by the first addition
  if (!rhs.has_data) return *this; 

  if (!has_data) {
    layer          = rhs.layer;
    side           = rhs.side;
    neff_electrons = rhs.neff_electrons;
    phi_integral   = rhs.phi_integral;
    time_integral  = rhs.time_integral;
    phi_bin_lo     = rhs.phi_bin_lo;
    phi_bin_hi     = rhs.phi_bin_hi;
    time_bin_lo    = rhs.time_bin_lo;
    time_bin_hi    = rhs.time_bin_hi;
    nphibins       = rhs.nphibins;
    hasPhiBins     = rhs.hasPhiBins;
    hasTimeBins    = rhs.hasTimeBins;
    layerGeom      = rhs.layerGeom;
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
    hasPhiBins = true;

    if (rhs.time_bin_lo < time_bin_lo) time_bin_lo = rhs.time_bin_lo;
    if (rhs.time_bin_hi > time_bin_hi) time_bin_hi = rhs.time_bin_hi;
    hasTimeBins = true;

    neff_electrons += rhs.neff_electrons;
    phi_integral   += rhs.phi_integral;
    time_integral  += rhs.time_integral;
  }
  set_has_data();
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
  has_data       = false;
}

TpcClusterBuilder::PairCluskeyCluster TpcClusterBuilder::build(MapHitsetkeyUInt& cluster_cnt) const {
  TrkrClusterv4* cluster = new TrkrClusterv4();

  int phi_size = phi_bin_hi-phi_bin_lo+1;
  if (phi_size > CHAR_MAX || phi_size < 0) return { 0, nullptr };
  cluster->setPhiSize  ( static_cast<char>(phi_size) );

  int Z_size = time_bin_hi-time_bin_lo+1;
  if (Z_size > CHAR_MAX || Z_size < 0)  return { 0, nullptr };
  cluster->setZSize    ( static_cast<char>(Z_size));

  cluster->setAdc      ( neff_electrons );
  cout << " FIXME a-2 " << endl;

  // copy logic from packages/tpc/TpcClusterizer.cc ~line 333
  const double clusphi { phi_integral  / neff_electrons };
  const double clust   { time_integral / neff_electrons };
  const int    phi_pad_number        = layerGeom->get_phibin(clusphi);
  const auto   phibins               = layerGeom->get_phibins();
  const unsigned int pads_per_sector = phibins / 12;
  const unsigned int sector          = phi_pad_number / pads_per_sector;
  TrkrDefs::hitsetkey hitsetkey = TpcDefs::genHitSetKey(layer, sector, side);
  
  const double radius = layerGeom->get_radius();  // returns center of layer
  const float  clusx  = radius * cos(clusphi);
  const float  clusy  = radius * sin(clusphi);

  cout << " FIXME a-2 " << endl;
  const unsigned short NTBins = (unsigned short) layerGeom->get_zbins();
  const double m_tdriftmax = AdcClockPeriod * NTBins / 2.0;  
  const double zdriftlength = clust * tGeometry->get_drift_velocity();
  const double z_sign = (side==0) ? -1 : 1;
  const double clusz  = z_sign * (m_tdriftmax * tGeometry->get_drift_velocity()-zdriftlength);

  cout << " FIXME a-1 " << endl;
  Acts::Vector3 global(clusx, clusy, clusz);
  TrkrDefs::subsurfkey subsurfkey = 0;
  Surface surface = tGeometry->get_tpc_surface_from_coords( hitsetkey, global, subsurfkey);

	// check if  surface not found and we can't track with it. 
  cout << " FIXME a0 " << endl;
  if(!surface) return {0, nullptr};
  cout << " FIXME a1 " << endl;

  global *= Acts::UnitConstants::cm;
  cout << " FIXME a2 " << endl;
  Acts::Vector3 local = surface->transform(tGeometry->geometry().getGeoContext()).inverse() * global;
  cout << " FIXME a3 " << endl;
  local /= Acts::UnitConstants::cm;     
  cout << " FIXME a4 " << endl;
  cluster->setLocalX (local(0));
  cout << " FIXME a5 " << endl;
  cluster->setLocalY (clust);
  cout << " FIXME a6 " << endl;

  // get the radius from the layerGeom, and make xyz global position
  // then follow along from TpcClusterizer to make local coordinates
  // + [ ] set the subsurf key, (get this from getting the surface from TpcCoordinates...)
  // the hitsetkey's have to match for reco and truth TrkrClusters to match
  //     (necessary but not sufficient

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

void TpcClusterBuilder::set_has_data() {
  has_data = ( 
      neff_electrons > 0. 
      && phi_integral != 0.
      && time_integral > 0.
      && hasPhiBins
      && hasTimeBins 
      && layerGeom != nullptr 
      && layer >= 0 
      && layer <  55
  );
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

ActsGeometry* TpcClusterBuilder::tGeometry = nullptr;
