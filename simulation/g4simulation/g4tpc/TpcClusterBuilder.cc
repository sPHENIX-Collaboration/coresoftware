#include "TpcClusterBuilder.h"
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TpcDefs.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <algorithm>
#include <Acts/Surfaces/Surface.hpp>

class TpcClusterBuilder;
using std::cout, std::endl;

TpcClusterBuilder& TpcClusterBuilder::operator+=(const TpcClusterBuilder& rhs) {
  // only builders with the same hitsetkey (and therefore side, layerGeom, etc...) are added
  int rhs_phi_bin_lo = rhs.phi_bin_lo;
  int rhs_phi_bin_hi = rhs.phi_bin_hi;

  // each individual set of phi_bin_{lo,hi} are already checked to not wrap-around;
  // it remains only to check if the the rhs vs *this side are wrapped relative to each other;
  // if they are, then wrap the smaller one around to the larger one
  if ( (phi_bin_lo - rhs_phi_bin_lo) > (nphibins/2)) {
    rhs_phi_bin_lo += nphibins;
    rhs_phi_bin_hi += nphibins;
  } else if ( (rhs_phi_bin_lo - phi_bin_lo) > (nphibins/2)) {
    phi_bin_lo += nphibins;
    phi_bin_hi += nphibins;
  }

  if (rhs_phi_bin_lo < phi_bin_lo) phi_bin_lo = rhs_phi_bin_lo;
  if (rhs_phi_bin_hi > phi_bin_hi) phi_bin_hi = rhs_phi_bin_hi;

  // extend the timing bins
  if (rhs.time_bin_lo < time_bin_lo) time_bin_lo = rhs.time_bin_lo;
  if (rhs.time_bin_hi > time_bin_hi) time_bin_hi = rhs.time_bin_hi;

  // add electrons
  neff_electrons += rhs.neff_electrons;
  phi_integral   += rhs.phi_integral;
  time_integral  += rhs.time_integral;

  return *this;
}

/* TpcClusterBuilder::PairCluskeyCluster TpcClusterBuilder::build(MapHitsetkeyUInt& cluster_cnt) const { */
TrkrCluster* TpcClusterBuilder::build() const {
  int phi_size = phi_bin_hi-phi_bin_lo+1;
  if (phi_size > CHAR_MAX || phi_size < 0) return nullptr;

  int Z_size = time_bin_hi-time_bin_lo+1;
  if (Z_size > CHAR_MAX || Z_size < 0)  return nullptr;

  TrkrClusterv4* cluster = new TrkrClusterv4();
  cluster->setPhiSize ( static_cast<char>(phi_size) );
  cluster->setZSize   ( static_cast<char>(Z_size)   );
  cluster->setAdc     ( neff_electrons );

  // copy logic from packages/tpc/TpcClusterizer.cc ~line 333
  const double clusphi { phi_integral  / neff_electrons };
  const double clust   { time_integral / neff_electrons };
  const double radius = layerGeom->get_radius();  // returns center of layer
  const float  clusx  = radius * cos(clusphi);
  const float  clusy  = radius * sin(clusphi);

  const unsigned short NTBins = (unsigned short) layerGeom->get_zbins();
  const double m_tdriftmax = AdcClockPeriod * NTBins / 2.0;  
  const double zdriftlength = clust * tGeometry->get_drift_velocity();
  const double z_sign = (side==0) ? -1 : 1;
  const double clusz  = z_sign * (m_tdriftmax * tGeometry->get_drift_velocity()-zdriftlength);

  Acts::Vector3 global(clusx, clusy, clusz);
  TrkrDefs::subsurfkey subsurfkey = 0;
  Surface surface = tGeometry->get_tpc_surface_from_coords( hitsetkey, global, subsurfkey);

	// check if  surface not found and we can't track with it. 
  if(!surface) {
    delete cluster;
    return nullptr;
  };

  global *= Acts::UnitConstants::cm;
  Acts::Vector3 local = surface->transform(tGeometry->geometry().getGeoContext()).inverse() * global;
  local /= Acts::UnitConstants::cm;     
  cluster->setLocalX (local(0));
  cluster->setLocalY (clust);
  return cluster;
}

void TpcClusterBuilder::set_has_data() {
  has_data = ( 
      neff_electrons   > 0. 
      && phi_integral  > 0.
      && time_integral > 0.
      && hasPhiBins
      && hasTimeBins 
      && layerGeom != nullptr 
      && layer >= 0 
      && layer <  55
      && nphibins > 0
  );
  if (!has_data) return;

  // Make the hitsetkey
  // copy logic from packages/tpc/TpcClusterizer.cc ~line 333
  const double clusphi { phi_integral  / neff_electrons };
  const int    phi_pad_number        = layerGeom->get_phibin(clusphi);
  const auto   phibins               = layerGeom->get_phibins();
  const unsigned int pads_per_sector = phibins / 12;
  const unsigned int sector          = phi_pad_number / pads_per_sector;
  hitsetkey = TpcDefs::genHitSetKey(layer, sector, side);
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
