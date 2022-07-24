#include "TrkrTruthCentroidBuilder.h"
#include <cmath>

void TrkrTruthCentroidBuilder::add_electron(float phi, float z, float energy) {
  if (is_empty)  {
    is_empty = false;
    near_phi_boundary = (abs(phi) > 0.75*M_PI);
  }
  if (near_phi_boundary) phi += M_PI;
  if (phi >  M_PI) phi -= 2*M_PI;
  sum_phi   += phi*        energy;
  sumSq_phi += pow(phi,2.) * energy;
  sum_z     += z           * energy;
  sumSq_z   += pow(z,2.)   * energy;
  sumE      += energy;
}

std::pair<float,float> TrkrTruthCentroidBuilder::calc_mean_var(const double& Hj, const double& Ej, const double& Wj) {
    const float mean = Hj/Wj;
    const float var = std::sqrt(Ej/Wj-mean*mean);
    return { mean , var };
}

std::array<float,5> TrkrTruthCentroidBuilder::build_centroid() {
    auto phi_meanstd   = calc_mean_var(sum_phi,   sumSq_phi, sumE);
    if (near_phi_boundary) {
      phi_meanstd.first -= M_PI;
      if (phi_meanstd.first < -M_PI) phi_meanstd.first += 2*M_PI;
    }
    auto z_meanstd = calc_mean_var(sum_z, sumSq_z, sumE);
    return { phi_meanstd.first, phi_meanstd.second, z_meanstd.first, z_meanstd.second, static_cast<float>(sumE) };
}

void TrkrTruthCentroidBuilder::reset() {
  sum_phi = 0; sumSq_phi = 0; sum_z = 0; sumSq_z = 0; is_empty=true; near_phi_boundary=false;
}
