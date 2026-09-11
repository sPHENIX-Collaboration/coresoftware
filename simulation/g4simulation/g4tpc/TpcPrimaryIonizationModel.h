// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4TPC_TPCPRIMARYIONIZATIONMODEL_H
#define G4TPC_TPCPRIMARYIONIZATIONMODEL_H

#include <vector>

/**
 * Fast, path-length-based, MIP-calibrated TPC primary-ionization model for
 * non-electron charged particles. PHG4TpcElectronDrift retains its established
 * energy-deposit response for electrons and positrons.
 *
 * The model supplies the mean primary-cluster density for a gas mixture and
 * samples the number of electrons in one primary cluster. Random-number
 * generation remains the responsibility of the caller so that this class is
 * deterministic and independently testable.
 */
class TpcPrimaryIonizationModel
{
 public:
  struct GasFractions
  {
    double argon{0.0};
    double cf4{0.0};
    double isobutane{0.0};
  };

  static constexpr unsigned int maximum_cluster_size = 384;

  explicit TpcPrimaryIonizationModel(const GasFractions &fractions);

  [[nodiscard]] double primary_clusters_per_cm() const
  {
    return m_primaryClustersPerCm;
  }

  [[nodiscard]] double mean_primary_clusters(const double path_length_cm) const;
  [[nodiscard]] unsigned int sample_cluster_size(const double uniform_01) const;
  [[nodiscard]] double mean_cluster_size() const;

 private:
  double m_primaryClustersPerCm{0.0};
  std::vector<double> m_clusterSizeCdf;
};

#endif  // G4TPC_TPCPRIMARYIONIZATIONMODEL_H
