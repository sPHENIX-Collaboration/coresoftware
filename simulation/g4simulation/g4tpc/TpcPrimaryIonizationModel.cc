#include "TpcPrimaryIonizationModel.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <stdexcept>
#include <vector>

namespace
{
  constexpr double argon_primary_clusters_per_cm = 23.0;
  constexpr double cf4_primary_clusters_per_cm = 51.0;
  constexpr double isobutane_primary_clusters_per_cm = 84.0;

  // Measured single-cluster electron multiplicities for Ar and CH4:
  // H. Fischle, J. Heintze, B. Schmidt, Nucl. Instrum. Meth. A 301
  // (1991) 202-214, doi:10.1016/0168-9002(91)90460-8.
  constexpr std::array<double, 19> argon_cluster_probabilities = {
      0.656, 0.150, 0.064, 0.035, 0.0225, 0.0155, 0.0105,
      0.0081, 0.0061, 0.0049, 0.0039, 0.0030, 0.0025, 0.0020,
      0.0016, 0.0012, 0.00095, 0.00075, 0.00063};

  constexpr std::array<double, 19> methane_cluster_probabilities = {
      0.786, 0.120, 0.034, 0.016, 0.0095, 0.0060, 0.0044,
      0.0034, 0.0027, 0.0021, 0.0017, 0.0013, 0.0010, 0.0008,
      0.0006, 0.00050, 0.00042, 0.00037, 0.00033};

  template <std::size_t Size>
  std::vector<double> build_cluster_probabilities(
      const std::array<double, Size> &measured)
  {
    std::vector<double> probabilities(
        TpcPrimaryIonizationModel::maximum_cluster_size + 1, 0.0);

    double measured_sum = 0.0;
    for (std::size_t index = 0; index < measured.size(); ++index)
    {
      probabilities[index + 1] = measured[index];
      measured_sum += measured[index];
    }

    if (!(measured_sum > 0.0 && measured_sum <= 1.0))
    {
      throw std::logic_error("invalid measured TPC cluster-size probabilities");
    }

    // The published table is used through n=19. Allocate its remaining
    // probability to a truncated 1/n^2 tail through the explicit cutoff.
    constexpr unsigned int tail_begin = Size + 1;
    double inverse_square_sum = 0.0;
    for (unsigned int size = tail_begin;
         size <= TpcPrimaryIonizationModel::maximum_cluster_size;
         ++size)
    {
      inverse_square_sum += 1.0 / static_cast<double>(size * size);
    }

    const double tail_normalization =
        (1.0 - measured_sum) / inverse_square_sum;
    for (unsigned int size = tail_begin;
         size <= TpcPrimaryIonizationModel::maximum_cluster_size;
         ++size)
    {
      probabilities[size] =
          tail_normalization / static_cast<double>(size * size);
    }

    return probabilities;
  }
}  // namespace

TpcPrimaryIonizationModel::TpcPrimaryIonizationModel(
    const GasFractions &fractions)
{
  const std::array<double, 3> values = {
      fractions.argon, fractions.cf4, fractions.isobutane};
  for (const double fraction : values)
  {
    if (!std::isfinite(fraction) || fraction < 0.0)
    {
      throw std::invalid_argument(
          "TPC primary-ionization gas fractions must be finite and nonnegative");
    }
  }

  const double fraction_sum =
      fractions.argon + fractions.cf4 + fractions.isobutane;
  if (!std::isfinite(fraction_sum) || std::abs(fraction_sum - 1.0) > 1.e-6)
  {
    throw std::invalid_argument(
        "TPC primary-ionization gas fractions must sum to one");
  }

  const double argon_weight =
      fractions.argon * argon_primary_clusters_per_cm;
  const double cf4_weight = fractions.cf4 * cf4_primary_clusters_per_cm;
  const double isobutane_weight =
      fractions.isobutane * isobutane_primary_clusters_per_cm;
  m_primaryClustersPerCm = argon_weight + cf4_weight + isobutane_weight;

  if (!(m_primaryClustersPerCm > 0.0))
  {
    throw std::invalid_argument(
        "TPC primary-ionization gas mixture has zero cluster density");
  }

  const auto argon_probabilities =
      build_cluster_probabilities(argon_cluster_probabilities);
  const auto methane_probabilities =
      build_cluster_probabilities(methane_cluster_probabilities);

  // No validated CF4 or isobutane table is supplied by this first model.
  // CF4 therefore uses the Ar shape and isobutane uses CH4 as an explicit
  // proxy. Each component is weighted by its contribution to the primary
  // cluster density rather than by volume fraction alone.
  const double argon_shape_weight =
      (argon_weight + cf4_weight) / m_primaryClustersPerCm;
  const double methane_shape_weight =
      isobutane_weight / m_primaryClustersPerCm;

  m_clusterSizeCdf.resize(maximum_cluster_size + 1, 0.0);
  double cumulative_probability = 0.0;
  for (unsigned int size = 1; size <= maximum_cluster_size; ++size)
  {
    cumulative_probability +=
        argon_shape_weight * argon_probabilities[size] +
        methane_shape_weight * methane_probabilities[size];
    m_clusterSizeCdf[size] = cumulative_probability;
  }
  m_clusterSizeCdf[maximum_cluster_size] = 1.0;
}

double TpcPrimaryIonizationModel::mean_primary_clusters(
    const double path_length_cm) const
{
  if (!std::isfinite(path_length_cm) || path_length_cm < 0.0)
  {
    throw std::invalid_argument(
        "TPC primary-ionization path length must be finite and nonnegative");
  }
  return m_primaryClustersPerCm * path_length_cm;
}

unsigned int TpcPrimaryIonizationModel::sample_cluster_size(
    const double uniform_01) const
{
  if (!std::isfinite(uniform_01) || uniform_01 < 0.0 || uniform_01 >= 1.0)
  {
    throw std::invalid_argument(
        "TPC cluster-size random variate must lie in [0, 1)");
  }

  const auto found = std::lower_bound(
      std::next(m_clusterSizeCdf.begin()),
      m_clusterSizeCdf.end(),
      uniform_01);
  return static_cast<unsigned int>(
      std::distance(m_clusterSizeCdf.begin(), found));
}

double TpcPrimaryIonizationModel::mean_cluster_size() const
{
  double previous_cdf = 0.0;
  double mean = 0.0;
  for (unsigned int size = 1; size <= maximum_cluster_size; ++size)
  {
    const double probability = m_clusterSizeCdf[size] - previous_cdf;
    mean += static_cast<double>(size) * probability;
    previous_cdf = m_clusterSizeCdf[size];
  }
  return mean;
}
