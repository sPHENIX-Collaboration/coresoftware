#ifndef QA_TPCCALIBRATIONQA_CONDITIONSUTILS_H
#define QA_TPCCALIBRATIONQA_CONDITIONSUTILS_H
// ROOT-independent numerical definitions used by the CDB snapshot exporter.
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>
namespace TpcConditionsDetail
{
  inline constexpr double referenceKeff = 1.55;
  inline constexpr double par0 = .185;
  inline constexpr double par1 = .815;
  inline constexpr double referenceSouth = 3.36138;
  inline constexpr double referenceNorth = 3.20357;

  inline double missing() { return std::numeric_limits<double>::quiet_NaN(); }
  inline const std::vector<std::string> &r1Channels(int side)
  {
    static const std::vector<std::string> south = {"S_01_R1_G4_IMon", "S_04_R1_G4_IMon",
                                                   "S_05_R1_G4_IMon", "S_10_R1_G4_IMon"};
    static const std::vector<std::string> north = {
        "N_03_R1_G4_IMon", "N_04_R1_G4_IMon", "N_06_R1_G4_IMon", "N_07_R1_G4_IMon",
        "N_08_R1_G4_IMon", "N_10_R1_G4_IMon", "N_11_R1_G4_IMon", "N_12_R1_G4_IMon"};
    return side == 0 ? south : north;
  }
  inline double median(std::vector<double> values)
  {
    if (values.empty() ||
        !std::all_of(values.begin(), values.end(), [](double x)
                     { return std::isfinite(x); }))
    {
      return missing();
    }
    std::sort(values.begin(), values.end());
    const auto n = values.size();
    return n % 2 ? values[n / 2] : .5 * (values[n / 2 - 1] + values[n / 2]);
  }
  struct Stats
  {
    std::size_t n{0};
    double mean{missing()}, rms{missing()}, low{missing()}, high{missing()};
  };
  inline Stats summarize(const std::vector<double> &values)
  {
    std::vector<double> valid;
    for (double x : values)
    {
      if (std::isfinite(x))
      {
        valid.push_back(x);
      }
    }
    Stats out;
    out.n = valid.size();
    if (valid.empty())
    {
      return out;
    }
    out.mean = std::accumulate(valid.begin(), valid.end(), 0.) / valid.size();
    double sum = 0.;
    for (double x : valid)
    {
      sum += (x - out.mean) * (x - out.mean);
    }
    out.rms = std::sqrt(sum / valid.size());
    out.low = *std::min_element(valid.begin(), valid.end());
    out.high = *std::max_element(valid.begin(), valid.end());
    return out;
  }
  inline double fallbackKeff(double current, int side, double scale = referenceKeff,
                             double offset = par0, double slope = par1, double south = referenceSouth,
                             double north = referenceNorth)
  {
    const double reference = side == 0 ? south : north;
    if (!std::isfinite(current) || current <= 0 || !std::isfinite(reference) || reference <= 0)
    {
      return missing();
    }
    return scale * (offset + slope * current / reference);
  }
}  // namespace TpcConditionsDetail
#endif
