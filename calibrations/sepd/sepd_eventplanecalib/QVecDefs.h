#ifndef QVECDEFS_H
#define QVECDEFS_H

#include <array>
#include <string>
#include <format>
#include <vector>

namespace QVecShared
{
  static constexpr size_t CENT_BINS = 8;
  static constexpr std::array<int, 3> HARMONICS = {2, 3, 4};
  static constexpr int SEPD_CHANNELS = 744;

  enum class ChannelStatus : int
  {
    Good = 0,
    Dead = 1,
    Hot = 2,
    Cold = 3
  };

  enum class Subdetector : size_t
  {
    S = 0,
    N = 1,
    NS = 2,
    Count = 3
  };

  enum class QComponent
  {
    X,
    Y
  };

  struct QVec
  {
    double x{0.0};
    double y{0.0};
  };

 /**
   * @brief Centralized helper to generate standard histogram names for the sEPD calibration.
   * * Standardizes the naming convention: h_sEPD_Q_{det}_{var}_{n}{suffix}_avg
   * * @param det The detector arm ("S" for South, "N" for North).
   * @param var The physics variable or moment (e.g., "x", "y", "xx", "xy").
   * @param n The harmonic order (e.g., 2, 3, 4).
   * @param suffix Optional pass-specific suffix (e.g., "_corr", "_corr2").
   * @return A formatted std::string representing the ROOT histogram name.
   */
  inline std::string get_hist_name(const std::string& det, const std::string& var, int n, const std::string& suffix = "")
  {
    return std::format("h_sEPD_Q_{}_{}_{}{}_avg", det, var, n, suffix);
  }
}  // namespace QVecShared

#endif // QVECDEFS_H
