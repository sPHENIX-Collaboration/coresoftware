#ifndef RHOUECALIB_RHOUECALIBMAKER_H
#define RHOUECALIB_RHOUECALIBMAKER_H

#include "RhoUECalibBase.h"

#include <array>
#include <string>

class PHCompositeNode;
class TH2D;
class TH3D;

/// \class RhoUECalibMaker
///
/// \brief Step 3 of the rho UE calibration: sums for the eta-shape calibration of the
/// per-layer rho, read by TowerRhoToBackground
///
/// Runs after the three DetermineTowerRho (MULT; with the tower shape map if the
/// analysis uses it: the calibration must be made with the rho it is applied to).
/// Per selected event, layer and eta strip:
///   x = rho * cosh(eta_corr)   UE per tower energy predicted by rho
///   y = E_strip / N_good       mean energy of the good towers of the strip
/// summed with weight N_good per (layer, ieta, z bin, Q' bin): sw, sx, sy, sxx, sxy, syy
/// and the event count n (TH3D <sum>_<L>: ieta x z edges x Q' edges), and per class the
/// sums of z, Q' and events (TH2D zsum, qsum, nsum) for the interpolation nodes.
/// RhoUECalibPayload::MakeRhoUEPayload fits w = sy / sx and the affine y = alpha x + beta.
class RhoUECalibMaker : public RhoUECalibBase
{
 public:
  explicit RhoUECalibMaker(const std::string &outfile, const std::string &name = "RhoUECalibMaker")
    : RhoUECalibBase(name, outfile)
  {
  }

  void set_rho_nodes(const std::string &cemc, const std::string &hcalin, const std::string &hcalout)
  {
    m_rho_nodes = {cemc, hcalin, hcalout};
  }

  static constexpr int NSUM = 7;
  static const std::array<std::string, NSUM> &sum_names();  // sw sx sy sxx sxy syy n

 protected:
  int fill(PHCompositeNode *topNode) override;
  void book() override;

 private:
  std::array<std::string, NLAYER> m_rho_nodes{"TowerRho_CEMC_MULT", "TowerRho_HCALIN_MULT", "TowerRho_HCALOUT_MULT"};
  std::array<std::array<TH3D *, NSUM>, NLAYER> m_sums{};
  TH2D *m_zsum{nullptr};
  TH2D *m_qsum{nullptr};
  TH2D *m_nsum{nullptr};
};

#endif
