#ifndef RHOUECALIB_TOWERSHAPEMAPMAKER_H
#define RHOUECALIB_TOWERSHAPEMAPMAKER_H

#include "RhoUECalibBase.h"

#include <array>
#include <string>

class PHCompositeNode;
class TH1D;
class TH3F;

/// \class TowerShapeMapMaker
///
/// \brief Step 2 of the rho UE calibration: sums for the tower shape map (s-map) of the
/// DetermineTowerRho seed jets
///
/// For every selected event and good tower of the retowered EMCal, IHCal and OHCal:
/// E_T = E / cosh(eta_corr) summed per (layer, ieta, iphi, class), with the number of
/// entries (good towers only), separately for even and odd selected events (the two
/// halves give the statistical precision of the map). Histograms sum_<L>_h<0|1> and
/// cnt_<L>_h<0|1>: TH3F (ieta 24, iphi 64, class = izbin * n_Q' + iQ'); h_nevt per class.
/// RhoUECalibPayload::MakeTowerShapePayload turns the merged sums into the payload of
/// DetermineTowerRho::set_tower_shape_path.
class TowerShapeMapMaker : public RhoUECalibBase
{
 public:
  explicit TowerShapeMapMaker(const std::string &outfile, const std::string &name = "TowerShapeMapMaker")
    : RhoUECalibBase(name, outfile)
  {
  }

 protected:
  int fill(PHCompositeNode *topNode) override;
  void book() override;

 private:
  std::array<std::array<TH3F *, 2>, NLAYER> m_sum{};
  std::array<std::array<TH3F *, 2>, NLAYER> m_cnt{};
  TH1D *m_nevt{nullptr};
};

#endif
