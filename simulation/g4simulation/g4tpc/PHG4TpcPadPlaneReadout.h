#ifndef G4TPC_PHG4TPCPADPLANEREADOUT_H
#define G4TPC_PHG4TPCPADPLANEREADOUT_H

#include "PHG4TpcPadPlane.h"
#include "TpcClusterBuilder.h"

#include <g4main/PHG4HitContainer.h>

#include <gsl/gsl_rng.h>

#include <array>
#include <climits>
#include <cmath>
#include <string>  // for string
#include <vector>

class PHCompositeNode;
class PHG4TpcCylinderGeomContainer;
class PHG4TpcCylinderGeom;
class TH2;
class TNtuple;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class PHG4TpcPadPlaneReadout : public PHG4TpcPadPlane
{
 public:
  PHG4TpcPadPlaneReadout(const std::string &name = "PHG4TpcPadPlaneReadout");

  ~PHG4TpcPadPlaneReadout() override;

  int InitRun(PHCompositeNode *topNode) override;

  void UseGain(const int flagToUseGain);
  void ReadGain();

  void SetDriftVelocity(double vd) override { drift_velocity = vd; }
  void SetReadoutTime(float t) override { extended_readout_time = t; }
  // otherwise warning of inconsistent overload since only one MapToPadPlane methow is overridden
  using PHG4TpcPadPlane::MapToPadPlane;

  void MapToPadPlane(TpcClusterBuilder &tpc_clustbuilder, TrkrHitSetContainer *single_hitsetcontainer, TrkrHitSetContainer *hitsetcontainer, TrkrHitTruthAssoc * /*hittruthassoc*/, const double x_gem, const double y_gem, const double t_gem, const unsigned int side, PHG4HitContainer::ConstIterator hiter, TNtuple * /*ntpad*/, TNtuple * /*nthit*/) override;

  void SetDefaultParameters() override;
  void UpdateInternalParameters() override;

 private:
  //  void populate_rectangular_phibins(const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_zigzag_phibins(const unsigned int side, const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_tbins(const double t, const std::array<double, 2> &cloud_sig_tt, std::vector<int> &adc_tbin, std::vector<double> &adc_tbin_share);

  double check_phi(const unsigned int side, const double phi, const double radius);

  PHG4TpcCylinderGeomContainer *GeomContainer = nullptr;
  PHG4TpcCylinderGeom *LayerGeom = nullptr;

  double neffelectrons_threshold = std::numeric_limits<double>::signaling_NaN();

  std::array<double, 3> MinRadius{};
  std::array<double, 3> MaxRadius{};

  static constexpr int NSides = 2;
  static constexpr int NSectors = 12;
  static const int NRSectors = 3;

  double sigmaT = std::numeric_limits<double>::signaling_NaN();
  std::array<double, 2> sigmaL{};
  std::array<double, 3> PhiBinWidth{};
  double drift_velocity = 8.0e-03;  // default value, override from macro
  float extended_readout_time = 0;  // ns
  int NTBins = std::numeric_limits<int>::max();
  int m_NHits = 0;
  // Using Gain maps is turned off by default
  int m_flagToUseGain = 0;
  // gaussian sampling
  static constexpr double _nsigmas = 5;

  double averageGEMGain = std::numeric_limits<double>::signaling_NaN();

  std::array<std::array<std::vector<double>, NRSectors>, NSides> sector_min_Phi_sectors;
  std::array<std::array<std::vector<double>, NRSectors>, NSides> sector_max_Phi_sectors;

  // return random distribution of number of electrons after amplification of GEM for each initial ionizing electron
  double getSingleEGEMAmplification();
  gsl_rng *RandomGenerator = nullptr;

  std::array<TH2 *, 2> h_gain{nullptr};
};

#endif
