#ifndef G4TPC_PHG4TPCPADPLANEREADOUT_H
#define G4TPC_PHG4TPCPADPLANEREADOUT_H

#include "PHG4TpcPadPlane.h"

#include <g4main/PHG4HitContainer.h>

#include <gsl/gsl_rng.h>

#include <array>
#include <climits>
#include <cmath>
#include <string>  // for string
#include <vector>

class PHCompositeNode;
//class PHG4CellContainer;
class PHG4CylinderCellGeomContainer;
class PHG4CylinderCellGeom;
class TNtuple;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class PHG4TpcPadPlaneReadout : public PHG4TpcPadPlane
{
 public:
  PHG4TpcPadPlaneReadout(const std::string &name = "PHG4TpcPadPlaneReadout");

  ~PHG4TpcPadPlaneReadout() override;

  void SetDriftVelocity(double vd) override { drift_velocity = vd; }

  int CreateReadoutGeometry(PHCompositeNode * /*topNode*/, PHG4CylinderCellGeomContainer *seggeo) override;

  // otherwise warning of inconsistent overload since only one MapToPadPlane methow is overridden
  using PHG4TpcPadPlane::MapToPadPlane;

  void MapToPadPlane(TrkrHitSetContainer *single_hitsetcontainer, TrkrHitSetContainer *hitsetcontainer, TrkrHitTruthAssoc * /*hittruthassoc*/, const double x_gem, const double y_gem, const double t_gem, const unsigned int side, PHG4HitContainer::ConstIterator hiter, TNtuple * /*ntpad*/, TNtuple * /*nthit*/) override;

  void SetDefaultParameters() override;
  void UpdateInternalParameters() override;

 private:
  //  void populate_rectangular_phibins(const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_zigzag_phibins(const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_tbins(const double t, const std::array<double, 2> &cloud_sig_tt, std::vector<int> &adc_tbin, std::vector<double> &adc_tbin_share);

  std::string seggeonodename;

  PHG4CylinderCellGeomContainer *GeomContainer = nullptr;
  PHG4CylinderCellGeom *LayerGeom = nullptr;

  double rad_gem = NAN;
  double output_radius = 0;

  static const unsigned int print_layer = 18;

  double neffelectrons_threshold = NAN;

  std::array<int, 3> MinLayer;
  std::array<double, 3> MinRadius;
  std::array<double, 3> MaxRadius;
  std::array<double, 3> Thickness;
  double MaxZ = NAN;
  double MinT = NAN;
  double MaxT = NAN;
  double sigmaT = NAN;
  std::array<double, 2> sigmaL;
  std::array<double, 3> PhiBinWidth;
  double ZBinWidth = NAN;
  double TBinWidth = NAN;
  double drift_velocity = 8.0e-03;  // default value, override from macro
  double tpc_adc_clock = NAN;

  int NTBins = INT_MAX;
  std::array<int, 3> NPhiBins;
  std::array<int, 3> NTpcLayers;
  int m_NHits = 0;

  // gaussian sampling
  static constexpr double _nsigmas = 5;

  double averageGEMGain = NAN;

  std::vector<int> adc_tbin;
  std::vector<int> pad_phibin;
  std::vector<double> pad_phibin_share;
  std::vector<double> adc_tbin_share;

  // return random distribution of number of electrons after amplification of GEM for each initial ionizing electron
  double getSingleEGEMAmplification();
  gsl_rng *RandomGenerator;
};

#endif
