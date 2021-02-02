#ifndef G4TPC_PHG4TPCPADPLANEREADOUT_H
#define G4TPC_PHG4TPCPADPLANEREADOUT_H

#include "PHG4TpcPadPlane.h"

#include <g4main/PHG4HitContainer.h>

#include <gsl/gsl_rng.h>

#include <array>
#include <climits>
#include <cmath>
#include <string>                     // for string
#include <vector>

class PHCompositeNode;
class PHG4CellContainer;
class PHG4CylinderCellGeomContainer;
class PHG4CylinderCellGeom;
class TF1;
class TNtuple;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class PHG4TpcPadPlaneReadout : public PHG4TpcPadPlane
{
 public:
  PHG4TpcPadPlaneReadout(const std::string &name = "PHG4TpcPadPlaneReadout");

  virtual ~PHG4TpcPadPlaneReadout();

  int CreateReadoutGeometry(PHCompositeNode *topNode, PHG4CylinderCellGeomContainer *seggeo);

  void MapToPadPlane(PHG4CellContainer *g4cells, const double x_gem, const double y_gem, const double t_gem, PHG4HitContainer::ConstIterator hiter, TNtuple *ntpad, TNtuple *nthit);

  void MapToPadPlane(TrkrHitSetContainer *single_hitsetcontainer, TrkrHitSetContainer *hitsetcontainer, TrkrHitTruthAssoc *hittruthassoc, const double x_gem, const double y_gem, const double t_gem, PHG4HitContainer::ConstIterator hiter, TNtuple *ntpad, TNtuple *nthit);

  void SetDefaultParameters();
  void UpdateInternalParameters();

  private:

  void populate_rectangular_phibins(const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_zigzag_phibins(const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_zbins(const double z, const std::array<double,2>& cloud_sig_zz, std::vector<int> &adc_zbin, std::vector<double> &adc_zbin_share);

  std::string seggeonodename;

  PHG4CylinderCellGeomContainer *GeomContainer = nullptr;
  PHG4CylinderCellGeom *LayerGeom = nullptr;

  double rad_gem = NAN;
  double output_radius = 0;

  static const unsigned int print_layer = 18;

  double neffelectrons_threshold = NAN;

  std::array<int,3> MinLayer;
  std::array<double,3> MinRadius;
  std::array<double,3> MaxRadius;
  std::array<double,3> Thickness;
  double MinZ = NAN;
  double MaxZ = NAN;
  double sigmaT = NAN;
  std::array<double,2> sigmaL;
  std::array<double,3> PhiBinWidth;
  double ZBinWidth = NAN;
  double tpc_drift_velocity = NAN;
  double tpc_adc_clock = NAN;

  int NZBins = INT_MAX;
  std::array<int,3> NPhiBins;
  std::array<int,3> NTpcLayers;
  int zigzag_pads = INT_MAX;
  int hit = 0;

  // gaussian sampling
  static constexpr double _nsigmas = 5;

  double averageGEMGain = NAN;

  std::vector<int> adc_zbin;
  std::vector<int> pad_phibin;
  std::vector<double> pad_phibin_share;
  std::vector<double> adc_zbin_share;

  // return random distribution of number of electrons after amplification of GEM for each initial ionizing electron
  double getSingleEGEMAmplification();
  gsl_rng *RandomGenerator;

};

#endif
