#ifndef G4TPC_PHG4TPCPADPLANEREADOUT_H
#define G4TPC_PHG4TPCPADPLANEREADOUT_H

#include "PHG4TpcPadPlane.h"

#include <g4main/PHG4HitContainer.h>

#include <vector>

class TNtuple;
class TF1;
class PHG4CellContainer;
class PHG4CylinderCellGeomContainer;
class PHG4CylinderCellGeom;
class TrkrHitSetContainer;
class TrkrHitTruthContainer;

class PHG4TpcPadPlaneReadout : public PHG4TpcPadPlane
{
 public:
  PHG4TpcPadPlaneReadout(const std::string &name = "PHG4TpcPadPlaneReadout");
  virtual ~PHG4TpcPadPlaneReadout();

  int CreateReadoutGeometry(PHCompositeNode *topNode, PHG4CylinderCellGeomContainer *seggeo);

  void MapToPadPlane(PHG4CellContainer *g4cells, const double x_gem, const double y_gem, const double t_gem, PHG4HitContainer::ConstIterator hiter, TNtuple *ntpad, TNtuple *nthit);

  void MapToPadPlane(TrkrHitSetContainer *hitsetcontainer, TrkrHitTruthAssoc *hittruthassoc, const double x_gem, const double y_gem, const double t_gem, PHG4HitContainer::ConstIterator hiter, TNtuple *ntpad, TNtuple *nthit);

  void populate_rectangular_phibins(const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_zigzag_phibins(const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_zbins(const double z, const double cloud_sig_zz[2], std::vector<int> &adc_zbin, std::vector<double> &adc_zbin_share);

  void SetDefaultParameters();
  void UpdateInternalParameters();

 protected:
  std::string seggeonodename;

  TF1 *fcharge;
  TF1 *fpad[10];

  PHG4CylinderCellGeomContainer *GeomContainer;
  PHG4CylinderCellGeom *LayerGeom;

  double rad_gem;
  double output_radius;

  double neffelectrons_threshold;

  int MinLayer[3];
  int MaxLayer[3];
  double MinRadius[3];
  double MaxRadius[3];
  double Thickness[3];
  double MinZ;
  double MaxZ;
  double sigmaT;
  double sigmaL[2];
  double PhiBinWidth[3];
  double ZBinWidth;
  double tpc_drift_velocity;
  double tpc_adc_clock;

  int NZBins;
  int NPhiBins[3];

  int NTpcLayers[3];
  int tpc_region;
  int zigzag_pads;
  int hit;

  std::vector<int> adc_zbin;
  std::vector<int> pad_phibin;
  std::vector<double> pad_phibin_share;
  std::vector<double> adc_zbin_share;
};

#endif
