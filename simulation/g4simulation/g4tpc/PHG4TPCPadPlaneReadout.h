#ifndef PHG4TPCPadPlaneReadout_h
#define PHG4TPCPadPlaneReadout_h

#include "PHG4TPCPadPlane.h"
#include <g4main/PHG4HitContainer.h>

#include <vector>

class TNtuple;
class TF1;
class PHG4CellContainer;
class PHG4CylinderCellGeomContainer;
class PHG4CylinderCellGeom;

class PHG4TPCPadPlaneReadout: public PHG4TPCPadPlane
{
public:
  PHG4TPCPadPlaneReadout(const std::string& name = "ReadoutPadPlane");
  virtual ~PHG4TPCPadPlaneReadout(){}

  int CreateReadoutGeometry(PHCompositeNode *topNode, PHG4CylinderCellGeomContainer *seggeo);

  void MapToPadPlane(PHG4CellContainer *g4cells, const double x_gem, const double y_gem, const double t_gem, PHG4HitContainer::ConstIterator hiter, TNtuple *ntpad, TNtuple *nthit);

  void populate_rectangular_phibins(const unsigned int layernum, const double phi,  const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_zigzag_phibins(const unsigned int layernum, const double phi,  const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_zbins( const double z,  const double cloud_sig_zz[2], std::vector<int> &adc_zbin, std::vector<double> &adc_zbin_share);

  void SetDefaultParameters();
  void UpdateInternalParameters();

protected:

  std::string seggeonodename;

  TF1 *fcharge;
  TF1 *fpad[10];

  //double max_active_radius;
  //double min_active_radius;
  //double rbinwidth;
  //double phibinwidth;
  //double tbinwidth;

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
