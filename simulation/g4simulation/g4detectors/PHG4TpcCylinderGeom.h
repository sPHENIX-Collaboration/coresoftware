// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4TPCCYLINDERGEOM_H
#define G4DETECTORS_PHG4TPCCYLINDERGEOM_H

//#include <phool/PHObject.h>
#include "PHG4CylinderGeom.h"

#include <cmath>
#include <iostream>  // for cout, ostream
#include <string>
#include <utility>  // for pair
#include <array>

class PHG4TpcCylinderGeom : public PHG4CylinderGeom
{
 public:
  PHG4TpcCylinderGeom(){};

  ~PHG4TpcCylinderGeom(){};

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;

  int get_layer() const override { return layer; }
  double get_radius() const override { return radius; }
  double get_thickness() const override { return thickness; }
  int get_binning() const { return binning; }
  int get_zbins() const;
  int get_phibins() const;
  double get_zmin() const override;
  double get_phistep() const;
  double get_phimin() const;
  double get_zstep() const;
  int get_etabins() const;
  double get_etastep() const;
  double get_etamin() const;

  virtual std::pair<double, double> get_zbounds(const int ibin) const;
  virtual std::pair<double, double> get_phibounds(const int ibin) const;
  virtual std::pair<double, double> get_etabounds(const int ibin) const;
  virtual double get_etacenter(const int ibin) const;
  virtual double get_zcenter(const int ibin) const;
  virtual double get_phicenter(const int ibin) const;
  virtual double get_phicenter_new(const int ibin) const;
  virtual double get_phi(const float ibin) const;

  virtual int get_etabin(const double eta) const;
  virtual int get_zbin(const double z) const;
  virtual int get_phibin(const double phi, int side = 0) const;
  virtual int get_phibin_new(const double phi) const;
  
  virtual int find_phibin(const double phi, int side = 0) const;

  void set_layer(const int i) override { layer = i; }
  void set_binning(const int i) { binning = i; }
  void set_radius(const double r) override { radius = r; }
  void set_thickness(const double t) override { thickness = t; }
  void set_zbins(const int i);
  void set_zmin(const double z) override;
  void set_zstep(const double z);
  void set_phibins(const int i);
  void set_phistep(const double phi);
  void set_phimin(const double phi);
  void set_etabins(const int i);
  void set_etamin(const double z);
  void set_etastep(const double z);

  static const int NSides = 2;

  void set_r_bias(std::array<std::vector<double>, NSides > dr){sector_R_bias = dr;}
  void set_phi_bias(std::array<std::vector<double>, NSides > dphi){sector_Phi_bias = dphi;}

  //void set_sector_phi(const double s_phi){sector_phi = s_phi;}
  void set_sector_min_phi(std::array<std::vector<double>, NSides >  s_min_phi){sector_min_Phi = s_min_phi;}
  void set_sector_max_phi(std::array<std::vector<double>, NSides >  s_max_phi){sector_max_Phi = s_max_phi;}

  std::array<std::vector<double>, NSides > get_sector_min_phi(){return sector_min_Phi;}
  std::array<std::vector<double>, NSides > get_sector_max_phi(){return sector_max_Phi;}



 protected:
  void check_binning_method(const int i) const;
  void check_binning_method_eta(const std::string& src = "") const;
  void check_binning_method_phi(const std::string& src = "") const;
  std::string methodname(const int i) const;
  int layer = -999;
  int binning = 0;
  double radius = NAN;
  int nzbins = -1;
  double zmin = NAN;
  double zstep = NAN;
  int nphibins = -1;
  double phimin = -M_PI;
  double phistep = NAN;
  double thickness = NAN;

  std::array<std::vector<double>, NSides > sector_R_bias;
  std::array<std::vector<double>, NSides > sector_Phi_bias;
  std::array<std::vector<double>, NSides > sector_min_Phi;
  std::array<std::vector<double>, NSides > sector_max_Phi;


  ClassDefOverride(PHG4TpcCylinderGeom, 2)
};

#endif
