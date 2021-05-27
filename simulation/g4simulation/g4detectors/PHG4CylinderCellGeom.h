// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELLGEOM_H
#define G4DETECTORS_PHG4CYLINDERCELLGEOM_H

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <string>
#include <utility>           // for pair

class PHG4CylinderCellGeom: public PHObject
{
 public:
  PHG4CylinderCellGeom();

  ~PHG4CylinderCellGeom() override {}

// from PHObject
  void identify(std::ostream& os = std::cout) const override;

  int get_layer() const {return layer;}
  double get_radius() const {return radius;}
  double get_thickness() const {return thickness;}
  int get_binning() const {return binning;}
  int get_zbins() const;
  int get_phibins() const;
  double get_zmin() const;
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

  virtual int get_etabin(const double eta) const;
  virtual int get_zbin(const double z) const;
  virtual int get_phibin(const double phi) const;

   void set_layer(const int i) {layer = i;}
   void set_binning(const int i) {binning = i;}
   void set_radius(const double r) {radius = r;}
   void set_thickness(const double t) {thickness = t;}
   void set_zbins(const int i);
   void set_zmin(const double z);
   void set_zstep(const double z);
   void set_phibins(const int i);
   void set_phistep(const double phi);
   void set_phimin(const double phi);
   void set_etabins(const int i);
   void set_etamin(const double z);
   void set_etastep(const double z);
  
 protected:
  void check_binning_method(const int i) const;
  void check_binning_method_eta(const std::string & src = "") const;
  void check_binning_method_phi(const std::string & src = "") const;
  std::string methodname(const int i) const;
  int layer;
  int binning;
  double radius;
  int nzbins;
  double zmin;
  double zstep;
  int nphibins;
  double phimin;
  double phistep;
  double thickness;

  ClassDefOverride(PHG4CylinderCellGeom,1)
};

#endif
