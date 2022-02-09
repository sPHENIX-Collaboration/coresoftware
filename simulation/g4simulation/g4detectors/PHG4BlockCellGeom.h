// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKCELLGEOM_H
#define G4DETECTORS_PHG4BLOCKCELLGEOM_H

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream
#include <string>
#include <utility>  // for pair

class PHG4BlockCellGeom : public PHObject
{
 public:
  PHG4BlockCellGeom();

  ~PHG4BlockCellGeom() override {}

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;

  int get_layer() const { return _layer; }
  int get_binning() const { return _binning; }

  int get_zbins() const;
  double get_zstep() const;
  double get_zmin() const;

  int get_xbins() const;
  double get_xstep() const;
  double get_xmin() const;

  int get_etabins() const;
  double get_etastep() const;
  double get_etamin() const;

  std::pair<double, double> get_zbounds(const int ibin) const;
  std::pair<double, double> get_etabounds(const int ibin) const;
  std::pair<double, double> get_xbounds(const int ibin) const;

  double get_etacenter(const int ibin) const;
  double get_zcenter(const int ibin) const;
  double get_xcenter(const int ibin) const;

  int get_etabin(const double eta) const;
  int get_zbin(const double z) const;
  int get_xbin(const double x) const;

  void set_layer(const int i) { _layer = i; }
  void set_binning(const int i) { _binning = i; }

  void set_zbins(const int i);
  void set_zmin(const double z);
  void set_zstep(const double z);

  void set_xbins(const int i);
  void set_xstep(const double x);
  void set_xmin(const double x);

  void set_etabins(const int i);
  void set_etamin(const double z);
  void set_etastep(const double z);

 protected:
  void check_binning_method(const int i) const;
  void check_binning_method_eta(const std::string& src = "") const;
  void check_binning_method_x(const std::string& src = "") const;
  std::string methodname(const int i) const;

  int _layer;
  int _binning;

  int _nzbins;
  double _zmin;
  double _zstep;

  int _nxbins;
  double _xmin;
  double _xstep;

  ClassDefOverride(PHG4BlockCellGeom, 1)
};

#endif
