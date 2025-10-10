// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4TPCGEOMV1_H
#define G4DETECTORS_PHG4TPCGEOMV1_H

#include "PHG4TpcGeom.h"

#include <array>
#include <cmath>
#include <iostream>  // for cout, ostream
#include <string>
#include <utility>  // for pair

class PHG4TpcGeomv1 : public PHG4TpcGeom
{
 public:
  PHG4TpcGeomv1() = default;

  ~PHG4TpcGeomv1() override = default;

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;

  int get_layer() const override { return layer; }
  double get_radius() const override { return radius; }
  double get_thickness() const override { return thickness; }
  int get_binning() const override { return binning; }
  int get_zbins() const override;
  int get_phibins() const override;
  double get_zmin() const override;
  double get_phistep() const override;
  double get_phimin() const override;
  double get_zstep() const override;
  int get_etabins() const override;
  double get_etastep() const override;
  double get_etamin() const override;

  double get_max_driftlength() const override { return max_driftlength; }
  double get_CM_halfwidth() const override  { return CM_halfwidth; }
  double get_adc_clock() const override  { return  adc_clock; } // default sim value
  double get_extended_readout_time() const override  { return extended_readout_time; }
  double get_drift_velocity_sim() const override  { return drift_velocity_sim; }
  
  std::pair<double, double> get_zbounds(const int ibin) const override;
  std::pair<double, double> get_phibounds(const int ibin) const override;
  std::pair<double, double> get_etabounds(const int ibin) const override;
  double get_etacenter(const int ibin) const override;
  double get_zcenter(const int ibin) const override;
  double get_phicenter(const int ibin, const int side = 0) const override;
  double get_phicenter_new(const int ibin) const override;
  double get_phi(const float ibin, const int side = 0) const override;
  
  int get_etabin(const double eta) const override;
  int get_zbin(const double z) const override;
  int get_phibin(const double phi, int side = 0) const override;
  int get_phibin_new(const double phi) const override;
  
  float get_pad_float(const double phi, int side = 0) const override;
  float get_tbin_float(const double z) const override;
  int find_phibin(const double phi, int side = 0) const override;

  void set_layer(const int i) override { layer = i; }
  void set_binning(const int i) override { binning = i; }
  void set_radius(const double r) override { radius = r; }
  void set_thickness(const double t) override { thickness = t; }
  void set_zbins(const int i) override;
  void set_zmin(const double z) override;
  void set_zstep(const double z) override;
  void set_phibins(const int i) override;
  void set_phistep(const double phi) override;
  void set_phimin(const double phi) override;
  void set_etabins(const int i) override;
  void set_etamin(const double z) override;
  void set_etastep(const double z) override;
  // capture the z geometry related setup parameters
  void set_max_driftlength(const double val) override { max_driftlength = val; }
  void set_CM_halfwidth(const double val) override { CM_halfwidth = val; }
  void set_adc_clock(const double val) override { adc_clock = val; }
  void set_extended_readout_time(const double val) override { extended_readout_time = val; }
  void set_drift_velocity_sim(const double val) override { drift_velocity_sim = val; }
  
  static const int NSides = 2;

  void set_r_bias(const std::array<std::vector<double>, NSides> &dr) override { sector_R_bias = dr; }
  void set_phi_bias(const std::array<std::vector<double>, NSides> &dphi) override { sector_Phi_bias = dphi; }

  void set_sector_min_phi(const std::array<std::vector<double>, NSides> &s_min_phi) override
  {
    sector_min_Phi = s_min_phi;
  }
  void set_sector_max_phi(const std::array<std::vector<double>, NSides> &s_max_phi) override
  {
    sector_max_Phi = s_max_phi;
  }

 const std::array<std::vector<double>, NSides> &get_sector_min_phi() override
  {
    return sector_min_Phi;
  }
  const std::array<std::vector<double>, NSides> &get_sector_max_phi() override
  {
    return sector_max_Phi;
  }

 protected:
  void check_binning_method(const int i) const;
  void check_binning_method_eta(const std::string& src = "") const;
  void check_binning_method_phi(const std::string& src = "") const;
  std::string methodname(const int i) const;
  int layer{-999};
  int binning{0};
  double radius{std::numeric_limits<double>::quiet_NaN()};
  int nzbins{-1};
  double zmin{std::numeric_limits<double>::quiet_NaN()};
  double zstep{std::numeric_limits<double>::quiet_NaN()};
  int nphibins{-1};
  double phimin{-M_PI};
  double phistep{std::numeric_limits<double>::quiet_NaN()};
  double thickness{std::numeric_limits<double>::quiet_NaN()};

  double max_driftlength{std::numeric_limits<double>::quiet_NaN()};
  double CM_halfwidth{std::numeric_limits<double>::quiet_NaN()};
  double adc_clock{std::numeric_limits<double>::quiet_NaN()};
  double extended_readout_time{std::numeric_limits<double>::quiet_NaN()};
  double drift_velocity_sim{std::numeric_limits<double>::quiet_NaN()};
  
  std::array<std::vector<double>, NSides> sector_R_bias;
  std::array<std::vector<double>, NSides> sector_Phi_bias;
  std::array<std::vector<double>, NSides> sector_min_Phi;
  std::array<std::vector<double>, NSides> sector_max_Phi;

  // streamer
  friend std::ostream& operator<<(std::ostream&, const PHG4TpcGeomv1&);

  ClassDefOverride(PHG4TpcGeomv1, 1)
};

#endif
