// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4TPCGEOM_H
#define G4DETECTORS_PHG4TPCGEOM_H

#include <phool/PHObject.h>

#include <phool/phool.h>

#include <cmath>
#include <iostream>  // for cout, ostream
#include <array>

class PHParameters;

class PHG4TpcGeom : public PHObject
{
 public:
  ~PHG4TpcGeom() override {}

  // from PHObject
  void identify(std::ostream &os = std::cout) const override;
  
  virtual int get_layer() const
  {
    PHOOL_VIRTUAL_WARN("get_layer()");
    return -99999;
  }
  virtual double get_radius() const
  {
    PHOOL_VIRTUAL_WARN("get_radius()");
    return NAN;
  }
  virtual double get_thickness() const
  {
    PHOOL_VIRTUAL_WARN("get_thickness()");
    return NAN;
  }
  virtual double get_zmin() const
  {
    PHOOL_VIRTUAL_WARN("get_zmin()");
    return NAN;
  }
  virtual double get_zmax() const
  {
    PHOOL_VIRTUAL_WARN("get_zmax()");
    return NAN;
  }
  virtual int get_binning() const
  {
    PHOOL_VIRTUAL_WARN("get_binning()");
    return -99999;
  }
  virtual int get_zbins() const
  {
    PHOOL_VIRTUAL_WARN("get_zbins()");
    return -99999;
  }
  virtual int get_phibins() const
  {
    PHOOL_VIRTUAL_WARN("get_phibins()");
    return -99999;
  }
  virtual double get_phistep() const
  {
    PHOOL_VIRTUAL_WARN("get_phistep()");
    return NAN;
  }
  virtual double get_phimin() const
  {
    PHOOL_VIRTUAL_WARN("get_phimin()");
    return NAN;
  }
  virtual double get_zstep() const
  {
    PHOOL_VIRTUAL_WARN("get_zstep()");
    return NAN;
  }
  virtual int get_etabins() const
  {
    PHOOL_VIRTUAL_WARN("get_etabins()");
    return -99999;
  }
  virtual double get_etastep() const
  {
    PHOOL_VIRTUAL_WARN("get_etastep()");
    return NAN;
  }
  virtual double get_etamin() const
  {
    PHOOL_VIRTUAL_WARN("get_etamin()");
    return NAN;
  }
  virtual double get_max_driftlength() const
  {
    PHOOL_VIRTUAL_WARN("get_max_driftlength()");
    return NAN;
  }
  virtual double get_CM_halfwidth() const
  {
    PHOOL_VIRTUAL_WARN("get_CM_halfwidth()");
    return NAN;
  }
  virtual double get_adc_clock() const
  {
    PHOOL_VIRTUAL_WARN("get_adc_clock()");
    return NAN;
  }
  virtual double get_extended_readout_time() const
  {
    PHOOL_VIRTUAL_WARN("get_extended_readout_time()");
    return NAN;
  }
  virtual double get_drift_velocity_sim() const
  {
    PHOOL_VIRTUAL_WARN("get_drift_velocity_sim()");
    return NAN;
  }
  virtual std::pair<double, double> get_zbounds(const int) const
  {
    PHOOL_VIRTUAL_WARN("get_zbounds()");
    return std::make_pair(NAN,NAN);
  }
  virtual std::pair<double, double> get_phibounds(const int) const
   {
    PHOOL_VIRTUAL_WARN("get_zbounds()");
    return std::make_pair(NAN,NAN);
  }
  virtual std::pair<double, double> get_etabounds(const int) const
   {
    PHOOL_VIRTUAL_WARN("get_zbounds()");
    return std::make_pair(NAN,NAN);
  }
  virtual double get_etacenter(const int) const
  {
    PHOOL_VIRTUAL_WARN("get_etacenter()");
    return NAN;
  }
  virtual double get_zcenter(const int) const
  {
    PHOOL_VIRTUAL_WARN("get_zcenter()");
    return NAN;
  }    
  virtual double get_phicenter(const int, const int=0) const
  {
    PHOOL_VIRTUAL_WARN("get_phicenter()");
    return NAN;
  }
  virtual double get_phicenter_new(const int) const
  {
    PHOOL_VIRTUAL_WARN("get_phicenter_new()");
    return NAN;
  }
  virtual double get_phi(const float, const int) const
  {
    PHOOL_VIRTUAL_WARN("get_phi()");
    return NAN;
  }
  virtual int get_etabin(const double) const
  {
    PHOOL_VIRTUAL_WARN("get_etabin()");
    return -99999;
  }
  virtual int get_zbin(const double) const
  {
    PHOOL_VIRTUAL_WARN("get_zbin()");
    return -99999;
  }
  virtual int get_phibin(const double, int=0) const
  {
    PHOOL_VIRTUAL_WARN("get_phibin()");
    return -99999;
  }
  virtual int get_phibin_new(const double) const
  {
    PHOOL_VIRTUAL_WARN("get_phibin_new()");
    return -99999;
  }
  virtual float get_pad_float(const double, int) const
  {
    PHOOL_VIRTUAL_WARN("get_pad_float()");
    return -99999;
  }
  virtual float get_tbin_float(const double) const
  {
    PHOOL_VIRTUAL_WARN("get_tbin_float()");
    return -99999;
  }
  virtual int find_phibin(const double, int) const
  {
    PHOOL_VIRTUAL_WARN("get_phibin()");
    return -99999;
  }

  virtual const std::array<std::vector<double>, 2> &get_sector_min_phi();
  virtual const std::array<std::vector<double>, 2> &get_sector_max_phi();

  virtual void set_sector_min_phi(const std::array<std::vector<double>, 2>&)
  {
    PHOOL_VIRTUAL_WARN("set_sector_min_phi(const std::array<std::vector<double>, 2>&)");
  }
  virtual void set_sector_max_phi(const std::array<std::vector<double>, 2>&)
  {
    PHOOL_VIRTUAL_WARN("set_sector_max_phi(const std::array<std::vector<double>, 2>&)");
  }
  virtual void set_r_bias(const std::array<std::vector<double>, 2> &)
  {
    PHOOL_VIRTUAL_WARN("set_r_bias(const std::array<std::vector<double>, 2>&)");  
  }
  virtual void set_phi_bias(const std::array<std::vector<double>, 2> &)
  {
    PHOOL_VIRTUAL_WARN("set_phi_bias(const std::array<std::vector<double>, 2>&)");  
  }
  
  virtual void set_layer(const int) { PHOOL_VIRTUAL_WARN("set_layer(const int)"); }
  virtual void set_radius(const double) { PHOOL_VIRTUAL_WARN("set_radius(const double)"); }
  virtual void set_thickness(const double) { PHOOL_VIRTUAL_WARN("set_thickness(const double)"); }
  virtual void set_zmin(const double) { PHOOL_VIRTUAL_WARN("set_zmin(const double)"); }
  virtual void set_zmax(const double) { PHOOL_VIRTUAL_WARN("set_zmax(const double)"); }
  
  virtual void set_binning(const int) { PHOOL_VIRTUAL_WARN("set_binning(const int)"); }
  virtual void set_zbins(const int) { PHOOL_VIRTUAL_WARN("set_zbins(const int)"); }
  virtual void set_zstep(const double) { PHOOL_VIRTUAL_WARN("set_zstep(const double)"); }  
  virtual void set_phibins(const int) { PHOOL_VIRTUAL_WARN("set_phibins(const int)"); }

  virtual void set_phistep(const double) { PHOOL_VIRTUAL_WARN("set_phistep(const double)"); }  
  virtual void set_phimin(const double) { PHOOL_VIRTUAL_WARN("set_phimin(const double)"); }  


  virtual void set_etabins(const int) { PHOOL_VIRTUAL_WARN("set_etabins(const int)"); }
  virtual void set_etamin(const double) { PHOOL_VIRTUAL_WARN("set_etamin(const double)"); }
  virtual void set_etamax(const double) { PHOOL_VIRTUAL_WARN("set_etamax(const double)"); }
  virtual void set_etastep(const double) { PHOOL_VIRTUAL_WARN("set_etastep(const double)"); }

  virtual void set_max_driftlength(const double) { PHOOL_VIRTUAL_WARN("set_max_driftlength(const double)"); }
  virtual void set_CM_halfwidth(const double) { PHOOL_VIRTUAL_WARN("set_CM_halfwidth(const double)"); }
  virtual void set_adc_clock(const double) { PHOOL_VIRTUAL_WARN("set_adc_clock(const double)"); }
  virtual void set_extended_readout_time(const double) { PHOOL_VIRTUAL_WARN("set_extended_readout_time(const double)"); }
  virtual void set_drift_velocity_sim(const double) { PHOOL_VIRTUAL_WARN("set_drift_velocity_sim(const double)"); }
  
  //! load parameters from PHParameters, which interface to Database/XML/ROOT files
  virtual void ImportParameters(const PHParameters & /*param*/) { return; }

 protected:
  PHG4TpcGeom() {}

  ClassDefOverride(PHG4TpcGeom, 1)
};

#endif
