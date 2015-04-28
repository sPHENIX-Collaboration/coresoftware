// $$Id: RawTowerGeomv2.h,v 1.2 2014/10/29 16:55:26 pinkenbu Exp $$

/*!
 * \file RawTowerGeomv2.h
 * \brief generic tower description with variable eta binning
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2014/10/29 16:55:26 $$
 */

#ifndef RawTowerGeomv2_H__
#define RawTowerGeomv2_H__

#include "RawTowerGeom.h"

#include <map>
#include <string>

//! generic tower description with variable eta binning
class RawTowerGeomv2 : public RawTowerGeom
{
public:
  RawTowerGeomv2();
  RawTowerGeomv2(RawTowerGeom *geo);

  virtual
  ~RawTowerGeomv2()
  {
  }

  void
  identify(std::ostream& os = std::cout) const;
  double
  get_radius() const
  {
    return radius;
  }
  double
  get_thickness() const
  {
    return thickness;
  }
  int
  get_phibins() const
  {
    return nphibins;
  }
  double
  get_phistep() const
  {
    return phistep;
  }
  double
  get_phimin() const
  {
    return phimin;
  }
  int
  get_etabins() const
  {
    return eta_bound_map.size();
  }
  virtual double
  get_etastep() const
  {
    PHOOL_VIRTUAL_WARN(
        "get_etastep() is not supported in RawTowerGeomv2 (variable eta bins)");
    return NAN;
  }
  virtual double
  get_etamin() const
  {
    PHOOL_VIRTUAL_WARN(
        "get_etamin() is not supported in RawTowerGeomv2 (variable eta bins)");
    return NAN;
  }

  std::pair<double, double>
  get_phibounds(const int ibin) const;
  std::pair<double, double>
  get_etabounds(const int ibin) const;
  double
  get_etacenter(const int ibin) const;
  double
  get_phicenter(const int ibin) const;

  int
  get_etabin(const double eta) const;
  int
  get_phibin(const double phi) const;

  void
  set_radius(const double r)
  {
    radius = r;
  }
  void
  set_thickness(const double t)
  {
    thickness = t;
  }
  void
  set_phibins(const int i)
  {
    nphibins = i;
  }
  void
  set_phistep(const double phi)
  {
    phistep = phi;
  }
  void
  set_phimin(const double phi)
  {
    phimin = phi;
  }
  void
  set_etabins(const int i);

  virtual void
  set_etamin(const double z)
  {
    PHOOL_VIRTUAL_WARN(
        "set_etamin(const double) is not supported in RawTowerGeomv2 (variable eta bins)");
  }
  virtual void
  set_etastep(const double z)
  {
    PHOOL_VIRTUAL_WARN(
        "set_etastep(const double) is not supported in RawTowerGeomv2 (variable eta bins)");
  }
  void
  set_etabounds(const int ibin, const std::pair<double, double> & bounds);

protected:
  double radius;
  double thickness;
  int nphibins;
  double phimin;
  double phistep;

  typedef std::pair<double, double> bound_t;
  typedef std::vector<bound_t> bound_map_t;

  bound_map_t eta_bound_map;

ClassDef(RawTowerGeomv2,1)
};

#endif
