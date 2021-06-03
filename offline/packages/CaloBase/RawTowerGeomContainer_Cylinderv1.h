#ifndef CALOBASE_RAWTOWERGEOMCONTAINER_CYLINDERV1_H
#define CALOBASE_RAWTOWERGEOMCONTAINER_CYLINDERV1_H

#include "RawTowerGeomContainerv1.h"

#include "RawTowerDefs.h"

#include <iostream>
#include <utility>
#include <vector>

/*! \class RawTowerGeomContainer_Cylinderv1
 \brief With additional description to conveniently use in central calorimeter with eta-phi bins
 */
class RawTowerGeomContainer_Cylinderv1 : public RawTowerGeomContainerv1
{
 public:
  RawTowerGeomContainer_Cylinderv1(
      RawTowerDefs::CalorimeterId caloid = RawTowerDefs::NONE);
  ~RawTowerGeomContainer_Cylinderv1() override { Reset(); }

  void
  identify(std::ostream& os = std::cout) const override;

  void Reset() override;

  double
  get_radius() const override
  {
    return radius;
  }

  double
  get_thickness() const override
  {
    return thickness;
  }

  int get_phibins() const override
  {
    return phi_bound_map.size();
  }
  int get_etabins() const override
  {
    return eta_bound_map.size();
  }

  std::pair<double, double>
  get_phibounds(const int ibin) const override;
  std::pair<double, double>
  get_etabounds(const int ibin) const override;

  double
  get_etacenter(const int ibin) const override;

  double
  get_phicenter(const int ibin) const override;

  int get_etabin(const double eta) const override;
  int get_phibin(const double phi) const override;

  void
  set_radius(const double r) override
  {
    radius = r;
  }
  void
  set_thickness(const double t) override
  {
    thickness = t;
  }
  void
  set_phibins(const int i) override;
  void
  set_etabins(const int i) override;

  void
  set_etabounds(const int ibin, const std::pair<double, double>& bounds) override;
  void
  set_phibounds(const int ibin, const std::pair<double, double>& bounds) override;

 protected:
  double radius;
  double thickness;

  typedef std::pair<double, double> bound_t;
  typedef std::vector<bound_t> bound_map_t;

  bound_map_t eta_bound_map;
  bound_map_t phi_bound_map;

  ClassDefOverride(RawTowerGeomContainer_Cylinderv1, 1)
};

#endif
