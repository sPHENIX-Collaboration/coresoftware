#ifndef RawTowerGeomContainer_Cylinderv1_H__
#define RawTowerGeomContainer_Cylinderv1_H__

#include "RawTowerGeomContainerv1.h"

/*! \class RawTowerGeomContainer_Cylinderv1
 \brief With additional description to conveniently use in central calorimeter with eta-phi bins
 */
class RawTowerGeomContainer_Cylinderv1 : public RawTowerGeomContainerv1
{

public:
  RawTowerGeomContainer_Cylinderv1(
      RawTowerDefs::CalorimeterId caloid  = RawTowerDefs::NONE);
  virtual ~RawTowerGeomContainer_Cylinderv1() {Reset();}

  virtual void
  identify(std::ostream& os = std::cout) const;

  virtual void Reset();

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
    return phi_bound_map.size();
  }
  int
  get_etabins() const
  {
    return eta_bound_map.size();
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
  set_phibins(const int i);
  void
  set_etabins(const int i);


  void
  set_etabounds(const int ibin, const std::pair<double, double> & bounds);
  void
  set_phibounds(const int ibin, const std::pair<double, double> & bounds);

protected:

  double radius;
  double thickness;

  typedef std::pair<double, double> bound_t;
  typedef std::vector<bound_t> bound_map_t;

  bound_map_t eta_bound_map;
  bound_map_t phi_bound_map;

ClassDef(RawTowerGeomContainer_Cylinderv1,1)
};

#endif /* RawTowerGeomContainer_Cylinderv1_H__ */
