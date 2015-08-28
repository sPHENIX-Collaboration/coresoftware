// $Id: $                                                                                             

/*!
 * \file PHG4CylinderCellGeomSpacalv1.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHG4CYLINDERCELLGEOMSPACALV1_H_
#define PHG4CYLINDERCELLGEOMSPACALV1_H_

#include <PHG4CylinderCellGeom.h>
#include <vector>

/*!
 * \brief PHG4CylinderCellGeom_Spacalv1
 */
class PHG4CylinderCellGeom_Spacalv1 : public PHG4CylinderCellGeom
{
public:
  PHG4CylinderCellGeom_Spacalv1();
  virtual
  ~PHG4CylinderCellGeom_Spacalv1();

  virtual std::pair<double, double> get_zbounds(const int ibin) const;
  virtual std::pair<double, double> get_etabounds(const int ibin) const;

  virtual double get_etacenter(const int ibin) const;
  virtual double get_zcenter(const int ibin) const;

  virtual int get_etabin(const double eta) const;
  virtual int get_zbin(const double z) const;

  void
  set_zbounds(const int ibin, const std::pair<double, double> & bounds);
  void
  set_etabounds(const int ibin, const std::pair<double, double> & bounds);

  typedef std::pair<double, double> bound_t;
  typedef std::vector<bound_t> bound_map_t;

  const bound_map_t & get_eta_bound_map() const
    {
      return eta_bound_map;
    }

  void
  set_eta_bound_map(const bound_map_t & etaBoundMap)
  {
    eta_bound_map = etaBoundMap;
  }

  const bound_map_t &
  get_bound_map() const
  {
    return z_bound_map;
  }

  void
  set_bound_map(const bound_map_t &  boundMap)
  {
    z_bound_map = boundMap;
  }

protected:

  void map_consistency_check();

  bound_map_t z_bound_map;
  bound_map_t eta_bound_map;

  ClassDef(PHG4CylinderCellGeom_Spacalv1,1)

};

#endif /* PHG4CYLINDERCELLGEOMSPACALV1_H_ */
