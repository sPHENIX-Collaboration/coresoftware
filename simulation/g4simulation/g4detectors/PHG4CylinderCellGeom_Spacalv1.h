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

#include "PHG4CylinderCellGeom.h"
#include <map>

/*!
 * \brief PHG4CylinderCellGeom_Spacalv1
 */
class PHG4CylinderCellGeom_Spacalv1 : public PHG4CylinderCellGeom
{
public:
  PHG4CylinderCellGeom_Spacalv1();
  virtual
  ~PHG4CylinderCellGeom_Spacalv1();
  void
  identify(std::ostream& os = std::cout) const;

  virtual std::pair<double, double>
  get_zbounds(const int ibin) const;
  virtual std::pair<double, double>
  get_etabounds(const int ibin) const;

  virtual double
  get_etacenter(const int ibin) const;
  virtual double
  get_zcenter(const int ibin) const;

  virtual int
  get_etabin(const double eta) const;
  virtual int
  get_zbin(const double z) const;

  void
  set_zbounds(const int ibin, const std::pair<double, double> & bounds);
  void
  set_etabounds(const int ibin, const std::pair<double, double> & bounds);

  typedef std::pair<double, double> bound_t;
  typedef std::map<int, bound_t> bound_map_t;

  const bound_map_t &
  get_eta_bound_map() const
  {
    map_consistency_check();
    return eta_bound_map;
  }

  void
  set_eta_bound_map(const bound_map_t & etaBoundMap)
  {
    eta_bound_map = etaBoundMap;
  }

  const bound_map_t &
  get_z_bound_map() const
  {
    map_consistency_check();
    return z_bound_map;
  }

  void
  set_z_bound_map(const bound_map_t & boundMap)
  {
    z_bound_map = boundMap;
  }

  //! map tower_z_ID -> eta_bin number
  typedef std::map<int, int> tower_z_ID_eta_bin_map_t;

  //! map tower_z_ID -> eta_bin number
  const tower_z_ID_eta_bin_map_t &
  get_tower_z_ID_eta_bin_map() const
  {
    return tower_z_ID_eta_bin_map;
  }

  virtual int
  get_etabin_block(const int tower_z_ID) const;

  //! map tower_z_ID -> eta_bin number for blocks
  void
  set_tower_z_ID_eta_bin_map(const tower_z_ID_eta_bin_map_t & m)
  {
    tower_z_ID_eta_bin_map = m;
  }

protected:

  void
  map_consistency_check() const;

  bound_map_t z_bound_map;
  bound_map_t eta_bound_map;

  //! map tower_z_ID -> eta_bin number for blocks
  tower_z_ID_eta_bin_map_t tower_z_ID_eta_bin_map;

ClassDef(PHG4CylinderCellGeom_Spacalv1,2)

};

#endif /* PHG4CYLINDERCELLGEOMSPACALV1_H_ */
