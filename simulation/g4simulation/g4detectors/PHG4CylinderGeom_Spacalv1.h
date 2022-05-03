// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $$Id: PHG4CylinderGeom_Spacalv1.h,v 1.2 2014/08/12 03:49:12 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2014/08/12 03:49:12 $$
 */
#ifndef G4DETECTORS_PHG4CYLINDERGEOMSPACALV1_H
#define G4DETECTORS_PHG4CYLINDERGEOMSPACALV1_H

#include "PHG4CylinderGeomv2.h"

#include <iostream>  // for cout, ostream
#include <map>
#include <string>

class PHParameters;

class PHG4CylinderGeom_Spacalv1 : public PHG4CylinderGeomv2
{
 public:
  /** @name Ctor DTor and IDs
   */
  ///@{
  PHG4CylinderGeom_Spacalv1();

  ~PHG4CylinderGeom_Spacalv1() override
  {
    sector_map.clear();
  }

  // from PHObject
  void identify(std::ostream &os = std::cout) const override;

  // from TObject
  void Print(Option_t *option = "") const override;

  virtual void SetDefault();

  //! load parameters from PHParameters, which interface to Database/XML/ROOT files
  void ImportParameters(const PHParameters &param) override;

  ///@}

  /** @name Set Cylinder Geometry
   */
  ///@{
  double
  get_max_radius() const
  {
    return get_radius() + get_thickness();
  }

  double
  get_half_radius() const
  {
    return get_radius() + get_thickness() / 2.;
  }

  double
  get_length() const
  {
    return get_zmax() - get_zmin();
  }

  double
  get_xpos() const
  {
    return xpos;
  }

  void
  set_xpos(double x)
  {
    this->xpos = x;
  }

  double
  get_ypos() const
  {
    return ypos;
  }

  void
  set_ypos(double y)
  {
    this->ypos = y;
  }

  double
  get_zpos() const
  {
    return zpos;
  }

  void
  set_zpos(double z)
  {
    this->zpos = z;
  }

  ///@}

  /** @name Azimuthal segments
   */
  ///@{
  virtual int
  get_azimuthal_n_sec() const;

  virtual double
  get_azimuthal_distance() const;

  virtual double
  get_z_distance() const;

  //! sector map sector_ID -> azimuthal rotation.
  typedef std::map<int, double> sector_map_t;

  //! sector map sector_ID -> azimuthal rotation.
  const sector_map_t &
  get_sector_map() const
  {
    return sector_map;
  }

  //! sector map sector_ID -> azimuthal rotation.
  sector_map_t &
  get_sector_map()
  {
    return sector_map;
  }

  //! load a default map that populate all the sectors
  void
  init_default_sector_map();

  ///@}

  /** @name Fiber geometry
   */
  ///@{
  double
  get_fiber_outer_r() const
  {
    return get_fiber_clading_thickness() + get_fiber_core_diameter() / 2;
  }

  double
  get_fiber_clading_thickness() const
  {
    return fiber_clading_thickness;
  }

  void
  set_fiber_clading_thickness(double fiberCladingThickness)
  {
    fiber_clading_thickness = fiberCladingThickness;
  }

  double
  get_fiber_core_diameter() const
  {
    return fiber_core_diameter;
  }

  void
  set_fiber_core_diameter(double fiberCoreDiameter)
  {
    fiber_core_diameter = fiberCoreDiameter;
  }

  double
  get_fiber_distance() const
  {
    return fiber_distance;
  }

  void
  set_fiber_distance(double fiberDistance)
  {
    fiber_distance = fiberDistance;
  }

  ///@}

  /** @name Materials
   */
  ///@{
  std::string
  get_absorber_mat() const
  {
    return absorber_mat;
  }

  void
  set_absorber_mat(const std::string &absorberMat)
  {
    absorber_mat = absorberMat;
  }

  std::string
  get_fiber_clading_mat() const
  {
    return fiber_clading_mat;
  }

  void
  set_fiber_clading_mat(const std::string &fiberCladingMat)
  {
    fiber_clading_mat = fiberCladingMat;
  }

  std::string
  get_fiber_core_mat() const
  {
    return fiber_core_mat;
  }

  void
  set_fiber_core_mat(const std::string &fiberCoreMat)
  {
    fiber_core_mat = fiberCoreMat;
  }

  ///@}

  /** @name General options
   */
  ///@{

  double
  get_fiber_core_step_size() const
  {
    return get_fiber_core_diameter() / 10.;
  }

  enum config_t
  {

    //! fiber always placed radially
    kNonProjective = 0,
    //! alias of above, more explicit
    k1DProjectiveSpacal = kNonProjective,

    //! Fully projective spacal with 2D tapered modules
    kFullProjective_2DTaper = 2,

    //! Fully projective spacal with 2D tapered modules. To speed up construction, same-length fiber is used cross one tower
    kFullProjective_2DTaper_SameLengthFiberPerTower = 3,

    //! Fully projective spacal with 2D tapered modules and allow azimuthal tilts
    kFullProjective_2DTaper_Tilted = 4,

    //! Fully projective spacal with 2D tapered modules and allow azimuthal tilts. To speed up construction, same-length fiber is used cross one tower
    kFullProjective_2DTaper_Tilted_SameLengthFiberPerTower = 5,

    //! alias of above, the default 2D-projective SPACAL
    k2DProjectiveSpacal = kFullProjective_2DTaper_Tilted_SameLengthFiberPerTower,

    //! max allowed value, for boundary cross check
    kInvalidSpacalConfig
  };

  config_t
  get_config() const
  {
    return config;
  }

  void
  set_config(config_t c)
  {
    this->config = c;
  }

  bool
  is_virualize_fiber() const
  {
    return virualize_fiber;
  }

  virtual bool is_azimuthal_seg_visible() const
  {
    return false;
  }

  void
  set_virualize_fiber(bool virualizeFiber = true)
  {
    virualize_fiber = virualizeFiber;
  }

  int get_construction_verbose() const
  {
    return construction_verbose;
  }

  void
  set_construction_verbose(int constructionVerbose)
  {
    construction_verbose = constructionVerbose;
  }

  ///@}

 protected:
  std::string absorber_mat;
  std::string fiber_core_mat;
  std::string fiber_clading_mat;
  double xpos;
  double ypos;
  double zpos;
  double fiber_core_diameter;
  double fiber_clading_thickness;
  double fiber_distance;
  config_t config;
  bool virualize_fiber;
  int construction_verbose;

  //! sector map sector_ID -> azimuthal rotation.
  sector_map_t sector_map;

  ClassDefOverride(PHG4CylinderGeom_Spacalv1, 2)
};

#endif
