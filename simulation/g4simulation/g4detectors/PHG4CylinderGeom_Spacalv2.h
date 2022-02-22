// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $$Id: PHG4CylinderGeom_Spacalv2.h,v 1.3 2014/08/28 22:18:35 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2014/08/28 22:18:35 $$
 */
#ifndef G4DETECTORS_PHG4CYLINDERGEOMSPACALV2_H
#define G4DETECTORS_PHG4CYLINDERGEOMSPACALV2_H

#include "PHG4CylinderGeom_Spacalv1.h"

#include <cmath>
#include <iostream>  // for cout, ostream

class PHParameters;

class PHG4CylinderGeom_Spacalv2 : public PHG4CylinderGeom_Spacalv1
{
 public:
  PHG4CylinderGeom_Spacalv2();

  ~PHG4CylinderGeom_Spacalv2() override
  {
  }

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;

  // from TObject
  void Print(Option_t* option = "") const override;

  void SetDefault() override;

  //! load parameters from PHParameters, which interface to Database/XML/ROOT files
  void ImportParameters(const PHParameters& param) override;

  int get_azimuthal_n_sec() const override;

  virtual void set_azimuthal_n_sec(int azimuthalNSec);

  //! azimuthal tilt in rad
  double
  get_azimuthal_tilt() const
  {
    return azimuthal_tilt;
  }

  //! azimuthal tilt in rad
  void
  set_azimuthal_tilt(double azimuthalTilt)
  {
    azimuthal_tilt = azimuthalTilt;
  }

  bool is_azimuthal_seg_visible() const override;

  virtual void set_azimuthal_seg_visible(bool b = true);

  double get_polar_taper_ratio() const
  {
    return polar_taper_ratio;
  }

  void
  set_polar_taper_ratio(double polarTaperRatio)
  {
    polar_taper_ratio = polarTaperRatio;
  }

  double get_half_polar_taper_angle() const;

  double
  get_sec_azimuthal_width() const;

  double
  get_sec_depth() const
  {
    const double available_depth = get_thickness() - (sqrt((get_radius()) * (get_radius()) + (get_sec_azimuthal_width() / 2) * (get_sec_azimuthal_width() / 2)) - get_radius()) - get_assembly_spacing();
    if (available_depth < get_sec_azimuthal_width())
      return NAN;
    else
      return sqrt(
          available_depth * available_depth - get_sec_azimuthal_width() * get_sec_azimuthal_width());
  }

  double
  get_block_width() const
  {
    return get_sec_azimuthal_width() - 0.0001;
  }

  double
  get_block_depth() const
  {
    return sqrt(get_sec_depth() * get_sec_depth() - get_polar_taper_ratio() * get_polar_taper_ratio() * get_block_width() * get_block_width()) - get_assembly_spacing();
  }

  double get_assembly_spacing() const
  {
    return assembly_spacing;
  }

  void
  set_assembly_spacing(double assemblySpacing)
  {
    assembly_spacing = assemblySpacing;
  }

  //! regulated fiber distance in the tapering direction
  double
  get_reg_fiber_grid_distance_taper() const;

  //! regulated fiber distance in the non-tapering direction
  double
  get_reg_fiber_grid_distance_nontaper() const;

 protected:
  int azimuthal_n_sec;

  //! azimuthal tilt in rad
  double azimuthal_tilt;
  bool azimuthal_seg_visible;
  double polar_taper_ratio;
  double assembly_spacing;

  ClassDefOverride(PHG4CylinderGeom_Spacalv2, 2)
};

#endif
