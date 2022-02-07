// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $$Id: PHG4CylinderGeom_Spacalv3.h,v 1.3 2014/08/28 22:18:35 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2014/08/28 22:18:35 $$
 */
#ifndef G4DETECTORS_PHG4CYLINDERGEOMSPACALV3_H
#define G4DETECTORS_PHG4CYLINDERGEOMSPACALV3_H

#include "PHG4CylinderGeom_Spacalv2.h"

#include <iostream>  // for operator<<, basic_ostream::op...
#include <map>
#include <string>
#include <utility>  // std::pair, std::make_pair

class PHParameters;

class PHG4CylinderGeom_Spacalv3 : public PHG4CylinderGeom_Spacalv2
{
 public:
  PHG4CylinderGeom_Spacalv3();

  ~PHG4CylinderGeom_Spacalv3() override;

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;

  // from TObject
  void Print(Option_t* option = "") const override;

  void SetDefault() override;

  //! load parameters from PHParameters, which interface to Database/XML/ROOT files
  void ImportParameters(const PHParameters& param) override;

  double
  get_sidewall_outer_torr() const
  {
    return sidewall_outer_torr;
  }

  void
  set_sidewall_outer_torr(double sidewallOuterTorr)
  {
    sidewall_outer_torr = sidewallOuterTorr;
  }

  double
  get_sidewall_thickness() const
  {
    return sidewall_thickness;
  }

  void
  set_sidewall_thickness(double sidewallThickness)
  {
    sidewall_thickness = sidewallThickness;
  }

  std::string
  get_sidewall_mat() const
  {
    return sidewall_mat;
  }

  void
  set_sidewall_mat(const std::string& absorberMat)
  {
    sidewall_mat = absorberMat;
  }

  int get_max_phi_bin_in_sec() const
  {
    return max_phi_bin_in_sec;
  }

  void
  set_max_phi_bin_in_sec(int maxPhiBinInSec)
  {
    max_phi_bin_in_sec = maxPhiBinInSec;
  }

  const std::string& get_divider_mat() const
  {
    return divider_mat;
  }

  void set_divider_mat(const std::string& dividerMat)
  {
    divider_mat = dividerMat;
  }

  double get_divider_width() const
  {
    return divider_width;
  }

  void set_divider_width(double dividerWidth)
  {
    divider_width = dividerWidth;
  }

  class geom_tower
  {
   public:
    geom_tower();
    virtual ~geom_tower()
    {
    }

    int id;
    double pDz;

    double pDy1;
    double pDx1;
    double pDx2;
    double pDy2;
    double pDx3;
    double pDx4;

    double pTheta;
    double pPhi;
    double pAlp1;
    double pAlp2;

    double pRotationAngleX;
    double centralX;
    double centralY;
    double centralZ;

    double ModuleSkinThickness;

    //! number of fiber along final azimuthal direction
    int NFiberX;
    //! number of fiber along final polar direction
    int NFiberY;

    //! number of fiber along final azimuthal direction
    int NSubtowerX;
    //! number of fiber along final polar direction
    int NSubtowerY;

    //! fiber layout (2D index of 0...NFiberX/NFiberY) -> fiber_id
    int
    compose_fiber_id(int index_x, int index_y) const;
    //! fiber_id -> sub tower ID x/azimuthal direction: 0 ... NSubtowerX -1
    int
    get_sub_tower_ID_x(int fiber_id) const;
    //! fiber_id -> sub tower ID y/polar direction: 0 ... NSubtowerY -1
    int
    get_sub_tower_ID_y(int fiber_id) const;
    //! fiber_id -> fraction position in sub tower ID in the x/azimuthal direction, [0-1]
    double
    get_position_fraction_x_in_sub_tower(int fiber_id) const;
    //! fiber_id -> fraction position in sub tower ID in the y/polar direction, [0-1]
    double
    get_position_fraction_y_in_sub_tower(int fiber_id) const;

    //! height of light guide
    double LightguideHeight;
    //! edge length ratio, narrow / wide
    double LightguideTaperRatio;
    //! edge length ratio, narrow / wide
    std::string LightguideMaterial;

    virtual void
    identify(std::ostream& os = std::cout) const;

    //! read via PHParameters
    void
    ImportParameters(const PHParameters& param,
                     const std::string& param_prefix);

    ClassDef(PHG4CylinderGeom_Spacalv3::geom_tower, 3)
  };
  typedef std::map<int, geom_tower> tower_map_t;

  void
  load_demo_sector_tower_map1();
  void
  load_demo_sector_tower_map2();
  void
  load_demo_sector_tower_map4();

  const tower_map_t&
  get_sector_tower_map() const
  {
    return sector_tower_map;
  }

  //! get approximate radial position of tower
  double
  get_tower_radial_position(const geom_tower& tower) const;

  //! check that all towers has consistent sub-tower divider
  void
  subtower_consistency_check() const;
  //! sub-tower divider along the polar direction
  int get_n_subtower_eta() const;
  //! sub-tower divider along the azimuthal direction
  int get_n_subtower_phi() const;
  //! max tolerance needed for the light guide
  double
  get_max_lightguide_height() const;

  //
  //  void
  //  set_geom_super_tower_map(geom_super_tower_map_t geomSuperTowerMap)
  //  {
  //    geom_super_tower_map = geomSuperTowerMap;
  //  }

  //! compact ID of each fiber in 32bit PHG4Hit::set_scint_id(). Buffer the result for repeated use.
  class scint_id_coder
  {
   public:
    explicit scint_id_coder(int scint_id);
    scint_id_coder(int sector_id, int tower_id, int fiber_id);
    virtual ~scint_id_coder()
    {
    }

    virtual void
    identify(std::ostream& os = std::cout) const
    {
      os << "scint_id_coder with "
         << "scint_ID(" << scint_ID << ") = "
         << "sector_ID(" << sector_ID << "), "
         << "tower_ID(" << tower_ID
         << "), "
         << "fiber_ID(" << fiber_ID << ")" << std::endl;
    }

    int scint_ID;
    int sector_ID;
    int tower_ID;
    int fiber_ID;

    static const int kfiber_bit = 13;  // max 8192 fiber per tower
    static const int ktower_bit = 11;  // max 2048 towers per sector
    static const int ksector_bit = 8;  // max 256 sectors

    ClassDef(PHG4CylinderGeom_Spacalv3::scint_id_coder, 1)
  };

  //! convert tower_ID + sector ID to eta and z bins as in other cylindrical calorimeters
  //! @return: a std::pair of zbin and phibin number
  virtual std::pair<int, int>
  get_tower_z_phi_ID(const int tower_ID, const int sector_ID) const;

 protected:
  double sidewall_thickness;
  double sidewall_outer_torr;
  std::string sidewall_mat;
  int max_phi_bin_in_sec;

  tower_map_t sector_tower_map;

  //! wdith along the approximate radial direction
  double divider_width;
  //! material for divider
  std::string divider_mat;

  ClassDefOverride(PHG4CylinderGeom_Spacalv3, 4)
};

#endif
