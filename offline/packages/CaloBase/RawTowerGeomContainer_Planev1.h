#ifndef RawTowerGeomContainer_Planev1_H__
#define RawTowerGeomContainer_Planev1_H__

#include "RawTowerGeomContainerv1.h"

/*! \class RawTowerGeomContainer_Planev1
 \brief With additional description to conveniently use in central calorimeter with eta-phi bins
 */
class RawTowerGeomContainer_Planev1 : public RawTowerGeomContainerv1
{
 public:
  RawTowerGeomContainer_Planev1(
      RawTowerDefs::CalorimeterId caloid = RawTowerDefs::NONE);
  virtual ~RawTowerGeomContainer_Planev1() { Reset(); }
  virtual void
  identify(std::ostream& os = std::cout) const;

  virtual void Reset();

  //! Topology of calorimeter
  virtual RawTowerDefs::CalorimeterTopology get_caloimeter_topology() const
  {
    return RawTowerDefs::PLANAR_CALO_TOPOLOGY;
  }

  //! center x of a planar calorimeter from front to back
  double get_center_x() const
  {
    return _center_x;
  }

  //! center x of a planar calorimeter from front to back
  void set_center_x(double centerX)
  {
    _center_x = centerX;
  }

  //! center y of a planar calorimeter from front to back
  double get_center_y() const
  {
    return _center_y;
  }

  //! center y of a planar calorimeter from front to back
  void set_center_y(double centerY)
  {
    _center_y = centerY;
  }

  //! center z of a planar calorimeter from front to back
  double get_center_z() const
  {
    return _center_z;
  }

  //! center z of a planar calorimeter from front to back
  void set_center_z(double centerZ)
  {
    _center_z = centerZ;
  }

  //! azimuthal angle of the orientation of a planar calorimeter from front to back from front to back
  double get_phi() const
  {
    return _phi;
  }

  //! azimuthal angle of the orientation of a planar calorimeter from front to back from front to back
  void set_phi(double phi)
  {
    _phi = phi;
  }

  //! polar angle of the orientation of a planar calorimeter from front to back from front to back
  double get_theta() const
  {
    return _theta;
  }

  //! polar angle of the orientation of a planar calorimeter from front to back from front to back
  void set_theta(double theta)
  {
    _theta = theta;
  }

 protected:
  //! center x of a planar calorimeter from front to back
  double _center_x;

  //! center y of a planar calorimeter from front to back
  double _center_y;

  //! center z of a planar calorimeter from front to back
  double _center_z;

  //! polar angle of the orientation of a planar calorimeter from front to back from front to back
  double _theta;

  //! azimuthal angle of the orientation of a planar calorimeter from front to back from front to back
  double _phi;

  ClassDef(RawTowerGeomContainer_Planev1, 1)
};

#endif /* RawTowerGeomContainer_Planev1_H__ */
