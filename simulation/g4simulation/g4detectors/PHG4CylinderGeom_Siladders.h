#ifndef PHG4CylinderGeom_Siladders_H__
#define PHG4CylinderGeom_Siladders_H__

#include "PHG4CylinderGeom.h"

#include <phool/phool.h>
#include <cmath>

class PHG4CylinderGeom_Siladders: public PHG4CylinderGeom
  {
  public:
    PHG4CylinderGeom_Siladders();
    PHG4CylinderGeom_Siladders(
      const int    layer_,
      const double strip_x_,
      const double strip_y_,
      const double strip_z0_,
      const double strip_z1_,
      const int    nstrips_z_sensor0_,
      const int    nstrips_z_sensor1_,
      const int    nstrips_phi_cell_,
      const int    nladders_layer_,
      const double ladder_z0_,
      const double ladder_z1_,
      const double sensor_radius_,
      const double strip_x_offset_,
      const double offsetphi_,
      const double offsetrot_) :
        layer(layer_),
        strip_x(strip_x_),
        strip_y(strip_y_),
        strip_z0(strip_z0_),
        strip_z1(strip_z1_),
        nstrips_z_sensor0(nstrips_z_sensor0_),
        nstrips_z_sensor1(nstrips_z_sensor1_),
        nstrips_phi_cell(nstrips_phi_cell_),
        nladders_layer(nladders_layer_),
        ladder_z0(ladder_z0_),
        ladder_z1(ladder_z1_),
        sensor_radius(sensor_radius_),
        strip_x_offset(strip_x_offset_),
        offsetphi(offsetphi_),
	offsetrot(offsetrot_),
	radius(NAN)
    {
      // Type-A
      strip_z_[0]          = strip_z0;
      ladder_z_[0]         = ladder_z0;
      nstrips_z_sensor_[0] = nstrips_z_sensor0;

      // Type-B
      strip_z_[1]          = strip_z1;
      ladder_z_[1]         = ladder_z1;
      nstrips_z_sensor_[1] = nstrips_z_sensor1;

      dphi_ = 2.*M_PI/(double)nladders_layer;
    }

    virtual ~PHG4CylinderGeom_Siladders();

    void identify(std::ostream& os = std::cout) const;
    void set_layer(const int i)
    {
      layer = i;
    }

    int get_layer() const
      {
        return layer;
      }

    double get_radius() const
      {
        //return sensor_radius_inner;
	return sensor_radius;
      }

    bool load_geometry();
    void find_segment_center(const int segment_z_bin, const int segment_phi_bin, double location[]);
    void find_strip_center(  const int segment_z_bin, const int segment_phi_bin, const int strip_column, const int strip_index, double location[]);
    void find_strip_index_values(const int segment_z_bin, const double ypos, const double zpos,  int &strip_y_index, int &strip_z_index);
    void find_strip_center_localcoords(const int segment_z_bin, const int strip_y_index, const int strip_z_index, double location[]);

    double get_thickness() const
      {
        return strip_x;
      }

    double get_strip_y_spacing() const
      {
        return strip_y;
      }

    double get_strip_z_spacing() const
      {
        return strip_z0;
      }

    double get_strip_tilt() const
      {
        return 0.;
      }

    double get_strip_phi_tilt() const
      {
        return offsetrot;
      }

  protected:

    int layer;
    double strip_x;
    double strip_y;
    double strip_z0;
    double strip_z1;
    int nstrips_z_sensor0;
    int nstrips_z_sensor1;
    int nstrips_phi_cell;
    int nladders_layer;
    double ladder_z0;
    double ladder_z1;
    double sensor_radius;
    double strip_x_offset;
    double offsetphi;
    double offsetrot;

    double strip_z_[2];
    double ladder_z_[2];
    int nstrips_z_sensor_[2];
    double dphi_;

    double radius;

    ClassDef(PHG4CylinderGeom_Siladders,1)
  };

#endif
