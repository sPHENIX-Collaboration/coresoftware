#pragma once

#include <phool/PHObject.h>

#include <iostream>

class Tpc_PolyTrackVertex : public PHObject
{
 public:
  Tpc_PolyTrackVertex() = default;
  ~Tpc_PolyTrackVertex() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Tpc_PolyTrackVertex base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  virtual unsigned int get_track_id() const { return 0; }
  virtual unsigned int get_source_assembled_track_id() const { return 0; }
  virtual double get_dca2d() const { return 0.0; }
  virtual double get_z0() const { return 0.0; }
  virtual int get_pca_valid() const { return 0; }
  virtual double get_pca_x() const { return 0.0; }
  virtual double get_pca_y() const { return 0.0; }
  virtual double get_pca_z() const { return 0.0; }
  virtual double get_pca_radius() const { return 0.0; }
  virtual double get_pca_phi() const { return 0.0; }

  virtual void set_track_id(unsigned int) {}
  virtual void set_source_assembled_track_id(unsigned int) {}
  virtual void set_dca2d(double) {}
  virtual void set_z0(double) {}
  virtual void set_pca_valid(int) {}
  virtual void set_pca_x(double) {}
  virtual void set_pca_y(double) {}
  virtual void set_pca_z(double) {}
  virtual void set_pca_radius(double) {}
  virtual void set_pca_phi(double) {}

 private:
  ClassDefOverride(Tpc_PolyTrackVertex, 0)
};
