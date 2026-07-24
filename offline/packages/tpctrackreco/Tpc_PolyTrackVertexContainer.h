#pragma once

#include <phool/PHObject.h>

#include <iostream>

class Tpc_PolyTrackVertex;

class Tpc_PolyTrackVertexContainer : public PHObject
{
 public:
  Tpc_PolyTrackVertexContainer() = default;
  ~Tpc_PolyTrackVertexContainer() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Tpc_PolyTrackVertexContainer base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  virtual unsigned int size() const { return 0; }
  virtual void add_vertex(Tpc_PolyTrackVertex*) {}
  virtual const Tpc_PolyTrackVertex* get_vertex(unsigned int) const { return nullptr; }
  virtual Tpc_PolyTrackVertex* get_vertex(unsigned int) { return nullptr; }

  virtual int get_collision_vertex_valid() const { return 0; }
  virtual unsigned int get_collision_vertex_count() const { return 0; }
  virtual double get_collision_x(unsigned int = 0) const { return 0.0; }
  virtual double get_collision_y(unsigned int = 0) const { return 0.0; }
  virtual double get_collision_z(unsigned int = 0) const { return 0.0; }
  virtual double get_collision_z_rms(unsigned int = 0) const { return 0.0; }
  virtual unsigned int get_collision_ntracks(unsigned int = 0) const { return 0; }
  virtual unsigned int get_collision_min_clusters() const { return 0; }

  virtual void set_collision_vertex_valid(int) {}
  virtual void set_collision_min_clusters(unsigned int) {}
  virtual void clear_collision_vertices() {}
  virtual void add_collision_vertex(double, double, double, double, unsigned int) {}

 private:
  ClassDefOverride(Tpc_PolyTrackVertexContainer, 0)
};
