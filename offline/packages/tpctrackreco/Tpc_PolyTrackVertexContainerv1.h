#pragma once

#include "Tpc_PolyTrackVertexContainer.h"

#include <iostream>
#include <vector>

class Tpc_PolyTrackVertex;

class Tpc_PolyTrackVertexContainerv1 : public Tpc_PolyTrackVertexContainer
{
 public:
  Tpc_PolyTrackVertexContainerv1();
  ~Tpc_PolyTrackVertexContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override { return static_cast<unsigned int>(m_vertices.size()); }
  void add_vertex(Tpc_PolyTrackVertex* vtx) override { m_vertices.push_back(vtx); }
  const Tpc_PolyTrackVertex* get_vertex(unsigned int i) const override;
  Tpc_PolyTrackVertex* get_vertex(unsigned int i) override;

  int get_collision_vertex_valid() const override { return m_collision_vertex_valid; }
  unsigned int get_collision_vertex_count() const override { return static_cast<unsigned int>(m_collision_x.size()); }
  double get_collision_x(unsigned int i = 0) const override;
  double get_collision_y(unsigned int i = 0) const override;
  double get_collision_z(unsigned int i = 0) const override;
  double get_collision_z_rms(unsigned int i = 0) const override;
  unsigned int get_collision_ntracks(unsigned int i = 0) const override;
  unsigned int get_collision_min_clusters() const override { return m_collision_min_clusters; }

  void set_collision_vertex_valid(int v) override { m_collision_vertex_valid = v; }
  void set_collision_min_clusters(unsigned int v) override { m_collision_min_clusters = v; }
  void clear_collision_vertices() override;
  void add_collision_vertex(double x, double y, double z, double z_rms, unsigned int ntracks) override;

 private:
  std::vector<Tpc_PolyTrackVertex*> m_vertices;
  int m_collision_vertex_valid {0};
  std::vector<double> m_collision_x;
  std::vector<double> m_collision_y;
  std::vector<double> m_collision_z;
  std::vector<double> m_collision_z_rms;
  std::vector<unsigned int> m_collision_ntracks;
  unsigned int m_collision_min_clusters {0};

  ClassDefOverride(Tpc_PolyTrackVertexContainerv1, 1)
};
