#include "Tpc_PolyTrackVertexContainerv1.h"
#include "Tpc_PolyTrackVertex.h"

ClassImp(Tpc_PolyTrackVertexContainerv1)

Tpc_PolyTrackVertexContainerv1::Tpc_PolyTrackVertexContainerv1()
{
  Reset();
}

Tpc_PolyTrackVertexContainerv1::~Tpc_PolyTrackVertexContainerv1()
{
  Reset();
}

void Tpc_PolyTrackVertexContainerv1::identify(std::ostream& os) const
{
  os << "Tpc_PolyTrackVertexContainerv1 with " << m_vertices.size()
     << " vertices and " << m_collision_x.size()
     << " collision vertices" << std::endl;
}

void Tpc_PolyTrackVertexContainerv1::Reset()
{
  for (Tpc_PolyTrackVertex* vtx : m_vertices) delete vtx;
  m_vertices.clear();
  m_collision_vertex_valid = 0;
  clear_collision_vertices();
}

int Tpc_PolyTrackVertexContainerv1::isValid() const
{
  return 1;
}

PHObject* Tpc_PolyTrackVertexContainerv1::CloneMe() const
{
  Tpc_PolyTrackVertexContainerv1* copy = new Tpc_PolyTrackVertexContainerv1();
  copy->m_collision_vertex_valid = m_collision_vertex_valid;
  copy->m_collision_x = m_collision_x;
  copy->m_collision_y = m_collision_y;
  copy->m_collision_z = m_collision_z;
  copy->m_collision_z_rms = m_collision_z_rms;
  copy->m_collision_ntracks = m_collision_ntracks;
  copy->m_collision_min_clusters = m_collision_min_clusters;
  for (Tpc_PolyTrackVertex* vtx : m_vertices)
  {
    if (vtx) copy->m_vertices.push_back(static_cast<Tpc_PolyTrackVertex*>(vtx->CloneMe()));
  }
  return copy;
}

const Tpc_PolyTrackVertex* Tpc_PolyTrackVertexContainerv1::get_vertex(unsigned int i) const
{
  if (i >= m_vertices.size()) return nullptr;
  return m_vertices[i];
}

Tpc_PolyTrackVertex* Tpc_PolyTrackVertexContainerv1::get_vertex(unsigned int i)
{
  if (i >= m_vertices.size()) return nullptr;
  return m_vertices[i];
}

double Tpc_PolyTrackVertexContainerv1::get_collision_x(unsigned int i) const
{
  if (i >= m_collision_x.size()) return 0.0;
  return m_collision_x[i];
}

double Tpc_PolyTrackVertexContainerv1::get_collision_y(unsigned int i) const
{
  if (i >= m_collision_y.size()) return 0.0;
  return m_collision_y[i];
}

double Tpc_PolyTrackVertexContainerv1::get_collision_z(unsigned int i) const
{
  if (i >= m_collision_z.size()) return 0.0;
  return m_collision_z[i];
}

double Tpc_PolyTrackVertexContainerv1::get_collision_z_rms(unsigned int i) const
{
  if (i >= m_collision_z_rms.size()) return 0.0;
  return m_collision_z_rms[i];
}

unsigned int Tpc_PolyTrackVertexContainerv1::get_collision_ntracks(unsigned int i) const
{
  if (i >= m_collision_ntracks.size()) return 0;
  return m_collision_ntracks[i];
}

void Tpc_PolyTrackVertexContainerv1::clear_collision_vertices()
{
  m_collision_x.clear();
  m_collision_y.clear();
  m_collision_z.clear();
  m_collision_z_rms.clear();
  m_collision_ntracks.clear();
}

void Tpc_PolyTrackVertexContainerv1::add_collision_vertex(double x, double y, double z, double z_rms, unsigned int ntracks)
{
  m_collision_x.push_back(x);
  m_collision_y.push_back(y);
  m_collision_z.push_back(z);
  m_collision_z_rms.push_back(z_rms);
  m_collision_ntracks.push_back(ntracks);
}
