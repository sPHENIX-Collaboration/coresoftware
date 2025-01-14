#include "GlobalVertexv2.h"

GlobalVertexv2::GlobalVertexv2(const unsigned int id)
  : _id(id)
{
}

GlobalVertexv2::~GlobalVertexv2()
{
  GlobalVertexv2::Reset();
}

void GlobalVertexv2::Reset()
{
  for (auto & _vtx : _vtxs)
  {
    for (auto vertex : _vtx.second)
    {
      delete vertex;
    }
  }
  _vtxs.clear();
}

void GlobalVertexv2::identify(std::ostream& os) const
{
  os << "---GlobalVertexv2-----------------------" << std::endl;

  os << " list of vtx ids: " << std::endl;
  for (ConstVertexIter iter = begin_vertexes(); iter != end_vertexes(); ++iter)
  {
    os << "  Vertex type " << iter->first << " has " << iter->second.size()
       << " vertices associated to it" << std::endl;
    for (auto& vertex : iter->second)
    {
      vertex->identify();
    }
  }

  os << "-----------------------------------------------" << std::endl;

  return;
}

int GlobalVertexv2::isValid() const
{
  if (_vtxs.empty())
  {
    return 0;
  }
  return 1;
}
void GlobalVertexv2::insert_vtx(GlobalVertex::VTXTYPE type, const Vertex* vertex)
{
  auto it = _vtxs.find(type);
  if (it == _vtxs.end())
  {
    VertexVector vector;
    vector.push_back(vertex);
    _vtxs.insert(std::make_pair(type, vector));
    return;
  }

  it->second.push_back(vertex);
}
void GlobalVertexv2::clone_insert_vtx(GlobalVertex::VTXTYPE type, const Vertex* vertex)
{
  auto it = _vtxs.find(type);
  Vertex *clone = dynamic_cast<Vertex *> (vertex->CloneMe());
  if (it == _vtxs.end())
  {
    VertexVector vector;
    vector.push_back(clone);
    _vtxs.insert(std::make_pair(type, vector));
    return;
  }

  it->second.push_back(clone);
}
size_t GlobalVertexv2::count_vtxs(GlobalVertex::VTXTYPE type) const
{
  auto it = _vtxs.find(type);
  if (it == _vtxs.end())
  {
    return 0;
  }

  return it->second.size();
}

float GlobalVertexv2::get_t() const
{
  auto it = _vtxs.find(GlobalVertex::VTXTYPE::MBD);
  if (it == _vtxs.end())
  {
    return std::numeric_limits<float>::quiet_NaN();
  }

  return it->second[0]->get_t();
}

float GlobalVertexv2::get_t_err() const
{
  auto it = _vtxs.find(GlobalVertex::VTXTYPE::MBD);
  if (it == _vtxs.end())
  {
    return std::numeric_limits<float>::quiet_NaN();
  }

  return it->second[0]->get_t_err();
}
float GlobalVertexv2::get_x() const
{
  return get_position(0);
}
float GlobalVertexv2::get_y() const
{
  return get_position(1);
}
float GlobalVertexv2::get_z() const
{
  return get_position(2);
}

float GlobalVertexv2::get_position(unsigned int coor) const
{
  auto svtxit = _vtxs.find(GlobalVertex::VTXTYPE::SVTX);
  if (svtxit == _vtxs.end())
  {
    auto mbdit = _vtxs.find(GlobalVertex::VTXTYPE::MBD);
    if (mbdit == _vtxs.end())
    {
      return std::numeric_limits<float>::quiet_NaN();
    }
    return mbdit->second[0]->get_position(coor);
  }

  GlobalVertex::VertexVector trackvertices = svtxit->second;
  size_t mosttracks = 0;
  float pos = std::numeric_limits<float>::quiet_NaN();
  for (auto vertex : trackvertices)
  {
    if (vertex->size_tracks() > mosttracks)
    {
      mosttracks = vertex->size_tracks();
      pos = vertex->get_position(coor);
    }
  }

  return pos;
}

float GlobalVertexv2::get_chisq() const
{
  auto svtxit = _vtxs.find(GlobalVertex::VTXTYPE::SVTX);
  if (svtxit == _vtxs.end())
  {
    return std::numeric_limits<float>::quiet_NaN();
  }

  GlobalVertex::VertexVector trackvertices = svtxit->second;
  size_t mosttracks = 0;
  float chisq = std::numeric_limits<float>::quiet_NaN();
  for (auto vertex : trackvertices)
  {
    if (vertex->size_tracks() > mosttracks)
    {
      mosttracks = vertex->size_tracks();
      chisq = vertex->get_chisq();
    }
  }

  return chisq;
}

unsigned int GlobalVertexv2::get_ndof() const
{
  auto svtxit = _vtxs.find(GlobalVertex::VTXTYPE::SVTX);
  if (svtxit == _vtxs.end())
  {
    return std::numeric_limits<unsigned int>::max();
  }

  GlobalVertex::VertexVector trackvertices = svtxit->second;
  size_t mosttracks = 0;
  unsigned int ndf = std::numeric_limits<unsigned int>::max();
  for (auto vertex : trackvertices)
  {
    if (vertex->size_tracks() > mosttracks)
    {
      mosttracks = vertex->size_tracks();
      ndf = vertex->get_ndof();
    }
  }

  return ndf;
}

float GlobalVertexv2::get_error(unsigned int i, unsigned int j) const
{
  auto svtxit = _vtxs.find(GlobalVertex::VTXTYPE::SVTX);
  if (svtxit == _vtxs.end())
  {
    auto mbdit = _vtxs.find(GlobalVertex::VTXTYPE::MBD);
    if (mbdit == _vtxs.end())
    {
      return std::numeric_limits<float>::quiet_NaN();
    }
    // MBD only has z error defined
    if (i == 2 && j == 2)
    {
      return mbdit->second[0]->get_z_err();
    }
    else
    {
      return std::numeric_limits<float>::quiet_NaN();
    }
  }

  GlobalVertex::VertexVector trackvertices = svtxit->second;
  size_t mosttracks = 0;
  float err = std::numeric_limits<float>::quiet_NaN();
  for (auto vertex : trackvertices)
  {
    if (vertex->size_tracks() > mosttracks)
    {
      mosttracks = vertex->size_tracks();
      err = vertex->get_error(i, j);
    }
  }

  return err;
}
