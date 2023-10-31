#include "GlobalVertexv2.h"

#include <cmath>

GlobalVertexv2::GlobalVertexv2()
  : _id(std::numeric_limits<unsigned int>::max())
  , _bco(std::numeric_limits<unsigned int>::max())
{
}

GlobalVertexv2::GlobalVertexv2(const unsigned int id)
  : _id(id)
  ,_bco(std::numeric_limits<unsigned int>::max())
{
}

void GlobalVertexv2::identify(std::ostream& os) const
{
  os << "---GlobalVertexv2-----------------------" << std::endl;

  os << " list of vtx ids: " << std::endl;
  for(ConstVertexIter iter = begin_vertexes(); iter != end_vertexes(); ++iter)
    {
      os << "  Vertex type " << iter->first << " has " << iter->second.size() 
	 << " vertices associated to it" << std::endl;
      for(auto& vertex : iter->second)
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
  if( it == _vtxs.end())
    {
      VertexVector vector;
      vector.push_back(vertex);
      _vtxs.insert(std::make_pair(type, vector));
      return;
    }

  it->second.push_back(vertex);


}
size_t GlobalVertexv2::count_vtxs(GlobalVertex::VTXTYPE type) const
{
  auto it = _vtxs.find(type);
  if(it == _vtxs.end()) 
    {
      return 0;
    }

  return it->second.size();
}
