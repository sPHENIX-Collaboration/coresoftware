/**
 * @file trackbase/TrackVertexCrossingAssoc_v1.cc
 * @author Tony Frawley
 * @date April 2024
 * @brief TrackVertexCrossingAssoc implementation
 */

#include "TrackVertexCrossingAssoc_v1.h"
#include "TrkrDefs.h"

#include <ostream>  // for operator<<, endl, basic_ostream, ostream, basic_o...

namespace
{
  TrackVertexCrossingAssoc_v1::Map dummy_map;
}

//_________________________________________________________________________
void TrackVertexCrossingAssoc_v1::Reset()
{
  // delete all entries
  Map empty_trackmap;
  empty_trackmap.swap(_track_assoc_map);
  empty_trackmap.clear();
  _track_assoc_map.clear();

  Map empty_vertexmap;
  empty_vertexmap.swap(_vertex_assoc_map);
  empty_vertexmap.clear();
  _vertex_assoc_map.clear();
}

//_________________________________________________________________________
void TrackVertexCrossingAssoc_v1::identify(std::ostream& os) const
{
  std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator iter;
  os << "-----TrackVertexCrossingAssoc_v1-----" << std::endl;
  os << "Number of vertex associations: " << sizeVertices()  << "Number of track associations: " << sizeTracks() << std::endl;
  for (const auto& cross : _crossing_set)
    { 
      auto vtxit = _vertex_assoc_map.equal_range(cross);
      os << "Crossing " << cross << std::endl;
      os << "   vertex ID's:  ";
      for (auto itr = vtxit.first; itr != vtxit.second; ++itr)
	{
	  os << itr->second << "  ";
	} 
      os << std::endl;

      auto trit = _track_assoc_map.equal_range(cross);
      os << "   track ID's:  ";
      for (auto itr = trit.first; itr != trit.second; ++itr)
	{
	  os << itr->second << "  ";
	} 
      os << std::endl;

    }

      /*
     for (const auto& map_pair : _vertex_assoc_map)
  {
    os << "crossing " << map_pair.first << std::dec
       << "  vertexid: " << map_pair.second; 
  }
  os << "------------------------------" << std::endl;
      */
  return;
}

//_________________________________________________________________________
void TrackVertexCrossingAssoc_v1::addTrackAssoc(short int crossing, unsigned int trackid)
{
  // insert association
  _track_assoc_map.insert(std::make_pair(crossing, trackid));
  _crossing_set.insert(crossing);
}

//_________________________________________________________________________
void TrackVertexCrossingAssoc_v1::addVertexAssoc(short int crossing, unsigned int vertexid)
{
  // insert association
  _vertex_assoc_map.insert(std::make_pair(crossing, vertexid));
  _crossing_set.insert(crossing);
}

//_________________________________________________________________________
std::set<short int> TrackVertexCrossingAssoc_v1::getCrossings() const
{
  return _crossing_set;
}

//_________________________________________________________________________
TrackVertexCrossingAssoc_v1::ConstRange TrackVertexCrossingAssoc_v1::getTracks(short int crossing) const
{
  const auto range = _track_assoc_map.equal_range(crossing);
  return range;
}

//_________________________________________________________________________
TrackVertexCrossingAssoc_v1::ConstRange TrackVertexCrossingAssoc_v1::getVertices(short int crossing) const
{
  const auto range = _vertex_assoc_map.equal_range(crossing);
  return range;
}

//_________________________________________________________________________
unsigned int TrackVertexCrossingAssoc_v1::sizeTracks() const
{
  unsigned int size = 0;
  size = _track_assoc_map.size();

  return size;
}

//_________________________________________________________________________
unsigned int TrackVertexCrossingAssoc_v1::sizeVertices() const
{
  unsigned int size = 0;
  size = _vertex_assoc_map.size();

  return size;
}

