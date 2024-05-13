/**
 * @file trackbase/TrackVertexCrossingAssoc.h
 * @author Tony Frawley
 * @date April 2024
 * @brief Base class for associating tracks and vertices to the bunch crossing they originated from
 */
#ifndef TRACKBASE_TRACKVERTEXCROSSINGASSOC_H
#define TRACKBASE_TRACKVERTEXCROSSINGASSOC_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <climits>
#include <iostream>  // for cout, ostream
#include <map>
#include <set>
#include <utility>  // for pair

/**
 * @brief Base class for associating clusters to the hits that went into them
 *
 * Store the associations between clusters and the hits that went into them.
 */
class TrackVertexCrossingAssoc : public PHObject
{
 public:
  using Map = std::multimap<short int, unsigned int>;
  using ConstIterator = Map::const_iterator;
  using ConstRange = std::pair<Map::const_iterator, Map::const_iterator>;

  void Reset() override;

  /**
   * @brief Add association between tracks, vertices and crossings
   * @param[in] ckey track ID, vertex ID
   * @param[in] bunch crossing number
   */
  virtual std::set<short int> getCrossings() const = 0;
  virtual void addTrackAssoc(short int crossing, unsigned int trackid) = 0;
  virtual void addVertexAssoc(short int crossing, unsigned int vertexid) = 0;

  virtual ConstRange getTracks(short int crossing) const = 0;
  virtual ConstRange getVertices(short int crossing) const = 0;

  virtual unsigned int sizeTracks() const { return 0; }
  virtual unsigned int sizeVertices() const { return 0; }

 protected:
  TrackVertexCrossingAssoc() = default;

 private:
  ClassDefOverride(TrackVertexCrossingAssoc, 1);
};

#endif  // TRACKBASE_TRACKVERTREXCROSSINGASSOC_H
