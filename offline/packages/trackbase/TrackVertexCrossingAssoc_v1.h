#ifndef TRACKBASE_TRACKVERTEXCROSSINGASSOCV1_H
#define TRACKBASE_TRACKVERTEXCROSSINGASSOCV1_H
/**

 * @file trackbase/TrackVertexCrossingAssoc_v1.h
 * @author Tony Frawley
 * @date March 2022
 * @brief Version 1 of class for associating clusters to the bunch crossing that created them
 */

#include "TrackVertexCrossingAssoc.h"
#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream
#include <map>
#include <set>
#include <utility>  // for pair

/**
 * @brief Class for associating clusters to the bunch crossing that created them
 *
 * Store the associations between clusters and beam crossings.
 */
class TrackVertexCrossingAssoc_v1 : public TrackVertexCrossingAssoc
{
 public:
  TrackVertexCrossingAssoc_v1() = default;

  void Reset() override;

  void identify(std::ostream &os = std::cout) const override;

  void addTrackAssoc(short int, unsigned int) override;
  void addVertexAssoc(short int, unsigned int) override;

  std::set<short int> getCrossings() const override;
  ConstRange getTracks(short int crossing) const override;
  ConstRange getVertices(short int crossing) const override;

  unsigned int sizeTracks(void) const override;
  unsigned int sizeVertices(void) const override;

 private:
  Map _track_assoc_map;
  Map _vertex_assoc_map;
  std::set<short int> _crossing_set;

  ClassDefOverride(TrackVertexCrossingAssoc_v1, 1);
};

#endif  // TRACKBASE_TRACKVERTEXCROSSINGASSOCV1_H
