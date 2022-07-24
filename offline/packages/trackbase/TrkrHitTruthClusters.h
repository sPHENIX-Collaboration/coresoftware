#ifndef TRACKBASE_TRKRHITTRUTHCLUSTERS_H
#define TRACKBASE_TRKRHITTRUTHCLUSTERS_H

/**
 * @file trackbase/TrkrHitTruthClusters.h
 * @author D. Stewart
 * @date June 2022
 * @brief Keep track of mean and variance of phi, eta, and Z in clusters from truth hits
 */

#include "TrkrDefs.h"

#include <g4main/PHG4HitDefs.h>
#include <phool/PHObject.h>

#include <iostream>              // for cout, ostream
#include <map>
#include <utility>               // for pair

/**
 * @brief Association object for PHG4Clusters contributiong to TrkrHitTruthClusters
 *
 * Association object holding a vector of PHG4Clusters associated with given Truth Electrons
 * and associated streamed electrons
 *
 */
struct EnergyCentroid {
  float phi_ave, phi_stdev, z_ave, z_stdev, sum_E;
  EnergyCentroid( std::array<float,5> input) :
    phi_ave   {input[0]},
    phi_stdev {input[1]},
    z_ave     {input[2]},
    z_stdev   {input[3]},
    sum_E     {input[4]} {};
  EnergyCentroid() : phi_ave{0}, phi_stdev{0}, z_ave{0}, z_stdev{0}, sum_E{0.} {};
  void set_values( std::array<float,5> input) {
    phi_ave   =input[0];
    phi_stdev =input[1];
    z_ave     =input[2];
    z_stdev   =input[3];
    sum_E     =input[4];
  }
  ~EnergyCentroid() {};
};

class TrkrHitTruthClusters : public PHObject
{
  public:
  
  //! typedefs for convenience 
  static const int N_GEM_LAYERS = 55;
  using CentroidsFor1Track = std::array<EnergyCentroid, N_GEM_LAYERS>;
  using MMap = std::map< int /*track-id*/, CentroidsFor1Track>; /*55 possible energy centroids*/
  using Iterator = MMap::iterator;
  using ConstIterator = MMap::const_iterator;
  using Range = std::pair<Iterator, Iterator>;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;

  void Reset() override
  {}

  virtual void print_clusters (std::ostream &/*os*/ = std::cout) const 
  {}

  virtual CentroidsFor1Track& add_track_centroids(const int track_id) =0;

  //virtual void removeAssoc(const TrkrDefs::hitsetkey /*hitsetkey*/, const TrkrDefs::hitkey /*hitkey*/)
  //{}

  protected:
  //! ctor
  TrkrHitTruthClusters() = default;

  private:
  ClassDefOverride(TrkrHitTruthClusters, 1);
};

#endif //TRACKBASE_TRKRHITTRUTHCLUSTERS_H
