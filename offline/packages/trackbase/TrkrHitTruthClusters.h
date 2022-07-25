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
  short layer_id;
  float phi_ave, phi_stdev, z_ave, z_stdev, sum_E;
  EnergyCentroid( short _layer_id, std::array<float,5> input) :
    layer_id  {_layer_id},
    phi_ave   {input[0]},
    phi_stdev {input[1]},
    z_ave     {input[2]},
    z_stdev   {input[3]},
    sum_E     {input[4]} {};
  EnergyCentroid() : layer_id{0}, phi_ave{0}, phi_stdev{0}, z_ave{0}, z_stdev{0}, sum_E{0.} {};
  void set_values( short _layer_id, std::array<float,5> input) {
    layer_id = _layer_id;
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
  /* using CentroidsFor1Track = std::array<EnergyCentroid, N_GEM_LAYERS>; */
  using VecEC = std::vector<EnergyCentroid>; // the data 
  /* using VecECiter = VecEC::iterator; */
  /* using ConstVecECiter = VecEC::const_iterator; */
  using MMap = std::map<short, VecEC>; // trk-id, first and last index in CentroidVector
  using Iterator = MMap::iterator;
  using ConstIterator = MMap::const_iterator;
  using Range = std::pair<Iterator, Iterator>;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;

  void Reset() override
  {}

  virtual void print_clusters (std::ostream &/*os*/ = std::cout) const 
  {}

  virtual VecEC& get_new_centroids_vec (short track_id) =0;

  //virtual void removeAssoc(const TrkrDefs::hitsetkey /*hitsetkey*/, const TrkrDefs::hitkey /*hitkey*/)
  //{}

  protected:
  //! ctor
  TrkrHitTruthClusters() = default;

  private:
  ClassDefOverride(TrkrHitTruthClusters, 1);
};

#endif //TRACKBASE_TRKRHITTRUTHCLUSTERS_H
