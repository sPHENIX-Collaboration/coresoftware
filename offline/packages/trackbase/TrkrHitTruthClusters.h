#ifndef TRACKBASE_TRKRHITTRUTHCLUSTERS_H
#define TRACKBASE_TRKRHITTRUTHCLUSTERS_H

/**
 * @file trackbase/TrkrHitTruthClusters.h
 * @author D. Stewart
 * @date June 2022
 * @brief Keep track of mean and variance of phi, eta, and Z in clusters from truth hits
 */

#include "TrkrCluster.h"
#include "MapToPadPlanePassData.h"
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
class TrkrHitTruthClusters : public PHObject
{
  
  public:
  
  //! typedefs for convenience 
  using Key           = std::pair<short,short>;
  using Entry         = std::pair<Key, TrkrCluster*>;

  using Vector        = std::vector<Entry>;

  using ConstIterator = Vector::const_iterator;
  using ConstRange    = std::pair<ConstIterator,ConstIterator>;

  using Iterator      = Vector::iterator; // not inplemented with an accessor
  using Range         = std::pair<Iterator,Iterator>; // not implemented with an accessor

  void Reset() override {};

  virtual std::vector<short> getTrkIds   (short layer=-1) const =0; // default to values for all layers
  virtual std::vector<short> getLayerIds (short trkid=-1) const =0; // default to values for all tracks
  virtual bool        hasTrkId           (short trkid)              const =0;
  virtual bool        hasTrkIdLayerId    (short trkid, short layer) const =0;
  virtual bool        hasLayerId         (short layer=-1)           const =0;
  virtual ConstRange  getClusters        (short trackid=-1)         const =0; // will only iterate over range of trackid (if provided)
  virtual void        addTruthCluster    (short trkid, MapToPadPlanePassData& hit_data) =0;

  protected:
  //! ctor
  TrkrHitTruthClusters() = default;

  private:
  ClassDefOverride(TrkrHitTruthClusters, 1);
};

#endif //TRACKBASE_TRKRHITTRUTHCLUSTERS_H
