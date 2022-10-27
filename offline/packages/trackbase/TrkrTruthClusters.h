#ifndef TRACKBASE_TRKRTRUTHCLUSTERS_H
#define TRACKBASE_TRKRTRUTHCLUSTERS_H

/**
 * @file trackbase/TrkrTruthClusters.h
 * @author D. Stewart
 * @date June 2022
 * @brief Keep track of mean and variance of phi, eta, and Z in clusters from truth hits
 *
 * Update October 2022 to remove pure virtual functions
 */

#include "TrkrCluster.h"
#include "TrkrClusterContainer.h"

#include <phool/PHObject.h>

#include <iostream>              // for cout, ostream
#include <map>
#include <utility>               // for pair

/**
 * @brief Association object for PHG4Clusters contributiong to TrkrTruthClusters
 *
 * Association object holding a vector of PHG4Clusters associated with given Truth Electrons
 * and associated streamed electrons
 *
 */
class TrkrTruthClusters : public PHObject
{
  public:
  //! typedefs for convenience 
  using Entry         = std::pair<short, TrkrDefs::cluskey>; 
  // Provide ways to compare Entry & Entry,short,pair<short,short>
  using Vector        = std::vector<Entry>; // note that the entries are set in order of trkid, and then cluskey ordered by layerid
  using ConstIterator = Vector::const_iterator; 
  using ConstRange    = std::pair<ConstIterator,ConstIterator>;

  using Iterator      = Vector::iterator; // not inplemented in TrkrTruthClustersv1
  using Range         = std::pair<Iterator,Iterator>; // not implemented in TrkrTruthClusterv1

  void Reset() override {};

  virtual ConstRange  getClusters       (short trackid=-1)   const; // will only iterate over range of trackid (defaults to all tracks if none provided)
  virtual std::vector<short> getTrkIds  (short /*layer*/=-1) const { return {}; }; // default to values for all layers
  virtual std::vector<short> getLayers  (short /*trkid*/=-1) const { return {}; }; // default to values for all tracks
  virtual bool        hasTrkId          (short /*trkid*/)    const { return false; };
  virtual bool        hasTrkIdLayer     (short /*trkid*/, short /*layer*/) const { return false; };
  virtual bool        hasLayer          (short /*layer*/=-1) const { return false; };
  virtual void        addTruthCluster   (short /*trkid*/, TrkrDefs::cluskey) {};// 
  virtual std::ostream& print_clusters (TrkrClusterContainer*, std::ostream &/*os*/ = std::cout) const { return std::cout; };
  virtual TrkrDefs::cluskey   getCluskey(short /*trackid*/, short /*layer*/) const { return 0; };// return 0 if not present

  // Method's required for STL to search (and sort) vectors of Entry (std::pair<short,TrkrDefs::cluskey>)
  // Sort by Track ID
  struct CompTrkId {
    bool operator()(const Entry& lhs, const short  rhs) const { return lhs.first < rhs; }
    bool operator()(const short  lhs, const Entry& rhs) const { return lhs < rhs.first; }
  };
  // Sort by TPC Layer
  struct CompLayer {
    bool operator()(const Entry& lhs, const short  rhs) const 
    { return TrkrDefs::getLayer(lhs.second) < rhs; }
    bool operator()(const short  lhs, const Entry& rhs) const 
    { return lhs < TrkrDefs::getLayer(rhs.second); }
  };
  // Sort by Track ID AND TPC Layer
  struct CompTrkIdLayer {
    bool operator()(const Entry& lhs, const std::pair<short, short> rhs) const 
    { 
      if (lhs.first == rhs.first) { return TrkrDefs::getLayer(lhs.second) < rhs.second; }
      else return lhs.first < rhs.first;
    }
    bool operator()(const std::pair<short,short> lhs, const Entry& rhs) const 
    { 
      if (lhs.first == rhs.first) { return rhs.second < TrkrDefs::getLayer(rhs.second); }
      else return lhs.first < rhs.first;
    }
  };

  void identify(std::ostream &/*os*/) const override {};

  protected:
  //! ctor
  TrkrTruthClusters() = default;

  private:
  ClassDefOverride(TrkrTruthClusters, 1);
};

#endif //TRACKBASE_TRKRTRUTHCLUSTERS_H
