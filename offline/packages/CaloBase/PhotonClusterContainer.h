#ifndef CALOBASE_PHOTONCLUSTERCONTAINER_H
#define CALOBASE_PHOTONCLUSTERCONTAINER_H

#include "RawClusterDefs.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <utility>

class PhotonCluster; // interface

// Container owning photon cluster objects (PhotonClusterv1 or derivatives)
class PhotonClusterContainer : public PHObject
{
 public:
  typedef std::map<RawClusterDefs::keytype, PhotonCluster*> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  PhotonClusterContainer() = default;
  ~PhotonClusterContainer() override {}

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream& os = std::cout) const override;

  ConstIterator AddCluster(PhotonCluster* clus); // takes ownership

  PhotonCluster* getCluster(const RawClusterDefs::keytype id);
  const PhotonCluster* getCluster(const RawClusterDefs::keytype id) const;

  ConstRange getClusters() const;
  Range getClusters();
  const Map& getClustersMap() const { return m_clusters; }
  Map& getClustersMap() { return m_clusters; }

  unsigned int size() const { return m_clusters.size(); }

 protected:
  Map m_clusters;

  ClassDefOverride(PhotonClusterContainer, 1)
};

#endif
