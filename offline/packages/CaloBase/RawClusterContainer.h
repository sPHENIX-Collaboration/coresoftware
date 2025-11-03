#ifndef CALOBASE_RAWCLUSTERCONTAINER_H
#define CALOBASE_RAWCLUSTERCONTAINER_H

#include "RawClusterDefs.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <utility>

class RawCluster;

class RawClusterContainer : public PHObject
{
 public:
  typedef std::map<RawClusterDefs::keytype, RawCluster *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  RawClusterContainer() = default;
  ~RawClusterContainer() override = default;

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream &os = std::cout) const override;

  ConstIterator AddCluster(RawCluster *rawcluster);

  RawCluster *getCluster(const RawClusterDefs::keytype key);
  const RawCluster *getCluster(const RawClusterDefs::keytype key) const;

  //! return all clusters
  ConstRange getClusters() const;
  Range getClusters();
  const Map &getClustersMap() const { return _clusters; }
  Map &getClustersMap() { return _clusters; }

  unsigned int size() const { return _clusters.size(); }
  double getTotalEdep() const;

 protected:
  Map _clusters;

  ClassDefOverride(RawClusterContainer, 1)
};

#endif
