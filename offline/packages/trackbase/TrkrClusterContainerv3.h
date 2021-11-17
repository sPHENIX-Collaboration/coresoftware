#ifndef TRACKBASE_TRKRCLUSTERCONTAINERV3_H
#define TRACKBASE_TRKRCLUSTERCONTAINERV3_H

/**
 * @file trackbase/TrkrClusterContainerv3.h
 * @author D. McGlinchey, Hugo Pereira Da Costa
 * @date June 2018
 * @brief Cluster container object
 */

#include "TrkrClusterContainer.h"

#include <phool/PHObject.h>

class TrkrCluster;

/**
 * @brief Cluster container object
 */
class TrkrClusterContainerv3 : public TrkrClusterContainer
{
  public:

  TrkrClusterContainerv3() = default;

  void Reset() override;

  void identify(std::ostream &os = std::cout) const override;

  ConstIterator addCluster(TrkrCluster*) override;

  ConstIterator addClusterSpecifyKey(const TrkrDefs::cluskey, TrkrCluster*) override;

  void removeCluster(TrkrDefs::cluskey) override;

  void removeCluster(TrkrCluster*) override;

  Iterator findOrAddCluster(TrkrDefs::cluskey) override;

  ConstRange getClusters(TrkrDefs::hitsetkey) const override;

  Map* getClusterMap(TrkrDefs::hitsetkey) override;

  TrkrCluster* findCluster(TrkrDefs::cluskey) const override;

  unsigned int size(void) const override;

  private:
  
  std::map<TrkrDefs::hitsetkey, Map> m_clusmap;

  ClassDefOverride(TrkrClusterContainerv3, 1)

};

#endif //TRACKBASE_TrkrClusterContainerv3_H
