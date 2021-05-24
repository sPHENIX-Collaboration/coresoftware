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

  virtual void Reset() override;

  virtual void identify(std::ostream &os = std::cout) const override;

  virtual ConstIterator addCluster(TrkrCluster*) override;

  virtual ConstIterator addClusterSpecifyKey(const TrkrDefs::cluskey, TrkrCluster*) override;

  virtual void removeCluster(TrkrDefs::cluskey) override;

  virtual void removeCluster(TrkrCluster*) override;

  virtual Iterator findOrAddCluster(TrkrDefs::cluskey) override;

  virtual ConstRange getClusters(TrkrDefs::hitsetkey) const override;

  virtual Map* getClusterMap(TrkrDefs::hitsetkey) override;

  virtual TrkrCluster* findCluster(TrkrDefs::cluskey) const override;

  virtual unsigned int size(void) const override;

  private:
  
  std::map<TrkrDefs::hitsetkey, Map> m_clusmap;

  ClassDefOverride(TrkrClusterContainerv3, 1)

};

#endif //TRACKBASE_TrkrClusterContainerv3_H
