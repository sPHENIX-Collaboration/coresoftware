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

  virtual void Reset();

  virtual void identify(std::ostream &os = std::cout) const;

  virtual ConstIterator addCluster(TrkrCluster*);

  virtual ConstIterator addClusterSpecifyKey(const TrkrDefs::cluskey, TrkrCluster*);

  virtual void removeCluster(TrkrDefs::cluskey);

  virtual void removeCluster(TrkrCluster*);

  virtual Iterator findOrAddCluster(TrkrDefs::cluskey);

  virtual ConstRange getClusters(TrkrDefs::hitsetkey) const;

  virtual Map* getClusterMap(TrkrDefs::hitsetkey);

  virtual TrkrCluster* findCluster(TrkrDefs::cluskey) const;

  virtual unsigned int size(void) const;

  private:
  
  std::map<TrkrDefs::hitsetkey, Map> m_clusmap;

  ClassDef(TrkrClusterContainerv3, 1)

};

#endif //TRACKBASE_TrkrClusterContainerv3_H
