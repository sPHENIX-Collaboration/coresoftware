#ifndef TRACKBASE_TRKRCLUSTERCONTAINERV2_H
#define TRACKBASE_TRKRCLUSTERCONTAINERV2_H

/**
 * @file trackbase/TrkrClusterContainerv2.h
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
class TrkrClusterContainerv2 : public TrkrClusterContainer
{
  public:

  TrkrClusterContainerv2() = default;

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

  virtual unsigned int size() const;

  private:
  
  unsigned int max_layer = 57;

  unsigned int max_phisegment = 20;

  unsigned int max_zsegment = 15;

  Map m_clusmap[57][20][15];

  ClassDef(TrkrClusterContainerv2, 1)

};

#endif //TRACKBASE_TrkrClusterContainerv2_H
