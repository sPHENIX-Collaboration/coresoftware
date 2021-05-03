#ifndef TRACKBASE_TRKRCLUSTERCONTAINERV1_H
#define TRACKBASE_TRKRCLUSTERCONTAINERV1_H

/**
 * @file trackbase/TrkrClusterContainerv1.h
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
class TrkrClusterContainerv1 : public TrkrClusterContainer
{
  public:

  TrkrClusterContainerv1() = default;

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
  
  unsigned int max_layer = 57;

  unsigned int max_phisegment = 20;

  unsigned int max_zsegment = 15;

  Map m_clusmap[57][20][15];

  ClassDef(TrkrClusterContainerv1, 1)

};

#endif //TRACKBASE_TrkrClusterContainerv1_H
