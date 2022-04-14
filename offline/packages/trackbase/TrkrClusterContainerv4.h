#ifndef TRACKBASE_TRKRCLUSTERCONTAINERV4_H
#define TRACKBASE_TRKRCLUSTERCONTAINERV4_H

/**
 * @file trackbase/TrkrClusterContainerv4.h
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
class TrkrClusterContainerv4 : public TrkrClusterContainer
{
  public:

  TrkrClusterContainerv4() = default;

  void Reset() override;

  void identify(std::ostream &os = std::cout) const override;

  void addCluster(TrkrCluster*) override;

  void addClusterSpecifyKey(const TrkrDefs::cluskey, TrkrCluster*) override;

  void removeCluster(TrkrDefs::cluskey) override;

  void removeCluster(TrkrCluster*) override;

  ConstRange getClusters(TrkrDefs::hitsetkey) override;

  TrkrCluster* findCluster(TrkrDefs::cluskey) const override;

  unsigned int size(void) const override;

  private:
  
  //! the actual container
  std::map<TrkrDefs::hitsetkey, Vector> m_clusmap;

  // temporary map
  Map m_tmpmap;
  
  ClassDefOverride(TrkrClusterContainerv4, 1)
};

#endif //TRACKBASE_TrkrClusterContainerv4_H
