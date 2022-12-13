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

  void Reset() override;

  void identify(std::ostream& os = std::cout) const override;

  void addClusterSpecifyKey(const TrkrDefs::cluskey, TrkrCluster*) override;

  void removeCluster(TrkrDefs::cluskey) override;

  ConstRange getClusters(TrkrDefs::hitsetkey) override;

  TrkrCluster* findCluster(TrkrDefs::cluskey) const override;

  unsigned int size() const override;

 private:
  unsigned int max_layer = 57;

  unsigned int max_phisegment = 20;

  unsigned int max_zsegment = 15;

  Map m_clusmap[57][20][15];

  ClassDefOverride(TrkrClusterContainerv2, 1)
};

#endif  // TRACKBASE_TrkrClusterContainerv2_H
