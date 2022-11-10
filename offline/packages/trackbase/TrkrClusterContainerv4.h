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

  void identify(std::ostream& os = std::cout) const override;

  void addClusterSpecifyKey(const TrkrDefs::cluskey, TrkrCluster*) override;

  void removeCluster(TrkrDefs::cluskey) override;

  ConstRange getClusters() const override;  // deprecated

  ConstRange getClusters(TrkrDefs::hitsetkey) override;

  TrkrCluster* findCluster(TrkrDefs::cluskey) const override;

  HitSetKeyList getHitSetKeys() const override;

  HitSetKeyList getHitSetKeys(const TrkrDefs::TrkrId) const override;

  HitSetKeyList getHitSetKeys(const TrkrDefs::TrkrId, const uint8_t /* layer */) const override;

  unsigned int size(void) const override;

 private:
  /// convenient alias
  using Vector = std::vector<TrkrCluster*>;

  /// the actual container
  std::map<TrkrDefs::hitsetkey, Vector> m_clusmap;

  /// temporary map
  /**
   * the map is transient. It must not be written to the output.
   * To do this one adds //! after the declaration
   * see https://root.cern.ch/root/htmldoc/guides/users-guide/InputOutput.html for details
   */
  Map m_tmpmap;  //! transient. The temporary map does not get written to the output

  ClassDefOverride(TrkrClusterContainerv4, 1)
};

#endif  // TRACKBASE_TrkrClusterContainerv4_H
