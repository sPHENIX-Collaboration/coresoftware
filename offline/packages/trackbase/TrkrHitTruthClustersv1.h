#ifndef TRACKBASE_TRKRHITTRUTHCLUSTERSV1_H
#define TRACKBASE_TRKRHITTRUTHCLUSTERSV1_H

/**
 * @file trackbase/TrkrHitTruthClustersv1.h
 * @author D. STEWART
 * @date April 2022
 */
#include "TrkrHitTruthClusters.h"

/**
 * @brief Association object for PHG4Cells contributiong to TrkrHits
 *
 * Association object holding a multimap of PHG4Cells associated with a given TrkrHit
 */
class TrkrHitTruthClustersv1 : public TrkrHitTruthClusters
{
  
  public:

  TrkrHitTruthClustersv1() = default;
  void Reset() override;

  ConstRange getClusters         (short trkid=-1)           const override;
  std::vector<short> getTrkIds   (short layer=-1)           const override;
  std::vector<short> getLayerIds (short trkid=-1)           const override;
  bool       hasTrkId            (short)                    const override;
  bool       hasTrkIdLayerId     (short trkid, short layer) const override;
  bool       hasLayerId          (short layer=-1)           const override;
  void       addTruthCluster     (short trkid, MapToPadPlanePassData& hit_data) override;

  private:
  Vector m_data;
  
  ClassDefOverride(TrkrHitTruthClustersv1, 1);

};

#endif //TRACKBASE_TRKRHITTRUTHCLUSTERS_H

