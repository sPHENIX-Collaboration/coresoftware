#ifndef TRACKBASE_TRKRHITTRUTHCLUSTERSV1_H
#define TRACKBASE_TRKRHITTRUTHCLUSTERSV1_H

/**
 * @file trackbase/TrkrTruthClustersv1.h
 * @author D. STEWART
 * @date April 2022
 */
#include "TrkrTruthClusters.h"

/**
 * @brief Association object for PHG4Cells contributiong to Trkrs
 *
 * Association object holding a multimap of PHG4Cells associated with a given Trkr
 */
class TrkrTruthClustersv1 : public TrkrTruthClusters
{
  
  public:

  TrkrTruthClustersv1() = default;
  void Reset() override {};

  ConstRange getClusters       (short trkid=-1)           const override;
  std::vector<short> getTrkIds (short layer=-1)           const override;
  std::vector<short> getLayers (short trkid=-1)           const override;
  bool       hasTrkId          (short)                    const override;
  bool       hasTrkIdLayer     (short trkid, short layer) const override;
  TrkrDefs::cluskey getCluskey(short trackid, short layer) const override;
  bool       hasLayer          (short layer=-1)           const override;
  void       addTruthCluster   (short trkid, TrkrDefs::cluskey) override;
  std::ostream& print_clusters (TrkrClusterContainer*, std::ostream &os = std::cout) const override;

  private:
  Vector m_data;
  
  ClassDefOverride(TrkrTruthClustersv1, 1);
};

#endif //TRACKBASE_TRKRHITTRUTHCLUSTERS_H

