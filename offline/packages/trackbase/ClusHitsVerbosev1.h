#ifndef G4TRACKING_CLUSHITSVERBOSEV1_H
#define G4TRACKING_CLUSHITSVERBOSEV1_H

/**
 * @file trackbase/ClusHitsVerbosev1.h
 * @author D. Stewart
 * @date May 2023
 * @brief Ojbect to show hit locations in given clusters
 */

#include "ClusHitsVerbose.h"

/**
 * @brief Cluster container object
 */
class ClusHitsVerbosev1 : public ClusHitsVerbose
{
 public:
  ClusHitsVerbosev1() = default;

  void Reset() override;

  bool hasClusKey(TrkrDefs::cluskey) const override;
  Vector& phiBins(TrkrDefs::cluskey) override;
  Vector& zBins(TrkrDefs::cluskey) override;
  Vector& phiCutBins(TrkrDefs::cluskey) override;
  Vector& zCutBins(TrkrDefs::cluskey) override;
  Map& getMap() override { return m_data; };

  PairVector phiBins_pvecIE(TrkrDefs::cluskey) override;
  PairVector phiCutBins_pvecIE(TrkrDefs::cluskey) override;
  PairVector zBins_pvecIE(TrkrDefs::cluskey) override;
  PairVector zCutBins_pvecIE(TrkrDefs::cluskey) override;

  void addPhiHit(int _i, int _v) override { m_stage_phi.push_back({_i, _v}); };
  void addPhiCutHit(int _i, int _v) override { m_stage_phiCut.push_back({_i, _v}); };
  void addZHit(int _i, int _v) override { m_stage_z.push_back({_i, _v}); };
  void addZCutHit(int _i, int _v) override { m_stage_zCut.push_back({_i, _v}); };
  void push_hits(TrkrDefs::cluskey) override;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "ClusHitsVerbosev1 class" << std::endl;
  };

 private:
  // the data
  Map m_data{};
  Vector m_stage_phi{};
  Vector m_stage_z{};
  Vector m_stage_phiCut{};
  Vector m_stage_zCut{};

  PairVector pvecIE(TrkrDefs::cluskey, int which);
  Vector& vecBins(TrkrDefs::cluskey, int which);

  ClassDefOverride(ClusHitsVerbosev1, 1)
};

#endif  // G4TRACKING_CLUSHITSVERBOSEV1_H
