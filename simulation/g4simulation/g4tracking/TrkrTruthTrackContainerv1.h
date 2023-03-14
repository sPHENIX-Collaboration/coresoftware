#ifndef G4TRACKING_TRUTHTRACKCONTAINERV1_H
#define G4TRACKING_TRUTHTRACKCONTAINERV1_H

/**
 * @file trackbase/TrkrTruthTrackContainerv1.h
 * @author D. Stewart
 * @date September 2022
 * @brief TrkrTruthTrack container object
 */

#include "TrkrTruthTrackContainer.h"
class TrkrTruthTrack;
class PHG4TruthInfoContainer;

/**
 * @brief Cluster container object
 */
class TrkrTruthTrackContainerv1 : public TrkrTruthTrackContainer
{
 public:
  void Reset() override;
  /* void fillEmbeddedTrkIds(PHG4TruthInfoContainer*); // fill in all the embedded track ids at once */
  void addTruthTrack (TrkrTruthTrack*) override;
  TrkrTruthTrack* getTruthTrack(unsigned int trackid) override;
  TrkrTruthTrack* getTruthTrack(unsigned int trackid, PHG4TruthInfoContainer*) override;
  ConstRange getTruthTrackRange() const override;
  bool hasTrackid(unsigned int trackid) const override;
  Map& getMap() override;

  TrkrTruthTrackContainerv1() = default;

  void identify(std::ostream& os = std::cout) const override;
  int nhw_virt()  override { return 61; };

 private:
  // the data
  Map m_data {};

  ClassDefOverride(TrkrTruthTrackContainerv1, 1)
};

#endif  // G4TRACKING_TruthTrackContainerv1.h
