#ifndef TRACKBASE_TRUTHTRACKCONTAINERV1_H
#define TRACKBASE_TRUTHTRACKCONTAINERV1_H

/**
 * @file trackbase/TrkrTruthTrackContainerv1.h
 * @author D. Stewart
 * @date September 2022
 * @brief TrkrTruthTrack container object
 */

#include "TrkrTruthTrackContainer.h"
class TrkrTruthTrack;

/**
 * @brief Cluster container object
 */
class TrkrTruthTrackContainerv1 : public TrkrTruthTrackContainer
{
 public:
  void Reset() override;
  void addTruthTrack(TrkrTruthTrack*) override;
  ConstRange getTruthTrackRange() const override;
  bool hasTrackid(unsigned int trackid) const override;
  Vector& getTruthTracks() override;
  TrkrTruthTrack* getTruthTrack(unsigned int trackid) const override;

  TrkrTruthTrackContainerv1() = default;

  void identify(std::ostream& os = std::cout) const override;

 private:
  // the data
  Vector m_data;

  ClassDefOverride(TrkrTruthTrackContainerv1, 1)
};

#endif  // TRACKBASE_TruthTrackContainerv1.h
