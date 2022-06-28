#ifndef TRACKBASEHISTORIC_TRACKSEED_FASTSIM_V1_H
#define TRACKBASEHISTORIC_TRACKSEED_FASTSIM_V1_H


#include "TrackSeed_v1.h"


class TrackSeed_FastSim_v1 : public TrackSeed_v1
{


 public:
  TrackSeed_FastSim_v1() = default;
  TrackSeed_FastSim_v1(const TrackSeed& );
  ~TrackSeed_FastSim_v1() = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = TrackSeed_FastSim_v1(); }
  int isValid() const override { return 1; }
  void CopyFrom(const TrackSeed& ) override;
  void CopyFrom(TrackSeed* seed ) override { CopyFrom(*seed); }

  PHObject* CloneMe() const override { return new TrackSeed_FastSim_v1(*this); }

  unsigned int get_truth_track_id() const override
  { return m_truth_track_id; }

  unsigned int get_num_measurements() const override
  { return m_nmeas; }

  const HitIdMap& g4hit_ids() const override
  { return m_g4hit_ids; }

  bool empty_g4hit_id() const override
  { return m_g4hit_ids.empty(); }

  size_t size_g4hit_id() const override
  { return m_g4hit_ids.size(); }

  TrackSeed::HitIdConstIter begin_g4hit_id() const override
  { return m_g4hit_ids.begin(); }

  TrackSeed::HitIdConstIter end_g4hit_id() const override
  { return m_g4hit_ids.end(); }

  TrackSeed::HitIdConstIter find_g4hit_id(int volume) const override
  { return m_g4hit_ids.find(volume); }

  /// We need a separate function for truth tracks because the bend
  /// angle is already properly accounted for
  float get_phi(TrkrClusterContainer *clusters,
		ActsGeometry *tGeometry) const override;

  void set_truth_track_id(unsigned int truthTrackId) override
  { m_truth_track_id = truthTrackId; }

  void set_num_measurements(int nmeas) override
  { m_nmeas = nmeas; }

  void add_g4hit_id(int volume, PHG4HitDefs::keytype id) override
  { m_g4hit_ids[volume].insert(id); }

  size_t remove_g4hit_id(int volume, PHG4HitDefs::keytype id) override
  { return m_g4hit_ids[volume].erase(id); }

  size_t remove_g4hit_volume(int volume) override
  { return m_g4hit_ids.erase(volume); }

  TrackSeed::HitIdIter begin_g4hit_id() override
  { return m_g4hit_ids.begin(); }

  TrackSeed::HitIdIter end_g4hit_id() override
  { return m_g4hit_ids.end(); }

  TrackSeed::HitIdIter find_g4hit_id(int volume) override
  { return m_g4hit_ids.find(volume); }

  void clear_g4hit_id() override
  { return m_g4hit_ids.clear(); }


 private:
  unsigned int m_truth_track_id = std::numeric_limits<unsigned int>::max();
  unsigned int m_nmeas = 0;
  HitIdMap m_g4hit_ids;

  ClassDefOverride(TrackSeed_FastSim_v1, 1)
  
};


#endif
