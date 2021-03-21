/*
 * SvtxTrack_FastSim_v1.h
 */

#ifndef TRACKBASEHISTORIC_SVTXTRACKFASTSIMV1_H
#define TRACKBASEHISTORIC_SVTXTRACKFASTSIMV1_H

#include "SvtxTrack_FastSim.h"

// SvtxTrack_FastSim with recording of associated G4hits
class SvtxTrack_FastSim_v1 : public SvtxTrack_FastSim
{
 public:
  SvtxTrack_FastSim_v1();
  virtual ~SvtxTrack_FastSim_v1();

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const;
  void Reset() { *this = SvtxTrack_FastSim_v1(); }
  int isValid() const;
  PHObject* CloneMe() const { return new SvtxTrack_FastSim_v1(*this); }

  bool empty_g4hit_id() const { return _g4hit_ids.empty(); }
  size_t size_g4hit_id() const { return _g4hit_ids.size(); }
  void add_g4hit_id(int volume, PHG4HitDefs::keytype id) { _g4hit_ids[volume].insert(id); }
  SvtxTrack::HitIdIter begin_g4hit_id() { return _g4hit_ids.begin(); }
  SvtxTrack::HitIdConstIter begin_g4hit_id() const { return _g4hit_ids.begin(); }
  SvtxTrack::HitIdIter find_g4hit_id(int volume) { return _g4hit_ids.find(volume); }
  SvtxTrack::HitIdConstIter find_g4hit_id(int volume) const { return _g4hit_ids.find(volume); }
  SvtxTrack::HitIdIter end_g4hit_id() { return _g4hit_ids.end(); }
  SvtxTrack::HitIdConstIter end_g4hit_id() const { return _g4hit_ids.end(); }
  size_t remove_g4hit_id(int volume, PHG4HitDefs::keytype id) { return _g4hit_ids[volume].erase(id); }
  size_t remove_g4hit_volume(int volume) { return _g4hit_ids.erase(volume); }
  void clear_g4hit_id() { return _g4hit_ids.clear(); }
  const HitIdMap& g4hit_ids() const { return _g4hit_ids; }

 private:
  HitIdMap _g4hit_ids;  //< contained hit ids

  ClassDef(SvtxTrack_FastSim_v1, 1)
};

#endif /* __SVTXTRACK_FAST_SIMV1_H__ */
