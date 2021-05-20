/*
 * SvtxTrack_FastSim_v2.h
 */

#ifndef TRACKBASEHISTORIC_SVTXTRACKFASTSIMV2_H
#define TRACKBASEHISTORIC_SVTXTRACKFASTSIMV2_H

#include "SvtxTrack_v2.h"

// SvtxTrack_FastSim with recording of associated G4hits
class SvtxTrack_FastSim_v2 final: public SvtxTrack_v2
{
 public:

  //* constructor
  SvtxTrack_FastSim_v2() = default;

  //* base class copy constructor
  SvtxTrack_FastSim_v2( const SvtxTrack& );

  //* destructor
  virtual ~SvtxTrack_FastSim_v2() = default;

  // copy content from base class
  virtual void CopyFrom( const SvtxTrack& );
  virtual void CopyFrom( SvtxTrack* source )
  { CopyFrom( *source ); }

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const;
  void Reset() { *this = SvtxTrack_FastSim_v2(); }
  int isValid() const;

  PHObject* CloneMe() const
  { return new SvtxTrack_FastSim_v2(*this); }

  //!@name accessors
  //@{

  unsigned int get_truth_track_id() const
  { return _truth_track_id; }

  unsigned int get_num_measurements() const
  { return _nmeas; }

  const HitIdMap& g4hit_ids() const
  { return _g4hit_ids; }

  bool empty_g4hit_id() const
  { return _g4hit_ids.empty(); }

  size_t size_g4hit_id() const
  { return _g4hit_ids.size(); }

  SvtxTrack::HitIdConstIter begin_g4hit_id() const
  { return _g4hit_ids.begin(); }

  SvtxTrack::HitIdConstIter end_g4hit_id() const
  { return _g4hit_ids.end(); }

  SvtxTrack::HitIdConstIter find_g4hit_id(int volume) const
  { return _g4hit_ids.find(volume); }

  //@}


  //!@name modifiers
  //@{

  void set_truth_track_id(unsigned int truthTrackId)
  { _truth_track_id = truthTrackId; }

  void set_num_measurements(int nmeas)
  { _nmeas = nmeas; }

  void add_g4hit_id(int volume, PHG4HitDefs::keytype id)
  { _g4hit_ids[volume].insert(id); }

  size_t remove_g4hit_id(int volume, PHG4HitDefs::keytype id)
  { return _g4hit_ids[volume].erase(id); }

  size_t remove_g4hit_volume(int volume)
  { return _g4hit_ids.erase(volume); }

  SvtxTrack::HitIdIter begin_g4hit_id()
  { return _g4hit_ids.begin(); }

  SvtxTrack::HitIdIter end_g4hit_id()
  { return _g4hit_ids.end(); }

  SvtxTrack::HitIdIter find_g4hit_id(int volume)
  { return _g4hit_ids.find(volume); }

  void clear_g4hit_id()
  { return _g4hit_ids.clear(); }

  //@}

 private:

  unsigned int _truth_track_id = UINT_MAX;
  unsigned int _nmeas = 0;

  HitIdMap _g4hit_ids;

  ClassDefOverride(SvtxTrack_FastSim_v2, 1)
};

#endif /* __SVTXTRACK_FAST_SIMV1_H__ */
