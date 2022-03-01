/*
 * SvtxTrack_FastSim_v1.h
 */

#ifndef TRACKBASEHISTORIC_SVTXTRACKFASTSIMV1_H
#define TRACKBASEHISTORIC_SVTXTRACKFASTSIMV1_H

#include "SvtxTrack_FastSim.h"

// SvtxTrack_FastSim with recording of associated G4hits
class SvtxTrack_FastSim_v1 final: public SvtxTrack_FastSim
{
 public:

  //* constructor
  SvtxTrack_FastSim_v1() = default;

  //* base class copy constructor
  SvtxTrack_FastSim_v1( const SvtxTrack& );

  //* destructor
  ~SvtxTrack_FastSim_v1() override = default;

  // copy content from base class
  using PHObject::CopyFrom; // avoid warning for not implemented CopyFrom methods
  void CopyFrom( const SvtxTrack& ) override;
  void CopyFrom( SvtxTrack* source ) override
  { CopyFrom( *source ); }

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override;
  void Reset()  override
  { *this = SvtxTrack_FastSim_v1(); }
  
  int isValid() const override;
  
  PHObject* CloneMe() const  override
  { return new SvtxTrack_FastSim_v1(*this); }


  //!@name accessors
  //@{

  const HitIdMap& g4hit_ids() const override
  { return _g4hit_ids; }

  bool empty_g4hit_id() const override
  { return _g4hit_ids.empty(); }

  size_t size_g4hit_id() const override
  { return _g4hit_ids.size(); }

  SvtxTrack::HitIdConstIter begin_g4hit_id() const override
  { return _g4hit_ids.begin(); }

  SvtxTrack::HitIdConstIter end_g4hit_id() const override
  { return _g4hit_ids.end(); }

  SvtxTrack::HitIdConstIter find_g4hit_id(int volume) const override
  { return _g4hit_ids.find(volume); }

  //@}


  //!@name modifiers
  //@{

  void add_g4hit_id(int volume, PHG4HitDefs::keytype id) override
  { _g4hit_ids[volume].insert(id); }

  size_t remove_g4hit_id(int volume, PHG4HitDefs::keytype id) override
  { return _g4hit_ids[volume].erase(id); }

  size_t remove_g4hit_volume(int volume) override
  { return _g4hit_ids.erase(volume); }

  SvtxTrack::HitIdIter begin_g4hit_id() override
  { return _g4hit_ids.begin(); }

  SvtxTrack::HitIdIter end_g4hit_id() override
  { return _g4hit_ids.end(); }

  SvtxTrack::HitIdIter find_g4hit_id(int volume) override
  { return _g4hit_ids.find(volume); }

  void clear_g4hit_id() override
  { return _g4hit_ids.clear(); }

  //@}

 private:

  HitIdMap _g4hit_ids;

  ClassDefOverride(SvtxTrack_FastSim_v1, 1)
};

#endif /* __SVTXTRACK_FAST_SIMV1_H__ */
