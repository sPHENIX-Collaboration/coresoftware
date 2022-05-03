/*
 * SvtxTrack_FastSim.h
 *
 *  Created on: Jul 28, 2016
 *      Author: yuhw
 */

#ifndef TRACKBASEHISTORIC_SVTXTRACKFASTSIM_H
#define TRACKBASEHISTORIC_SVTXTRACKFASTSIM_H

#include "SvtxTrack_v1.h"

#include <iostream>        // for cout, ostream

class SvtxTrack_FastSim: public SvtxTrack_v1
{
 public:

  //* constructor
  SvtxTrack_FastSim() = default;

  //* destructor
  ~SvtxTrack_FastSim() override = default;

  //* base class copy constructor
  SvtxTrack_FastSim( const SvtxTrack& );

  // copy content from base class
  using PHObject::CopyFrom; // avoid warning for not implemented CopyFrom methods
  void CopyFrom( const SvtxTrack& ) override;
  void CopyFrom( SvtxTrack* source ) override
  { CopyFrom( *source ); }

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = SvtxTrack_FastSim(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new SvtxTrack_FastSim(*this); }

  //!@name accessors
  //@{

  unsigned int get_truth_track_id() const override
  { return _truth_track_id; }

  unsigned int get_num_measurements() const override
  { return _nmeas; }

  //@}

  //!@name modifiers
  //@{

  void set_truth_track_id(unsigned int truthTrackId) override
  { _truth_track_id = truthTrackId; }

  void set_num_measurements(int nmeas) override
  { _nmeas = nmeas; }

  //@}

  private:
  unsigned int _truth_track_id = UINT_MAX;
  unsigned int _nmeas = 0;

  ClassDefOverride(SvtxTrack_FastSim, 1)
};

#endif /* __SVTXTRACK_FAST_SIM_H__ */
