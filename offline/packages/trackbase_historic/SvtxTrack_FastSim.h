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
  virtual ~SvtxTrack_FastSim() = default;

  //* base class copy constructor
  SvtxTrack_FastSim( const SvtxTrack& );

  // copy content from base class
  virtual void CopyFrom( const SvtxTrack& );
  virtual void CopyFrom( SvtxTrack* source )
  { CopyFrom( *source ); }

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const;
  void Reset() { *this = SvtxTrack_FastSim(); }
  int isValid() const;
  PHObject* CloneMe() const { return new SvtxTrack_FastSim(*this); }

  //!@name accessors
  //@{

  unsigned int get_truth_track_id() const
  { return _truth_track_id; }

  unsigned int get_num_measurements() const
  { return _nmeas; }

  //@}

  //!@name modifiers
  //@{

  void set_truth_track_id(unsigned int truthTrackId)
  { _truth_track_id = truthTrackId; }

  void set_num_measurements(int nmeas)
  { _nmeas = nmeas; }

  //@}

  private:
  unsigned int _truth_track_id = UINT_MAX;
  unsigned int _nmeas = 0;

  ClassDefOverride(SvtxTrack_FastSim, 1)
};

#endif /* __SVTXTRACK_FAST_SIM_H__ */
