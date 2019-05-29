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

class SvtxTrack_FastSim : public SvtxTrack_v1
{
 public:
  SvtxTrack_FastSim();
  virtual ~SvtxTrack_FastSim();

  unsigned int get_truth_track_id() const
  {
    return _truth_track_id;
  }

  void set_truth_track_id(unsigned int truthTrackId)
  {
    _truth_track_id = truthTrackId;
  }

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const;
  void Reset() { *this = SvtxTrack_FastSim(); }
  int isValid() const;
  SvtxTrack* Clone() const { return new SvtxTrack_FastSim(*this); }

  void set_num_measurements(int nmeas) { _nmeas = nmeas; }
  unsigned int get_num_measurements() { return _nmeas; }

 private:
  unsigned int _truth_track_id;
  unsigned int _nmeas;

  ClassDef(SvtxTrack_FastSim, 1)
};

#endif /* __SVTXTRACK_FAST_SIM_H__ */
