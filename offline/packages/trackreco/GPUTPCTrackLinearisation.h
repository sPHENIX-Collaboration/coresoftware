// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCTrackLinearisation.h
/// \author Sergey Gorbunov, David Rohr

#ifndef GPUTPCTRACKLINEARISATION_H
#define GPUTPCTRACKLINEARISATION_H

#include "GPUTPCTrackParam.h"

/**
 * @class GPUTPCTrackLinearisation
 *
 * GPUTPCTrackLinearisation class describes the parameters which are used
 * to linearise the transport equations for the track trajectory.
 *
 * The class is used during track (re)fit, when the AliHLTTPCTrackParam track is only
 * partially fitted, and there is some apriory knowledge about trajectory.
 * This apriory knowledge is used to linearise the transport equations.
 *
 * In case the track is fully fitted, the best linearisation point is
 * the track trajectory itself (GPUTPCTrackLinearisation = AliHLTTPCTrackParam ).
 *
 */
class GPUTPCTrackLinearisation
{
 public:
  GPUTPCTrackLinearisation() : mSinPhi(0), mCosPhi(1), mDzDs(0), mQPt(0) {}
  GPUTPCTrackLinearisation(double SinPhi1, double CosPhi1, double DzDs1, double QPt1) : mSinPhi(SinPhi1), mCosPhi(CosPhi1), mDzDs(DzDs1), mQPt(QPt1) {}

  GPUTPCTrackLinearisation(const GPUTPCTrackParam& t);

  void Set(double SinPhi1, double CosPhi1, double DzDs1, double QPt1);

  double SinPhi() const { return mSinPhi; }
  double CosPhi() const { return mCosPhi; }
  double DzDs() const { return mDzDs; }
  double QPt() const { return mQPt; }

  double GetSinPhi() const { return mSinPhi; }
  double GetCosPhi() const { return mCosPhi; }
  double GetDzDs() const { return mDzDs; }
  double GetQPt() const { return mQPt; }

  void SetSinPhi(double v) { mSinPhi = v; }
  void SetCosPhi(double v) { mCosPhi = v; }
  void SetDzDs(double v) { mDzDs = v; }
  void SetQPt(double v) { mQPt = v; }

 private:
  double mSinPhi; // SinPhi
  double mCosPhi; // CosPhi
  double mDzDs;   // DzDs
  double mQPt;    // QPt
};

inline GPUTPCTrackLinearisation::GPUTPCTrackLinearisation(const GPUTPCTrackParam& t) : mSinPhi(t.SinPhi()), mCosPhi(0), mDzDs(t.DzDs()), mQPt(t.QPt())
{
  if (mSinPhi > GPUCA_MAX_SIN_PHI) {
    mSinPhi = GPUCA_MAX_SIN_PHI;
  } else if (mSinPhi < -GPUCA_MAX_SIN_PHI) {
    mSinPhi = -GPUCA_MAX_SIN_PHI;
  }
  mCosPhi = sqrt(1 - mSinPhi * mSinPhi);
  if (t.SignCosPhi() < 0) {
    mCosPhi = -mCosPhi;
  }
}

inline void GPUTPCTrackLinearisation::Set(double SinPhi1, double CosPhi1, double DzDs1, double QPt1)
{
  SetSinPhi(SinPhi1);
  SetCosPhi(CosPhi1);
  SetDzDs(DzDs1);
  SetQPt(QPt1);
}

#endif // GPUTPCTRACKLINEARISATION_H
