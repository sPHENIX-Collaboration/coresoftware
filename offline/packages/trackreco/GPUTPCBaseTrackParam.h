// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCBaseTrackParam.h
/// \author David Rohr, Sergey Gorbunov

#ifndef GPUTPCBASETRACKPARAM_H
#define GPUTPCBASETRACKPARAM_H

constexpr float GPUCA_MAX_SIN_PHI = 0.99;

class GPUTPCTrackParam;

/**
 * @class GPUTPCBaseTrackParam
 *
 * GPUTPCBaseTrackParam class contains track parameters
 * used in output of the GPUTPCTracker slice tracker.
 * This class is used for transfer between tracker and merger and does not contain the covariance matrice
 */
class GPUTPCBaseTrackParam
{
 public:
  float X() const { return mX; }
  float Y() const { return mP[0]; }
  float Z() const { return mP[1]; }
  float SinPhi() const { return mP[2]; }
  float DzDs() const { return mP[3]; }
  float QPt() const { return mP[4]; }
  float ZOffset() const { return mZOffset; }

  float GetX() const { return mX; }
  float GetY() const { return mP[0]; }
  float GetZ() const { return mP[1]; }
  float GetSinPhi() const { return mP[2]; }
  float GetDzDs() const { return mP[3]; }
  float GetQPt() const { return mP[4]; }
  float GetZOffset() const { return mZOffset; }

  float GetKappa(float Bz) const { return -mP[4] * Bz; }

  const float* Par() const { return mP; }
  const float* GetPar() const { return mP; }
  float GetPar(int i) const { return (mP[i]); }

  void SetPar(int i, float v) { mP[i] = v; }

  void SetX(float v) { mX = v; }
  void SetY(float v) { mP[0] = v; }
  void SetZ(float v) { mP[1] = v; }
  void SetSinPhi(float v) { mP[2] = v; }
  void SetDzDs(float v) { mP[3] = v; }
  void SetQPt(float v) { mP[4] = v; }
  void SetZOffset(float v) { mZOffset = v; }

 private:
  // WARNING, Track Param Data is copied in the GPU Tracklet Constructor element by element instead of using copy constructor!!!
  // This is neccessary for performance reasons!!!
  // Changes to Elements of this class therefore must also be applied to TrackletConstructor!!!
  float mX; // x position
  float mZOffset;
  float mP[5]; // 'active' track parameters: Y, Z, SinPhi, DzDs, q/Pt
};

#endif
