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

constexpr double GPUCA_MAX_SIN_PHI = 0.99;

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
  double X() const { return mX; }
  double Y() const { return mP[0]; }
  double Z() const { return mP[1]; }
  double SinPhi() const { return mP[2]; }
  double DzDs() const { return mP[3]; }
  double QPt() const { return mP[4]; }
  double ZOffset() const { return mZOffset; }

  double GetX() const { return mX; }
  double GetY() const { return mP[0]; }
  double GetZ() const { return mP[1]; }
  double GetSinPhi() const { return mP[2]; }
  double GetDzDs() const { return mP[3]; }
  double GetQPt() const { return mP[4]; }
  double GetZOffset() const { return mZOffset; }

  double GetKappa(double Bz) const { return -mP[4] * Bz; }

  const double* Par() const { return mP; }
  const double* GetPar() const { return mP; }
  double GetPar(int i) const { return (mP[i]); }

  void SetPar(int i, double v) { mP[i] = v; }

  void SetX(double v) { mX = v; }
  void SetY(double v) { mP[0] = v; }
  void SetZ(double v) { mP[1] = v; }
  void SetSinPhi(double v) { mP[2] = v; }
  void SetDzDs(double v) { mP[3] = v; }
  void SetQPt(double v) { mP[4] = v; }
  void SetZOffset(double v) { mZOffset = v; }

 private:
  // WARNING, Track Param Data is copied in the GPU Tracklet Constructor element by element instead of using copy constructor!!!
  // This is neccessary for performance reasons!!!
  // Changes to Elements of this class therefore must also be applied to TrackletConstructor!!!
  double mX; // x position
  double mZOffset;
  double mP[5]; // 'active' track parameters: Y, Z, SinPhi, DzDs, q/Pt
};

#endif
