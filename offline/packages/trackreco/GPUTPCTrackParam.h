// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCTrackParam.h
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#ifndef GPUTPCTRACKPARAM_H
#define GPUTPCTRACKPARAM_H

#include "GPUTPCBaseTrackParam.h"
#include <cmath>

class GPUTPCTrackLinearisation;
/**
 * @class GPUTPCTrackParam
 *
 * GPUTPCTrackParam class describes the track parametrisation
 * which is used by the GPUTPCTracker slice tracker.
 *
 */
class GPUTPCTrackParam
{
 public:
  struct GPUTPCTrackFitParam {
    float bethe, e, theta2, EP2, sigmadE2, k22, k33, k43, k44; // parameters
  };

   const GPUTPCBaseTrackParam& GetParam() const { return mParam; }
   void SetParam(const GPUTPCBaseTrackParam& v) { mParam = v; }
   void InitParam();

   float X() const { return mParam.X(); }
   float Y() const { return mParam.Y(); }
   float Z() const { return mParam.Z(); }
   float SinPhi() const { return mParam.SinPhi(); }
   float DzDs() const { return mParam.DzDs(); }
   float QPt() const { return mParam.QPt(); }
   float ZOffset() const { return mParam.ZOffset(); }
   float SignCosPhi() const { return mSignCosPhi; }
   float Chi2() const { return mChi2; }
   int NDF() const { return mNDF; }

   float Err2Y() const { return mC[0]; }
   float Err2Z() const { return mC[2]; }
   float Err2SinPhi() const { return mC[5]; }
   float Err2DzDs() const { return mC[9]; }
   float Err2QPt() const { return mC[14]; }

   float GetX() const { return mParam.GetX(); }
   float GetY() const { return mParam.GetY(); }
   float GetZ() const { return mParam.GetZ(); }
   float GetSinPhi() const { return mParam.GetSinPhi(); }
   float GetDzDs() const { return mParam.GetDzDs(); }
   float GetQPt() const { return mParam.GetQPt(); }
   float GetSignCosPhi() const { return mSignCosPhi; }
   float GetChi2() const { return mChi2; }
   int GetNDF() const { return mNDF; }

   float GetKappa(float Bz) const { return mParam.GetKappa(Bz); }
   float GetCosPhi() const { return mSignCosPhi * sqrt(1 - SinPhi() * SinPhi()); }

   float GetErr2Y() const { return mC[0]; }
   float GetErr2Z() const { return mC[2]; }
   float GetErr2SinPhi() const { return mC[5]; }
   float GetErr2DzDs() const { return mC[9]; }
   float GetErr2QPt() const { return mC[14]; }

   const float* Par() const { return mParam.Par(); }
   const float* Cov() const { return mC; }

   const float* GetPar() const { return mParam.GetPar(); }
   float GetPar(int i) const { return (mParam.GetPar(i)); }
   const float* GetCov() const { return mC; }
   float GetCov(int i) const { return mC[i]; }

   void SetPar(int i, float v) { mParam.SetPar(i, v); }
   void SetCov(int i, float v) { mC[i] = v; }

   void SetX(float v) { mParam.SetX(v); }
   void SetY(float v) { mParam.SetY(v); }
   void SetZ(float v) { mParam.SetZ(v); }
   void SetSinPhi(float v) { mParam.SetSinPhi(v); }
   void SetDzDs(float v) { mParam.SetDzDs(v); }
   void SetQPt(float v) { mParam.SetQPt(v); }
   void SetZOffset(float v) { mParam.SetZOffset(v); }
   void SetSignCosPhi(float v) { mSignCosPhi = v >= 0 ? 1 : -1; }
   void SetChi2(float v) { mChi2 = v; }
   void SetNDF(int v) { mNDF = v; }

   float GetDist2(const GPUTPCTrackParam& t) const;
   float GetDistXZ2(const GPUTPCTrackParam& t) const;

   float GetS(float x, float y, float Bz) const;

   void GetDCAPoint(float x, float y, float z, float& px, float& py, float& pz, float Bz) const;

   bool TransportToX(float x, float Bz, float maxSinPhi = GPUCA_MAX_SIN_PHI);
   bool TransportToXWithMaterial(float x, float Bz, float maxSinPhi = GPUCA_MAX_SIN_PHI);

   bool TransportToX(float x, GPUTPCTrackLinearisation& t0, float Bz, float maxSinPhi = GPUCA_MAX_SIN_PHI, float* DL = nullptr);

   bool TransportToX(float x, float sinPhi0, float cosPhi0, float Bz, float maxSinPhi = GPUCA_MAX_SIN_PHI);

   bool TransportToXWithMaterial(float x, GPUTPCTrackLinearisation& t0, GPUTPCTrackFitParam& par, float Bz, float maxSinPhi = GPUCA_MAX_SIN_PHI);

   bool TransportToXWithMaterial(float x, GPUTPCTrackFitParam& par, float Bz, float maxSinPhi = GPUCA_MAX_SIN_PHI);

   static float ApproximateBetheBloch(float beta2);
   static float BetheBlochGeant(float bg, float kp0 = 2.33f, float kp1 = 0.20f, float kp2 = 3.00f, float kp3 = 173e-9f, float kp4 = 0.49848f);
   static float BetheBlochSolid(float bg);
   static float BetheBlochGas(float bg);

   void CalculateFitParameters(GPUTPCTrackFitParam& par, float mass = 0.13957f);
   bool CorrectForMeanMaterial(float xOverX0, float xTimesRho, const GPUTPCTrackFitParam& par);

   bool Rotate(float alpha, float maxSinPhi = GPUCA_MAX_SIN_PHI);
   bool Rotate(float alpha, GPUTPCTrackLinearisation& t0, float maxSinPhi = GPUCA_MAX_SIN_PHI);
   bool Filter(float y, float z, float err2Y, float err2Z, float maxSinPhi = GPUCA_MAX_SIN_PHI, bool paramOnly = false);

   bool CheckNumericalQuality() const;

   void Print() const;

#ifndef GPUCA_GPUCODE
 private:
#endif //! GPUCA_GPUCODE
  GPUTPCBaseTrackParam
  mParam; // Track Parameters

 private:
  // WARNING, Track Param Data is copied in the GPU Tracklet Constructor element by element instead of using copy constructor!!!
  // This is neccessary for performance reasons!!!
  // Changes to Elements of this class therefore must also be applied to TrackletConstructor!!!
  float mC[15];      // the covariance matrix for Y,Z,SinPhi,..
  float mSignCosPhi; // sign of cosPhi
  float mChi2;       // the chi^2 value
  int mNDF;          // the Number of Degrees of Freedom
};

inline void GPUTPCTrackParam::InitParam()
{
  // Initialize Tracklet Parameters using default values
  SetSinPhi(0);
  SetDzDs(0);
  SetQPt(0);
  SetSignCosPhi(1);
  SetChi2(0);
  SetNDF(-3);
  SetCov(0, 1);
  SetCov(1, 0);
  SetCov(2, 1);
  SetCov(3, 0);
  SetCov(4, 0);
  SetCov(5, 1);
  SetCov(6, 0);
  SetCov(7, 0);
  SetCov(8, 0);
  SetCov(9, 1);
  SetCov(10, 0);
  SetCov(11, 0);
  SetCov(12, 0);
  SetCov(13, 0);
  SetCov(14, 1000.f);
  SetZOffset(0);
}

#endif // GPUTPCTRACKPARAM_H
