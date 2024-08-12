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
    double bethe, e, theta2, EP2, sigmadE2, k22, k33, k43, k44; // parameters
  };

   const GPUTPCBaseTrackParam& GetParam() const { return mParam; }
   void SetParam(const GPUTPCBaseTrackParam& v) { mParam = v; }
   void InitParam();

   double X() const { return mParam.X(); }
   double Y() const { return mParam.Y(); }
   double Z() const { return mParam.Z(); }
   double SinPhi() const { return mParam.SinPhi(); }
   double DzDs() const { return mParam.DzDs(); }
   double QPt() const { return mParam.QPt(); }
   double ZOffset() const { return mParam.ZOffset(); }
   double SignCosPhi() const { return mSignCosPhi; }
   double Chi2() const { return mChi2; }
   int NDF() const { return mNDF; }

   double Err2Y() const { return mC[0]; }
   double Err2Z() const { return mC[2]; }
   double Err2SinPhi() const { return mC[5]; }
   double Err2DzDs() const { return mC[9]; }
   double Err2QPt() const { return mC[14]; }

   double GetX() const { return mParam.GetX(); }
   double GetY() const { return mParam.GetY(); }
   double GetZ() const { return mParam.GetZ(); }
   double GetSinPhi() const { return mParam.GetSinPhi(); }
   double GetDzDs() const { return mParam.GetDzDs(); }
   double GetQPt() const { return mParam.GetQPt(); }
   double GetSignCosPhi() const { return mSignCosPhi; }
   double GetChi2() const { return mChi2; }
   int GetNDF() const { return mNDF; }

   double GetKappa(double Bz) const { return mParam.GetKappa(Bz); }
   double GetCosPhi() const { return mSignCosPhi * sqrt(1 - SinPhi() * SinPhi()); }

   double GetErr2Y() const { return mC[0]; }
   double GetErr2Z() const { return mC[2]; }
   double GetErr2SinPhi() const { return mC[5]; }
   double GetErr2DzDs() const { return mC[9]; }
   double GetErr2QPt() const { return mC[14]; }

   const double* Par() const { return mParam.Par(); }
   const double* Cov() const { return mC; }

   const double* GetPar() const { return mParam.GetPar(); }
   double GetPar(int i) const { return (mParam.GetPar(i)); }
   const double* GetCov() const { return mC; }
   double GetCov(int i) const { return mC[i]; }

   void SetPar(int i, double v) { mParam.SetPar(i, v); }
   void SetCov(int i, double v) { mC[i] = v; }

   void SetX(double v) { mParam.SetX(v); }
   void SetY(double v) { mParam.SetY(v); }
   void SetZ(double v) { mParam.SetZ(v); }
   void SetSinPhi(double v) { mParam.SetSinPhi(v); }
   void SetDzDs(double v) { mParam.SetDzDs(v); }
   void SetQPt(double v) { mParam.SetQPt(v); }
   void SetZOffset(double v) { mParam.SetZOffset(v); }
   void SetSignCosPhi(double v) { mSignCosPhi = v >= 0 ? 1 : -1; }
   void SetChi2(double v) { mChi2 = v; }
   void SetNDF(int v) { mNDF = v; }

   double GetDist2(const GPUTPCTrackParam& t) const;
   double GetDistXZ2(const GPUTPCTrackParam& t) const;

   double GetS(double x, double y, double Bz) const;

   void GetDCAPoint(double x, double y, double z, double& px, double& py, double& pz, double Bz) const;

   bool TransportToX(double x, double Bz, double maxSinPhi = GPUCA_MAX_SIN_PHI);
   bool TransportToXWithMaterial(double x, double Bz, double maxSinPhi = GPUCA_MAX_SIN_PHI);

   bool TransportToX(double x, GPUTPCTrackLinearisation& t0, double Bz, double maxSinPhi = GPUCA_MAX_SIN_PHI, double* DL = nullptr);

   bool TransportToX(double x, double sinPhi0, double cosPhi0, double Bz, double maxSinPhi = GPUCA_MAX_SIN_PHI);

   bool TransportToXWithMaterial(double x, GPUTPCTrackLinearisation& t0, GPUTPCTrackFitParam& par, double Bz, double maxSinPhi = GPUCA_MAX_SIN_PHI);

   bool TransportToXWithMaterial(double x, GPUTPCTrackFitParam& par, double Bz, double maxSinPhi = GPUCA_MAX_SIN_PHI);

   static double ApproximateBetheBloch(double beta2);
   static double BetheBlochGeant(double bg, double kp0 = 2.33f, double kp1 = 0.20f, double kp2 = 3.00f, double kp3 = 173e-9f, double kp4 = 0.49848f);
   static double BetheBlochSolid(double bg);
   double BetheBlochGas(double bg);

   void CalculateFitParameters(GPUTPCTrackFitParam& par, double mass = 0.13957f);
   bool CorrectForMeanMaterial(double xOverX0, double xTimesRho, const GPUTPCTrackFitParam& par);

   bool Rotate(double alpha, double maxSinPhi = GPUCA_MAX_SIN_PHI);
   bool Rotate(double alpha, GPUTPCTrackLinearisation& t0, double maxSinPhi = GPUCA_MAX_SIN_PHI);
   bool Filter(double y, double z, double err2Y, double err2Z, double maxSinPhi = GPUCA_MAX_SIN_PHI, bool paramOnly = false);

   bool CheckNumericalQuality() const;

   void Print() const;

   void setNeonFraction(double frac) { Ne_frac = frac; };
   void setArgonFraction(double frac) { Ar_frac = frac; };
   void setCF4Fraction(double frac) { CF4_frac = frac; };
   void setNitrogenFraction(double frac) { N2_frac = frac; };
   void setIsobutaneFraction(double frac) { isobutane_frac = frac; };

#ifndef GPUCA_GPUCODE
 private:
#endif //! GPUCA_GPUCODE
  GPUTPCBaseTrackParam
  mParam; // Track Parameters

 private:
  // WARNING, Track Param Data is copied in the GPU Tracklet Constructor element by element instead of using copy constructor!!!
  // This is neccessary for performance reasons!!!
  // Changes to Elements of this class therefore must also be applied to TrackletConstructor!!!
  double mC[15];      // the covariance matrix for Y,Z,SinPhi,..
  double mSignCosPhi; // sign of cosPhi
  double mChi2;       // the chi^2 value
  int mNDF;          // the Number of Degrees of Freedom


  //Gas parameters
  //Sources
  //https://doi.org/10.1016/0092-640X(84)90002-0
  //https://pdg.lbl.gov/2024/AtomicNuclearProperties/index.html
  //https://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf

  double Ne_frac = 0.00;
  double Ne_Rho = 0.839e-3; // g/cm3, density
  double Ne_RadLen = 28.93; // g/cm2, radiation length
  double Ne_x0 = 2.0735; // 1st Bethe-Bloch pole
  double Ne_x1 = 4.6421; //2nd Bethe-Bloch pole
  double Ne_mI = 137; //eV, mean excitation
  double Ne_mZA = 0.4955; // Mean atomic number / mass number
  double Ne_mA = 20.18; // Mean amass number

  double Ar_frac = 0.75;
  double Ar_Rho = 1.662e-3; // g/cm3, density
  double Ar_RadLen = 19.55; // g/cm2, radiation length
  double Ar_x0 = 1.7635; // 1st Bethe-Bloch pole
  double Ar_x1 = 4.4855; //2nd Bethe-Bloch pole
  double Ar_mI = 188; //eV, mean excitation
  double Ar_mZA = 0.4506; // Mean atomic number / mass number
  double Ar_mA = 39.948; // Mean amass number

  double CF4_frac = 0.20;
  double CF4_Rho = 3.78e-3; // g/cm3, density
  double CF4_RadLen = 33.99; // g/cm2, radiation length
  double CF4_x0 = 1.7635; // 1st Bethe-Bloch pole
  double CF4_x1 = 4.4855; //2nd Bethe-Bloch pole
  double CF4_mI = 115; //eV, mean excitation
  double CF4_mZA = 0.4772; // Mean atomic number / mass number
  double CF4_mA = 88.013; // Mean amass number

  double N2_frac = 0.00;
  double N2_Rho = 1.165e-3; // g/cm3, density
  double N2_RadLen = 37.99; // g/cm2, radiation length
  double N2_x0 = 1.7378; // 1st Bethe-Bloch pole
  double N2_x1 = 4.1323; //2nd Bethe-Bloch pole
  double N2_mI = 82; //eV, mean excitation
  double N2_mZA = 0.4998; // Mean atomic number / mass number
  double N2_mA = 14.007; // Mean amass number

  double isobutane_frac = 0.05;
  double isobutane_Rho = 2.59e-3; // g/cm3, density
  double isobutane_RadLen = 45.23; // g/cm2, radiation length
  double isobutane_x0 = 1.3788; // 1st Bethe-Bloch pole
  double isobutane_x1 = 3.7524; //2nd Bethe-Bloch pole
  double isobutane_mI = 48.3; //eV, mean excitation
  double isobutane_mZA = 0.595; // Mean atomic number / mass number
  double isobutane_mA = 57.146; // Mean amass number
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
