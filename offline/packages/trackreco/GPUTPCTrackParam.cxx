// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCTrackParam.cxx
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#include "GPUTPCTrackLinearisation.h"
#include "GPUTPCTrackParam.h"
#include <cmath>
#include <algorithm>

//
// Circle in XY:
//
// kCLight = 0.000299792458;
// Kappa = -Bz*kCLight*QPt;
// R  = 1/fstd::absf(Kappa);
// Xc = X - sin(Phi)/Kappa;
// Yc = Y + cos(Phi)/Kappa;
//

double GPUTPCTrackParam::GetDist2(const GPUTPCTrackParam& t) const
{
  // get squared distance between tracks

  double dx = GetX() - t.GetX();
  double dy = GetY() - t.GetY();
  double dz = GetZ() - t.GetZ();
  return dx * dx + dy * dy + dz * dz;
}

double GPUTPCTrackParam::GetDistXZ2(const GPUTPCTrackParam& t) const
{
  // get squared distance between tracks in X&Z

  double dx = GetX() - t.GetX();
  double dz = GetZ() - t.GetZ();
  return dx * dx + dz * dz;
}

double GPUTPCTrackParam::GetS(double x, double y, double Bz) const
{
  //* Get XY path length to the given point

  double k = GetKappa(Bz);
  double ex = GetCosPhi();
  double ey = GetSinPhi();
  x -= GetX();
  y -= GetY();
  double dS = x * ex + y * ey;
  if (std::abs(k) > 1.e-4f) {
    dS = atan2(k * dS, 1 + k * (x * ey - y * ex)) / k;
  }
  return dS;
}

void GPUTPCTrackParam::GetDCAPoint(double x, double y, double z, double& xp, double& yp, double& zp, double Bz) const
{
  //* Get the track point closest to the (x,y,z)

  double x0 = GetX();
  double y0 = GetY();
  double k = GetKappa(Bz);
  double ex = GetCosPhi();
  double ey = GetSinPhi();
  double dx = x - x0;
  double dy = y - y0;
  double ax = dx * k + ey;
  double ay = dy * k - ex;
  double a = std::sqrt(ax * ax + ay * ay);
  xp = x0 + (dx - ey * ((dx * dx + dy * dy) * k - 2 * (-dx * ey + dy * ex)) / (a + 1)) / a;
  yp = y0 + (dy + ex * ((dx * dx + dy * dy) * k - 2 * (-dx * ey + dy * ex)) / (a + 1)) / a;
  double s = GetS(x, y, Bz);
  zp = GetZ() + GetDzDs() * s;
  if (std::abs(k) > 1.e-2f) {
    double dZ = std::abs(GetDzDs() * 2*M_PI / k);
    if (dZ > .1f) {
      zp += round((z - zp) / dZ) * dZ;
    }
  }
}

//*
//* Transport routines
//*

bool GPUTPCTrackParam::TransportToX(double x, GPUTPCTrackLinearisation& t0, double Bz, double maxSinPhi, double* DL)
{
  //* Transport the track parameters to X=x, using linearization at t0, and the field value Bz
  //* maxSinPhi is the max. allowed value for |t0.SinPhi()|
  //* linearisation of trajectory t0 is also transported to X=x,
  //* returns 1 if OK
  //*

  double ex = t0.CosPhi();
  double ey = t0.SinPhi();
  double k = -t0.QPt() * Bz;
  double dx = x - X();

  double ey1 = k * dx + ey;
  double ex1;

  // check for intersection with X=x

  if (std::abs(ey1) > maxSinPhi) {
    return 0;
  }

  ex1 = std::sqrt(1 - ey1 * ey1);
  if (ex < 0) {
    ex1 = -ex1;
  }

  double dx2 = dx * dx;
  double ss = ey + ey1;
  double cc = ex + ex1;

  if (std::abs(cc) < 1.e-4f || std::abs(ex) < 1.e-4f || std::abs(ex1) < 1.e-4f) {
    return 0;
  }

  double tg = ss / cc; // tanf((phi1+phi)/2)

  double dy = dx * tg;
  double dl = dx * std::sqrt(1 + tg * tg);

  if (cc < 0) {
    dl = -dl;
  }
  double dSin = dl * k / 2;
  if (dSin > 1) {
    dSin = 1;
  }
  if (dSin < -1) {
    dSin = -1;
  }
  double dS = (std::abs(k) > 1.e-4f) ? (2 * asin(dSin) / k) : dl;
  double dz = dS * t0.DzDs();

  if (DL) {
    *DL = -dS * std::sqrt(1 + t0.DzDs() * t0.DzDs());
  }

  double cci = 1.f / cc;
  double exi = 1.f / ex;
  double ex1i = 1.f / ex1;

  double d[5] = {0, 0, GetPar(2) - t0.SinPhi(), GetPar(3) - t0.DzDs(), GetPar(4) - t0.QPt()};

  // double H0[5] = { 1,0, h2,  0, h4 };
  // double H1[5] = { 0, 1, 0, dS,  0 };
  // double H2[5] = { 0, 0, 1,  0, dxBz };
  // double H3[5] = { 0, 0, 0,  1,  0 };
  // double H4[5] = { 0, 0, 0,  0,  1 };

  double h2 = dx * (1 + ey * ey1 + ex * ex1) * exi * ex1i * cci;
  double h4 = dx2 * (cc + ss * ey1 * ex1i) * cci * cci * (-Bz);
  double dxBz = dx * (-Bz);

  t0.SetCosPhi(ex1);
  t0.SetSinPhi(ey1);

  SetX(X() + dx);
  SetPar(0, Y() + dy + h2 * d[2] + h4 * d[4]);
  SetPar(1, Z() + dz + dS * d[3]);
  SetPar(2, t0.SinPhi() + d[2] + dxBz * d[4]);

  double c00 = mC[0];
  double c10 = mC[1];
  double c11 = mC[2];
  double c20 = mC[3];
  double c21 = mC[4];
  double c22 = mC[5];
  double c30 = mC[6];
  double c31 = mC[7];
  double c32 = mC[8];
  double c33 = mC[9];
  double c40 = mC[10];
  double c41 = mC[11];
  double c42 = mC[12];
  double c43 = mC[13];
  double c44 = mC[14];

  mC[0] = c00 + h2 * h2 * c22 + h4 * h4 * c44 + 2 * (h2 * c20 + h4 * c40 + h2 * h4 * c42);

  mC[1] = c10 + h2 * c21 + h4 * c41 + dS * (c30 + h2 * c32 + h4 * c43);
  mC[2] = c11 + 2 * dS * c31 + dS * dS * c33;

  mC[3] = c20 + h2 * c22 + h4 * c42 + dxBz * (c40 + h2 * c42 + h4 * c44);
  mC[4] = c21 + dS * c32 + dxBz * (c41 + dS * c43);
  mC[5] = c22 + 2 * dxBz * c42 + dxBz * dxBz * c44;

  mC[6] = c30 + h2 * c32 + h4 * c43;
  mC[7] = c31 + dS * c33;
  mC[8] = c32 + dxBz * c43;
  mC[9] = c33;

  mC[10] = c40 + h2 * c42 + h4 * c44;
  mC[11] = c41 + dS * c43;
  mC[12] = c42 + dxBz * c44;
  mC[13] = c43;
  mC[14] = c44;

  return 1;
}

bool GPUTPCTrackParam::TransportToX(double x, double sinPhi0, double cosPhi0, double Bz, double maxSinPhi)
{
  //* Transport the track parameters to X=x, using linearization at phi0 with 0 curvature,
  //* and the field value Bz
  //* maxSinPhi is the max. allowed value for |t0.SinPhi()|
  //* linearisation of trajectory t0 is also transported to X=x,
  //* returns 1 if OK
  //*

  double ex = cosPhi0;
  double ey = sinPhi0;
  double dx = x - X();

  if (std::abs(ex) < 1.e-4f) {
    return 0;
  }
  double exi = 1.f / ex;

  double dxBz = dx * (-Bz);
  double dS = dx * exi;
  double h2 = dS * exi * exi;
  double h4 = .5f * h2 * dxBz;

  // double H0[5] = { 1,0, h2,  0, h4 };
  // double H1[5] = { 0, 1, 0, dS,  0 };
  // double H2[5] = { 0, 0, 1,  0, dxBz };
  // double H3[5] = { 0, 0, 0,  1,  0 };
  // double H4[5] = { 0, 0, 0,  0,  1 };

  double sinPhi = SinPhi() + dxBz * QPt();
  if (maxSinPhi > 0 && std::abs(sinPhi) > maxSinPhi) {
    return 0;
  }

  SetX(X() + dx);
  SetPar(0, GetPar(0) + dS * ey + h2 * (SinPhi() - ey) + h4 * QPt());
  SetPar(1, GetPar(1) + dS * DzDs());
  SetPar(2, sinPhi);

  double c00 = mC[0];
  double c10 = mC[1];
  double c11 = mC[2];
  double c20 = mC[3];
  double c21 = mC[4];
  double c22 = mC[5];
  double c30 = mC[6];
  double c31 = mC[7];
  double c32 = mC[8];
  double c33 = mC[9];
  double c40 = mC[10];
  double c41 = mC[11];
  double c42 = mC[12];
  double c43 = mC[13];
  double c44 = mC[14];

  // This is not the correct propagation!!! The max ensures the positional error does not decrease during track finding.
  // A different propagation method is used for the fit.
  mC[0] = std::max(c00, c00 + h2 * h2 * c22 + h4 * h4 * c44 + 2 * (h2 * c20 + h4 * c40 + h2 * h4 * c42));

  mC[1] = c10 + h2 * c21 + h4 * c41 + dS * (c30 + h2 * c32 + h4 * c43);
  mC[2] = std::max(c11, c11 + 2 * dS * c31 + dS * dS * c33);

  mC[3] = c20 + h2 * c22 + h4 * c42 + dxBz * (c40 + h2 * c42 + h4 * c44);
  mC[4] = c21 + dS * c32 + dxBz * (c41 + dS * c43);
  mC[5] = c22 + 2 * dxBz * c42 + dxBz * dxBz * c44;

  mC[6] = c30 + h2 * c32 + h4 * c43;
  mC[7] = c31 + dS * c33;
  mC[8] = c32 + dxBz * c43;
  mC[9] = c33;

  mC[10] = c40 + h2 * c42 + h4 * c44;
  mC[11] = c41 + dS * c43;
  mC[12] = c42 + dxBz * c44;
  mC[13] = c43;
  mC[14] = c44;

  return 1;
}

bool GPUTPCTrackParam::TransportToX(double x, double Bz, double maxSinPhi)
{
  //* Transport the track parameters to X=x
  GPUTPCTrackLinearisation t0(*this);
  return TransportToX(x, t0, Bz, maxSinPhi);
}

bool GPUTPCTrackParam::TransportToXWithMaterial(double x, GPUTPCTrackLinearisation& t0, GPUTPCTrackFitParam& par, double Bz, double maxSinPhi)
{
  //* Transport the track parameters to X=x  taking into account material budget

//  const double kRho = 1.025e-3f;  // 0.9e-3;
//  const double kRadLen = 29.532f; // 28.94;
  const double kRho = 2.27925e-3f;
  const double kRadLen = 14.403f;
  const double kRhoOverRadLen = kRho / kRadLen;
  double dl;

  if (!TransportToX(x, t0, Bz, maxSinPhi, &dl)) {
    return 0;
  }

  CorrectForMeanMaterial(dl * kRhoOverRadLen, dl * kRho, par);
  return 1;
}

bool GPUTPCTrackParam::TransportToXWithMaterial(double x, GPUTPCTrackFitParam& par, double Bz, double maxSinPhi)
{
  //* Transport the track parameters to X=x  taking into account material budget

  GPUTPCTrackLinearisation t0(*this);
  return TransportToXWithMaterial(x, t0, par, Bz, maxSinPhi);
}

bool GPUTPCTrackParam::TransportToXWithMaterial(double x, double Bz, double maxSinPhi)
{
  //* Transport the track parameters to X=x taking into account material budget

  GPUTPCTrackFitParam par;
  CalculateFitParameters(par);
  return TransportToXWithMaterial(x, par, Bz, maxSinPhi);
}

//*
//*  Multiple scattering and energy losses
//*
double GPUTPCTrackParam::BetheBlochGeant(double bg2, double kp0, double kp1, double kp2, double kp3, double kp4)
{
  //
  // This is the parameterization of the Bethe-Bloch formula inspired by Geant.
  //
  // bg2  - (beta*gamma)^2
  // kp0 - density [g/cm^3]
  // kp1 - density effect first junction point
  // kp2 - density effect second junction point
  // kp3 - mean excitation energy [GeV]
  // kp4 - mean Z/A
  //
  // The default values for the kp* parameters are for silicon.
  // The returned value is in [GeV/(g/cm^2)].
  //

  const double mK = 0.307075e-3f; // [GeV*cm^2/g]
  const double me = 0.511e-3f;    // [GeV/c^2]
  const double rho = kp0;
  const double x0 = kp1 * 2.303f;
  const double x1 = kp2 * 2.303f;
  const double mI = kp3;
  const double mZA = kp4;
  const double maxT = 2 * me * bg2; // neglecting the electron mass

  //*** Density effect
  double d2 = 0.f;
  const double x = 0.5f * log(bg2);
  const double lhwI = log(28.816f * 1e-9f * std::sqrt(rho * mZA) / mI);
  if (x > x1) {
    d2 = lhwI + x - 0.5f;
  } else if (x > x0) {
    const double r = (x1 - x) / (x1 - x0);
    d2 = lhwI + x - 0.5f + (0.5f - lhwI - x0) * r * r * r;
  }

  return mK * mZA * (1 + bg2) / bg2 * (0.5f * log(2 * me * bg2 * maxT / (mI * mI)) - bg2 / (1 + bg2) - d2);
}

double GPUTPCTrackParam::BetheBlochSolid(double bg)
{
  //------------------------------------------------------------------
  // This is an approximation of the Bethe-Bloch formula,
  // reasonable for solid materials.
  // All the parameters are, in fact, for Si.
  // The returned value is in [GeV]
  //------------------------------------------------------------------

  return BetheBlochGeant(bg);
}

double GPUTPCTrackParam::BetheBlochGas(double bg)
{
  //------------------------------------------------------------------
  // This is an approximation of the Bethe-Bloch formula,
  // reasonable for gas materials.
  // All the parameters are, in fact, for Ne.
  // The returned value is in [GeV]
  //------------------------------------------------------------------

//  const double rho = 0.9e-3f;
//  const double x0 = 2.f;
//  const double x1 = 4.f;
//  const double mI = 140.e-9f;
//  const double mZA = 0.49555f;

  // these are for the sPHENIX TPC gas mixture
  const double rho = 2.27925e-3f;
  const double x0 = 2.f;
  const double x1 = 4.f;
  const double mI = 14.e-9f;
  const double mZA = 0.47999f; 
  return BetheBlochGeant(bg, rho, x0, x1, mI, mZA);
}

double GPUTPCTrackParam::ApproximateBetheBloch(double beta2)
{
  //------------------------------------------------------------------
  // This is an approximation of the Bethe-Bloch formula with
  // the density effect taken into account at beta*gamma > 3.5
  // (the approximation is reasonable only for solid materials)
  //------------------------------------------------------------------
  if (beta2 >= 1) {
    return 0;
  }

  if (beta2 / (1 - beta2) > 3.5f * 3.5f) {
    return 0.153e-3f / beta2 * (log(3.5f * 5940) + 0.5f * log(beta2 / (1 - beta2)) - beta2);
  }
  return 0.153e-3f / beta2 * (log(5940 * beta2 / (1 - beta2)) - beta2);
}

void GPUTPCTrackParam::CalculateFitParameters(GPUTPCTrackFitParam& par, double mass)
{
  //*!

  double qpt = GetPar(4);
  if (mC[14] >= 1.f) {
    qpt = 1.f / 0.35f;
  }

  double p2 = (1.f + GetPar(3) * GetPar(3));
  double k2 = qpt * qpt;
  double mass2 = mass * mass;
  double beta2 = p2 / (p2 + mass2 * k2);

  double pp2 = (k2 > 1.e-8f) ? p2 / k2 : 10000; // impuls 2

  par.bethe = BetheBlochGas( pp2/mass2);
  //par.bethe = ApproximateBetheBloch(pp2 / mass2);
  par.e = std::sqrt(pp2 + mass2);
  par.theta2 = 14.1f * 14.1f / (beta2 * pp2 * 1e6f);
  par.EP2 = par.e / pp2;

  // Approximate energy loss fluctuation (M.Ivanov)

  const double knst = 0.07f; // To be tuned.
  par.sigmadE2 = knst * par.EP2 * qpt;
  par.sigmadE2 = par.sigmadE2 * par.sigmadE2;

  par.k22 = (1.f + GetPar(3) * GetPar(3));
  par.k33 = par.k22 * par.k22;
  par.k43 = 0;
  par.k44 = GetPar(3) * GetPar(3) * k2;
}

bool GPUTPCTrackParam::CorrectForMeanMaterial(double xOverX0, double xTimesRho, const GPUTPCTrackFitParam& par)
{
  //------------------------------------------------------------------
  // This function corrects the track parameters for the crossed material.
  // "xOverX0"   - X/X0, the thickness in units of the radiation length.
  // "xTimesRho" - is the product length*density (g/cm^2).
  //------------------------------------------------------------------

  double& mC22 = mC[5];
  double& mC33 = mC[9];
  double& mC40 = mC[10];
  double& mC41 = mC[11];
  double& mC42 = mC[12];
  double& mC43 = mC[13];
  double& mC44 = mC[14];

  // Energy losses************************

  double dE = par.bethe * xTimesRho;
  if (std::abs(dE) > 0.3f * par.e) {
    return 0; // 30% energy loss is too much!
  }
  double corr = (1.f - par.EP2 * dE);
  if (corr < 0.3f || corr > 1.3f) {
    return 0;
  }

  SetPar(4, GetPar(4) * corr);
  mC40 *= corr;
  mC41 *= corr;
  mC42 *= corr;
  mC43 *= corr;
  mC44 *= corr * corr;
  mC44 += par.sigmadE2 * std::abs(dE);

  // Multiple scattering******************

  double theta2 = par.theta2 * std::abs(xOverX0);
  mC22 += theta2 * par.k22 * (1.f - GetPar(2)) * (1.f + GetPar(2));
  mC33 += theta2 * par.k33;
  mC43 += theta2 * par.k43;
  mC44 += theta2 * par.k44;

  return 1;
}

//*
//* Rotation
//*
bool GPUTPCTrackParam::Rotate(double alpha, double maxSinPhi)
{
  //* Rotate the coordinate system in XY on the angle alpha

  double cA = cos(alpha);
  double sA = sin(alpha);
  double x = X(), y = Y(), sP = SinPhi(), cP = GetCosPhi();
  double cosPhi = cP * cA + sP * sA;
  double sinPhi = -cP * sA + sP * cA;

  if (std::abs(sinPhi) > maxSinPhi || std::abs(cosPhi) < 1.e-2f || std::abs(cP) < 1.e-2f) {
    return 0;
  }

  double j0 = cP / cosPhi;
  double j2 = cosPhi / cP;

  SetX(x * cA + y * sA);
  SetY(-x * sA + y * cA);
  SetSignCosPhi(cosPhi);
  SetSinPhi(sinPhi);

  // double J[5][5] = { { j0, 0, 0,  0,  0 }, // Y
  //                      {  0, 1, 0,  0,  0 }, // Z
  //                      {  0, 0, j2, 0,  0 }, // SinPhi
  //                    {  0, 0, 0,  1,  0 }, // DzDs
  //                    {  0, 0, 0,  0,  1 } }; // Kappa
  // cout<<"alpha="<<alpha<<" "<<x<<" "<<y<<" "<<sP<<" "<<cP<<" "<<j0<<" "<<j2<<endl;
  // cout<<"      "<<mC[0]<<" "<<mC[1]<<" "<<mC[6]<<" "<<mC[10]<<" "<<mC[4]<<" "<<mC[5]<<" "<<mC[8]<<" "<<mC[12]<<endl;
  mC[0] *= j0 * j0;
  mC[1] *= j0;
  mC[3] *= j0;
  mC[6] *= j0;
  mC[10] *= j0;

  mC[3] *= j2;
  mC[4] *= j2;
  mC[5] *= j2 * j2;
  mC[8] *= j2;
  mC[12] *= j2;

  if (cosPhi < 0) {
    SetSinPhi(-SinPhi());
    SetDzDs(-DzDs());
    SetQPt(-QPt());
    mC[3] = -mC[3];
    mC[4] = -mC[4];
    mC[6] = -mC[6];
    mC[7] = -mC[7];
    mC[10] = -mC[10];
    mC[11] = -mC[11];
  }

  // cout<<"      "<<mC[0]<<" "<<mC[1]<<" "<<mC[6]<<" "<<mC[10]<<" "<<mC[4]<<" "<<mC[5]<<" "<<mC[8]<<" "<<mC[12]<<endl;
  return 1;
}

bool GPUTPCTrackParam::Rotate(double alpha, GPUTPCTrackLinearisation& t0, double maxSinPhi)
{
  //* Rotate the coordinate system in XY on the angle alpha

  double cA = cos(alpha);
  double sA = sin(alpha);
  double x0 = X(), y0 = Y(), sP = t0.SinPhi(), cP = t0.CosPhi();
  double cosPhi = cP * cA + sP * sA;
  double sinPhi = -cP * sA + sP * cA;

  if (std::abs(sinPhi) > maxSinPhi || std::abs(cosPhi) < 1.e-2f || std::abs(cP) < 1.e-2f) {
    return 0;
  }

  // double J[5][5] = { { j0, 0, 0,  0,  0 }, // Y
  //                    {  0, 1, 0,  0,  0 }, // Z
  //                    {  0, 0, j2, 0,  0 }, // SinPhi
  //                  {  0, 0, 0,  1,  0 }, // DzDs
  //                  {  0, 0, 0,  0,  1 } }; // Kappa

  double j0 = cP / cosPhi;
  double j2 = cosPhi / cP;
  double d[2] = {Y() - y0, SinPhi() - sP};

  SetX(x0 * cA + y0 * sA);
  SetY(-x0 * sA + y0 * cA + j0 * d[0]);
  t0.SetCosPhi(cosPhi);
  t0.SetSinPhi(sinPhi);

  SetSinPhi(sinPhi + j2 * d[1]);

  mC[0] *= j0 * j0;
  mC[1] *= j0;
  mC[3] *= j0;
  mC[6] *= j0;
  mC[10] *= j0;

  mC[3] *= j2;
  mC[4] *= j2;
  mC[5] *= j2 * j2;
  mC[8] *= j2;
  mC[12] *= j2;

  return 1;
}

bool GPUTPCTrackParam::Filter(double y, double z, double err2Y, double err2Z, double maxSinPhi, bool paramOnly)
{
  //* Add the y,z measurement with the Kalman filter

  double c00 = mC[0], c11 = mC[2], c20 = mC[3], c31 = mC[7], c40 = mC[10];

  err2Y += c00;
  err2Z += c11;

  double z0 = y - GetPar(0), z1 = z - GetPar(1);

  if (err2Y < 1.e-8f || err2Z < 1.e-8f) {
    return 0;
  }

  double mS0 = 1.f / err2Y;
  double mS2 = 1.f / err2Z;

  // K = CHtS

  double k00, k11, k20, k31, k40;

  k00 = c00 * mS0;
  k20 = c20 * mS0;
  k40 = c40 * mS0;

  k11 = c11 * mS2;
  k31 = c31 * mS2;

  double sinPhi = GetPar(2) + k20 * z0;

  if (maxSinPhi > 0 && std::abs(sinPhi) >= maxSinPhi) {
    return 0;
  }

  SetPar(0, GetPar(0) + k00 * z0);
  SetPar(1, GetPar(1) + k11 * z1);
  SetPar(2, sinPhi);
  SetPar(3, GetPar(3) + k31 * z1);
  SetPar(4, GetPar(4) + k40 * z0);
  if (paramOnly) {
    return true;
  }

  mNDF += 2;
  mChi2 += mS0 * z0 * z0 + mS2 * z1 * z1;

  mC[0] -= k00 * c00;
  mC[3] -= k20 * c00;
  mC[5] -= k20 * c20;
  mC[10] -= k40 * c00;
  mC[12] -= k40 * c20;
  mC[14] -= k40 * c40;

  mC[2] -= k11 * c11;
  mC[7] -= k31 * c11;
  mC[9] -= k31 * c31;

  return 1;
}

bool GPUTPCTrackParam::CheckNumericalQuality() const
{
  //* Check that the track parameters and covariance matrix are reasonable

  bool ok = std::isfinite(GetX()) && std::isfinite(mSignCosPhi) && std::isfinite(mChi2);

  const double* c = Cov();
  for (int i = 0; i < 15; i++) {
    ok = ok && std::isfinite(c[i]);
  }
  for (int i = 0; i < 5; i++) {
    ok = ok && std::isfinite(Par()[i]);
  }

  if (c[0] <= 0 || c[2] <= 0 || c[5] <= 0 || c[9] <= 0 || c[14] <= 0) {
    ok = 0;
  }
  if (c[0] > 5.f || c[2] > 5.f || c[5] > 2.f || c[9] > 2
      //|| ( std::abs( QPt() ) > 1.e-2 && c[14] > 2. )
  ) {
    ok = 0;
  }

  if (std::abs(SinPhi()) > GPUCA_MAX_SIN_PHI) {
    ok = 0;
  }
  if (std::abs(QPt()) > 1.f / 0.05f) {
    ok = 0;
  }
  if (ok) {
    ok = ok && (c[1] * c[1] <= c[2] * c[0]) && (c[3] * c[3] <= c[5] * c[0]) && (c[4] * c[4] <= c[5] * c[2]) && (c[6] * c[6] <= c[9] * c[0]) && (c[7] * c[7] <= c[9] * c[2]) && (c[8] * c[8] <= c[9] * c[5]) && (c[10] * c[10] <= c[14] * c[0]) && (c[11] * c[11] <= c[14] * c[2]) &&
         (c[12] * c[12] <= c[14] * c[5]) && (c[13] * c[13] <= c[14] * c[9]);
  }
  return ok;
}

#if !defined(GPUCA_GPUCODE)
#include <iostream>
#endif

void GPUTPCTrackParam::Print() const
{
  //* print parameters

#if !defined(GPUCA_GPUCODE)
  std::cout << "track: x=" << GetX() << " c=" << GetSignCosPhi() << ", P= " << GetY() << " " << GetZ() << " " << GetSinPhi() << " " << GetDzDs() << " " << GetQPt() << std::endl;
  std::cout << "errs2: " << GetErr2Y() << " " << GetErr2Z() << " " << GetErr2SinPhi() << " " << GetErr2DzDs() << " " << GetErr2QPt() << std::endl;
#endif
}
