#include "BEmcRecFEMC.h"
#include "BEmcCluster.h"
#include "BEmcProfile.h"

#include <cmath>
#include <iostream>

BEmcRecFEMC::BEmcRecFEMC()
//  : _emcprof(nullptr)
{
  Name("BEmcRecFEMC");
  SetPlanarGeometry();
}

BEmcRecFEMC::~BEmcRecFEMC()
{
  // one can delete null pointers
  //  delete _emcprof;
}

void BEmcRecFEMC::LoadProfile(const std::string& fname)
{
  _emcprof = new BEmcProfile(fname);
}

void BEmcRecFEMC::GetImpactThetaPhi(float xg, float yg, float zg, float& theta, float& phi)
{
  theta = std::atan(std::sqrt(xg * xg + yg * yg) / std::fabs(zg - fVz));
  phi = std::atan2(yg, xg);
}

/*
float BEmcRecFEMC::GetProb(vector<EmcModule> HitList, float ecl, float xg, float yg, float zg, float& chi2, int& ndf)
{
  chi2 = 0;
  ndf = 0;
  float prob = -1;

  float theta = atan(sqrt(xg * xg + yg * yg) / fabs(zg - fVz));
  float phi = atan2(yg, xg);
  if (_emcprof != nullptr) prob = _emcprof->GetProb(&HitList, fNx, ecl, theta, phi);

  return prob;
}
*/

void BEmcRecFEMC::CorrectShowerDepth(float E, float xA, float yA, float zA, float& xC, float& yC, float& zC)
{
  // For ala PHENIX PbSc modules
  /*
  const float DZ = 18;   // in cm, tower half length
  const float D = 13.3;  // in cm, shower depth at 1 GeV relative to tower face; obtained from GEANT
  const float X0 = 2.2;  // in cm; obtained from GEANT (should be ~ rad length)
  */
  // For E684-based modules
  const float DZ = 8;    // in cm, tower half length
  const float D = 4.6;   // in cm, shower depth at 1 GeV relative to tower face; obtained from GEANT
  const float X0 = 0.8;  // in cm; obtained from GEANT (should be ~ rad length)

  float logE = log(0.1);
  if (E > 0.1) logE = std::log(E);
  float zV = zA - fVz;
  float cosT = std::fabs(zV) / std::sqrt(xA * xA + yA * yA + zV * zV);

  zC = (zA - DZ) + (D + X0 * logE) * cosT;
  //  zC = zA; // !!!!!

  xC = xA;
  yC = yA;
}

void BEmcRecFEMC::CorrectEnergy(float Energy, float /*x*/, float /*y*/,
                                float& Ecorr)
{
  Ecorr = Energy;
}

float BEmcRecFEMC::GetImpactAngle(float /*e*/, float /*x*/, float /*y*/)
// Get impact angle, (x,y) - position in Sector frame (cm)
{
  /*
  float xVert, yVert, zVert;
  float vx, vy, vz;

  GlobalToSector( fVx, fVy, fVz, &xVert, &yVert, &zVert );
  vz = -zVert;
  vy = y - yVert;
  vx = x - xVert;
  // From this point X coord in sector frame is Z coord in Global Coord System !!!
  *sinT = sqrt((vx*vx+vy*vy)/(vx*vx+vy*vy+vz*vz));
  */
  return 0;
}

void BEmcRecFEMC::CorrectECore(float Ecore, float x, float y, float& Ecorr)
{
  // Corrects the EM Shower Core Energy for attenuation in fibers,
  // long energy leakage and angle dependance
  //
  // (x,y) - shower CG in tower units (not projected anywhere!)

  //  Ecorr = Ecore;
  float ec, ec2, corr;
  const float par1 = 0.938;
  const float par2 = 0.50;
  const float par3 = 0.067;

  Ecorr = Ecore;
  if (Ecore < 0.01) return;

  float xA, yA, zA;
  Tower2Global(Ecore / 0.91, x, y, xA, yA, zA);  // 0.91 - average correction
  float tanT = std::sqrt(xA * xA + yA * yA) / std::fabs(zA - fVz);
  corr = par1 * (1 - tanT * tanT * tanT * (par2 + par3 * std::log(Ecore)));
  ec = Ecore / corr;

  //  CorrectEnergy( ec, x, y, &ec2);
  ec2 = ec;  // !!!!! CorrectEnergy must be implemented
  Ecorr = ec2;
}

void BEmcRecFEMC::CorrectPosition(float Energy, float x, float y,
                                  float& xc, float& yc)
{
  // Corrects the Shower Center of Gravity for the systematic shift due to
  // limited tower size
  //
  // Everything here is in tower units.
  // (x,y) - CG position, (xc,yc) - corrected position

  float xZero, yZero, bx, by;
  float t, x0, y0;
  int ix0, iy0;

  xc = x;
  yc = y;
  //  return;

  if (Energy < 0.01) return;

  float xA, yA, zA;
  Tower2Global(Energy, x, y, xA, yA, zA);
  zA -= fVz;
  //  float sinTx = xA / sqrt(xA * xA + zA * zA);
  //  float sinTy = yA / sqrt(yA * yA + zA * zA);
  float sinTy = xA / std::sqrt(xA * xA + zA * zA);  // x is second index in here
  float sinTx = yA / std::sqrt(yA * yA + zA * zA);
  float sin2Tx = sinTx * sinTx;
  float sin2Ty = sinTy * sinTy;

  if (sinTx > 0)
    xZero = -0.417 * sinTx - 1.500 * sin2Tx;
  else
    xZero = -0.417 * sinTx + 1.500 * sin2Tx;

  if (sinTy > 0)
    yZero = -0.417 * sinTy - 1.500 * sin2Ty;
  else
    yZero = -0.417 * sinTy + 1.500 * sin2Ty;

  t = 0.98 + 0.98 * std::sqrt(Energy);
  t *= 0.7;  // Temp correction
  bx = 0.17 - 0.009 * std::log(Energy) + t * sin2Tx;
  by = 0.17 - 0.009 * std::log(Energy) + t * sin2Ty;

  //  xZero *= 0.5;
  //  yZero *= 0.5;

  x0 = x + xZero;
  ix0 = EmcCluster::lowint(x0 + 0.5);

  if (EmcCluster::ABS(x0 - ix0) <= 0.5)
  {
    x0 = (ix0 - xZero) + bx * asinh(2. * (x0 - ix0) * sinh(0.5 / bx));
    xc = x0;
  }
  else
  {
    xc = x;
    std::cout << "????? Something wrong in BEmcRecFEMC::CorrectPosition: x = "
              << x << ",  dx = " << x0 - ix0 << std::endl;
  }

  y0 = y + yZero;
  iy0 = EmcCluster::lowint(y0 + 0.5);

  if (EmcCluster::ABS(y0 - iy0) <= 0.5)
  {
    y0 = (iy0 - yZero) + by * asinh(2. * (y0 - iy0) * sinh(0.5 / by));
    yc = y0;
  }
  else
  {
    yc = y;
    std::cout << "????? Something wrong in BEmcRecFEMC::CorrectPosition: y = "
              << y << ",  dy = " << y0 - iy0 << std::endl;
  }
}
