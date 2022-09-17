#include "BEmcRecEEMC.h"

#include "BEmcCluster.h"
#include "BEmcProfile.h"

#include <calobase/RawTowerDefs.h>

#include <cmath>
#include <iostream>

BEmcRecEEMC::BEmcRecEEMC()
{
  Name("BEmcRecEEMC");
  SetPlanarGeometry();
}

void BEmcRecEEMC::LoadProfile(const std::string& fname)
{
  //  std::cout << "Info from BEmcRecEEMC::LoadProfile(): no shower profile evaluation is defined yet for EEMC" << std::endl;
  _emcprof = new BEmcProfile(fname);
}

void BEmcRecEEMC::GetImpactThetaPhi(float xg, float yg, float zg, float& theta, float& phi)
{
  theta = std::atan(std::sqrt(xg * xg + yg * yg) / std::fabs(zg - fVz));
  phi = std::atan2(yg, xg);
}

/*
float BEmcRecEEMC::GetProb(vector<EmcModule> HitList, float ecl, float xg, float yg, float zg, float& chi2, int& ndf)
// ecl, xg, yg, zg not used here
{
  chi2 = 0;
  ndf = 0;
  float prob = -1;

  float theta = atan(sqrt(xg*xg + yg*yg)/fabs(zg-fVz));
  float phi = atan2(yg,xg);
  if( _emcprof != nullptr ) prob = _emcprof->GetProb(&HitList,fNx,ecl,theta,phi);

  return prob;
}
*/

void BEmcRecEEMC::CorrectShowerDepth(float E, float xA, float yA, float zA, float& xC, float& yC, float& zC)
{
  // DZ, D and X0 should be negative for -Z direction
  // Crystal:[-10., -6.1, -0.91], Sci_glass:[-20., -16.76, -2.5]
  // Original one:[-9, -6.1, -0.91] also for the lead tungsten

  int C_ID = BEmcRec::GetCalotype();
  //  std::cout << "Success pass data(ID): " << C_ID << std::endl;

  if ((C_ID == RawTowerDefs::EEMC) || (C_ID == RawTowerDefs::EEMC_crystal))
  {
    const float DZ = -0.5 * BEmcRec::GetScinSize();  // in cm, tower half length
    const float D = -6.1;                            // in cm, shower depth at 1 GeV relative to tower face; obtained from GEANT
    const float X0 = -0.91;                          // in cm; obtained from GEANT (should be ~ rad length)

    //      std::cout << "Success pass data(size): " << DZ << std::endl << std::endl;

    float logE = log(0.1);
    if (E > 0.1) logE = std::log(E);
    float zV = zA - fVz;
    float cosT = std::fabs(zV) / std::sqrt(xA * xA + yA * yA + zV * zV);

    zC = (zA - DZ) + (D + X0 * logE) * cosT;  //Only the shower depth corrected
    //  zC = zA; // !!!!!

    xC = xA;  // Keep the x and y the same. The x and y correction is in another code
    yC = yA;
  }
  else if (C_ID == RawTowerDefs::EEMC_glass)
  {
    const float DZ = -0.5 * BEmcRec::GetScinSize();  // in cm, tower half length
    const float D = -15.25;                          // in cm, shower depth at 1 GeV relative to tower face; obtained from GEANT
    const float X0 = -2.275;                         // in cm; obtained from GEANT (should be ~ rad length)

    //      std::cout << "Success pass data(size): " << DZ << std::endl << std::endl;

    float logE = log(0.1);
    if (E > 0.1) logE = std::log(E);
    float zV = zA - fVz;
    float cosT = std::fabs(zV) / std::sqrt(xA * xA + yA * yA + zV * zV);

    zC = (zA - DZ) + (D + X0 * logE) * cosT;  //Only the shower depth corrected
    //  zC = zA; // !!!!!

    xC = xA;  // Keep the x and y the same. The x and y correction is in another code
    yC = yA;
  }
}

void BEmcRecEEMC::CorrectEnergy(float Energy, float /*x*/, float /*y*/,
                                float& Ecorr)
{
  Ecorr = Energy;
}

void BEmcRecEEMC::CorrectECore(float Ecore, float /*x*/, float /*y*/, float& Ecorr)
{
  // Corrects the EM Shower Core Energy for attenuation in fibers,
  // long energy leakage and angle dependance
  //
  // (x,y) - shower CG in tower units (not projected anywhere!)

  Ecorr = Ecore;
}

float BEmcRecEEMC::GetImpactAngle(float /*e*/, float /*x*/, float /*y*/)
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

void BEmcRecEEMC::CorrectPosition(float Energy, float x, float y,
                                  float& xc, float& yc)
{
  // Corrects the Shower Center of Gravity for the systematic shift due to
  // the limited tower size
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
    xZero = -0.6 * sinTx - 1.500 * sin2Tx;
  else
    xZero = -0.6 * sinTx + 1.500 * sin2Tx;

  if (sinTy > 0)
    yZero = -0.6 * sinTy - 1.5 * sin2Ty;
  else
    yZero = -0.6 * sinTy + 1.5 * sin2Ty;

  yZero = -yZero;  // Because tower index in X decreases with increasing X in EEMC (y is actually x here!)

  t = 0.98 + 0.98 * std::sqrt(Energy);
  t *= 0.7;  // Temp correction
  bx = 0.20 - 0.009 * std::log(Energy) + t * sin2Tx;
  by = 0.20 - 0.009 * std::log(Energy) + t * sin2Ty;

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
    std::cout << "????? Something wrong in BEmcRecEEMC::CorrectPosition: x = "
              << x << ", dx = " << x0 - ix0 << std::endl;
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
    std::cout << "????? Something wrong in BEmcRecEEMC::CorrectPosition: y = "
              << y << ",  dy = " << y0 - iy0 << std::endl;
  }
}
