#include "BEmcRecEEMC.h"

#include "BEmcCluster.h"
#include "BEmcProfile.h"

#include <cmath>
#include <iostream>

using namespace std;

BEmcRecEEMC::BEmcRecEEMC() : _emcprof(nullptr) 
{
  Name("BEmcRecEEMC");
  SetPlanarGeometry();
}

BEmcRecEEMC::~BEmcRecEEMC()
{
// you can delete null pointers
  delete _emcprof;
}

void BEmcRecEEMC::LoadProfile(const std::string &fname)
{
  cout << "Info from BEmcRecEEMC::LoadProfile(): no shower profile evaluation is defined yet for EEMC" << endl;
}

float BEmcRecEEMC::GetProb(vector<EmcModule> HitList, float et, float xg, float yg, float zg, float& chi2, int& ndf)
// et, xg, yg, zg not used here
{
  chi2 = 0;
  ndf = 0;
  float prob = -1;
  return prob;
}

void BEmcRecEEMC::CorrectShowerDepth(float E, float xA, float yA, float zA, float& xC, float& yC, float& zC )
{
  // DZ, D and X0 should be negative for -Z direction
  const float DZ = -9;     // in cm, tower half length
  const float D = -6.1;    // in cm, shower depth at 1 GeV relative to tower face; obtained from GEANT
  const float X0 = -0.91;  // in cm; obtained from GEANT (should be ~ rad length)

  float logE = log(0.1);
  if (E > 0.1) logE = log(E);
  float zV = zA - fVz;
  float cosT = fabs(zV) / sqrt(xA * xA + yA * yA + zV * zV);

  zC = (zA - DZ) + (D + X0 * logE) * cosT;
  //  zC = zA; // !!!!!

  xC = xA;
  yC = yA;
}

void BEmcRecEEMC::CorrectEnergy(float Energy, float x, float y,
                                float* Ecorr)
{
  *Ecorr = Energy;
}

void BEmcRecEEMC::CorrectECore(float Ecore, float x, float y, float* Ecorr)
{
  // Corrects the EM Shower Core Energy for attenuation in fibers,
  // long energy leakage and angle dependance
  //
  // (x,y) - shower CG in tower units (not projected anywhere!)

  *Ecorr = Ecore;
}

float BEmcRecEEMC::GetImpactAngle(float e, float x, float y)
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
  float sinTy = xA / sqrt(xA * xA + zA * zA); // x is second index in here
  float sinTx = yA / sqrt(yA * yA + zA * zA);
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

  t = 0.98 + 0.98 * sqrt(Energy);
  t *= 0.7;  // Temp correction
  bx = 0.20 - 0.009 * log(Energy) + t * sin2Tx;
  by = 0.20 - 0.009 * log(Energy) + t * sin2Ty;

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
    cout << "????? Something wrong in BEmcRecEEMC::CorrectPosition: x = "
	 << x << ", dx = " << x0 - ix0 << endl;
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
    cout << "????? Something wrong in BEmcRecEEMC::CorrectPosition: y = "
	 << y << ",  dy = " << y0 - iy0 << endl;
  }
}
