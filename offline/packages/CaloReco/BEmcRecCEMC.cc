#include "BEmcRecCEMC.h"

#include "BEmcCluster.h"
#include "BEmcProfile.h"

#include <cmath>
#include <iostream>

BEmcRecCEMC::BEmcRecCEMC()
//  : _emcprof(nullptr)
{
  Name("BEmcRecCEMC");
  SetCylindricalGeometry();
}

BEmcRecCEMC::~BEmcRecCEMC()
{
  // you can delete null pointers
  //  delete _emcprof;
}

void BEmcRecCEMC::LoadProfile(const std::string& fname)
{
  //  std::cout << "Infor from BEmcRecCEMC::LoadProfile(): no external file used for shower profile evaluation in CEMC" << std::endl;
  _emcprof = new BEmcProfile(fname);
}

void BEmcRecCEMC::GetImpactThetaPhi(float xg, float yg, float zg, float& theta, float& phi)
{
  theta = 0;
  phi = 0;

  //  float theta = atan(sqrt(xg*xg + yg*yg)/fabs(zg-fVz));
  float rg = std::sqrt(xg * xg + yg * yg);
  float theta_twr;
  if (std::fabs(zg) <= 15)
    theta_twr = 0;
  else if (zg > 15)
    theta_twr = std::atan2(zg - 15, rg);
  else
    theta_twr = std::atan2(zg + 15, rg);
  float theta_tr = std::atan2(zg - fVz, rg);
  theta = std::fabs(theta_tr - theta_twr);
  //  phi = atan2(yg,xg);
}

/*
float BEmcRecCEMC::GetProb(vector<EmcModule> HitList, float ecl, float xg, float yg, float zg, float& chi2, int& ndf)
{
  chi2 = 0;
  ndf = 0;
  float prob = -1;

  //  float theta = atan(sqrt(xg*xg + yg*yg)/fabs(zg-fVz));
  float rg = sqrt(xg * xg + yg * yg);
  float theta_twr;
  if (fabs(zg) <= 15)
    theta_twr = 0;
  else if (zg > 15)
    theta_twr = atan2(zg - 15, rg);
  else
    theta_twr = atan2(zg + 15, rg);
  float theta_tr = atan2(zg - fVz, rg);
  float theta = fabs(theta_tr - theta_twr);

  float phi = atan2(yg, xg);
  if (_emcprof != nullptr) prob = _emcprof->GetProb(&HitList, fNx, ecl, theta, phi);

  return prob;
}
*/
/*
float BEmcRecCEMC::GetProb(vector<EmcModule> HitList, float et, float xg, float yg, float zg, float& chi2, int& ndf)
// et, xg, yg, zg not used here
{
  const float thresh = 0.01;
  const int DXY = 3;  // 2 is for 5x5 matrix; 3 for 7x7 matrix
  const int Nmax = 1000;
  float ee[Nmax];
  int iyy[Nmax];
  int izz[Nmax];

  int ich;
  vector<EmcModule>::iterator ph = HitList.begin();

  chi2 = 0;
  ndf = 0;

  int nn = 0;

  while (ph != HitList.end())
  {
    ee[nn] = ph->amp;
    if (ee[nn] > thresh)
    {
      ich = ph->ich;
      izz[nn] = ich % fNx;
      iyy[nn] = ich / fNx;
      nn++;
      if (nn >= Nmax)
      {
      std::cout << "BEmcRec::GetProb: Cluster size is too big. Skipping the rest of the towers" << std::endl;
        break;
      }
    }  // if( ee[nn]
    ++ph;
  }  // while( ph

  if (nn <= 0) return -1;

  int iy0 = -1, iz0 = -1;
  float emax = 0;

  for (int i = 0; i < nn; i++)
  {
    if (ee[i] > emax)
    {
      emax = ee[i];
      iy0 = iyy[i];
      iz0 = izz[i];
    }
  }

  if (emax <= 0) return -1;

  int id;
  float etot = 0;
  float sz = 0;
  float sy = 0;

  for (int idz = -DXY; idz <= DXY; idz++)
  {
    for (int idy = -DXY; idy <= DXY; idy++)
    {
      id = GetTowerID(iy0 + idy, iz0 + idz, nn, iyy, izz, ee);
      if (id >= 0)
      {
        etot += ee[id];
        sz += ee[id] * (iz0 + idz);
        sy += ee[id] * (iy0 + idy);
      }
    }
  }
  float zcg = sz / etot;  // Here cg allowed to be out of range
  float ycg = sy / etot;
  int iz0cg = int(zcg + 0.5);
  int iy0cg = int(ycg + 0.5);
  float ddz = fabs(zcg - iz0cg);
  float ddy = fabs(ycg - iy0cg);

  int isz = 1;
  if (zcg - iz0cg < 0) isz = -1;
  int isy = 1;
  if (ycg - iy0cg < 0) isy = -1;

  // 4 central towers: 43
  //                   12
  // Tower 1 - central one
  float e1, e2, e3, e4;
  e1 = e2 = e3 = e4 = 0;
  id = GetTowerID(iy0cg, iz0cg, nn, iyy, izz, ee);
  if (id >= 0) e1 = ee[id];
  id = GetTowerID(iy0cg, iz0cg + isz, nn, iyy, izz, ee);
  if (id >= 0) e2 = ee[id];
  id = GetTowerID(iy0cg + isy, iz0cg + isz, nn, iyy, izz, ee);
  if (id >= 0) e3 = ee[id];
  id = GetTowerID(iy0cg + isy, iz0cg, nn, iyy, izz, ee);
  if (id >= 0) e4 = ee[id];

  float e1t = (e1 + e2 + e3 + e4) / etot;
  float e2t = (e1 + e2 - e3 - e4) / etot;
  float e3t = (e1 - e2 - e3 + e4) / etot;
  float e4t = (e3) / etot;
  //  float e5t = (e2+e4)/etot;

  float rr = sqrt((0.5 - ddz) * (0.5 - ddz) + (0.5 - ddy) * (0.5 - ddy));

  float c1, c2, c11;

  float logE = log(etot);

  // e1 energy is the most effective for PID if properly tuned !
  // Discrimination power is very sensitive to paramter c1: the bigger it is
  // the better discrimination;
  c1 = 0.95;
  c2 = 0.0066364 * logE + 0.00466667;
  if (c2 < 0) c2 = 0;
  float e1p = c1 - c2 * rr * rr;
  c1 = 0.034 - 0.01523 * logE + 0.0029 * logE * logE;
  float err1 = c1;

  // For e2
  c1 = 0.00844086 + 0.00645359 * logE - 0.00119381 * logE * logE;
  if (etot > 15) c1 = 0.00844086 + 0.00645359 * log(15.) - 0.00119381 * log(15.) * log(15.);  // Const at etot>15GeV
  if (c1 < 0) c1 = 0;
  c2 = 3.9;                                                      // Fixed
  float e2p = sqrt(c1 + 0.25 * c2) - sqrt(c1 + c2 * ddy * ddy);  // =0 at ddy=0.5

  c1 = 0.0212333 + 0.0420473 / etot;
  c2 = 0.090;  // Fixed
  float err2 = c1 + c2 * ddy;
  if (ddy > 0.3) err2 = c1 + c2 * 0.3;  // Const at ddy>0.3

  // For e3
  c1 = 0.0107857 + 0.0056801 * logE - 0.000892016 * logE * logE;
  if (etot > 15) c1 = 0.0107857 + 0.0056801 * log(15.) - 0.000892016 * log(15.) * log(15.);  // Const at etot>15GeV
  if (c1 < 0) c1 = 0;
  c2 = 3.9;                                                      // Fixed
  float e3p = sqrt(c1 + 0.25 * c2) - sqrt(c1 + c2 * ddz * ddz);  // =0 at ddz=0.5

  //  c1 = 0.0200 + 0.042/etot;
  c1 = 0.0167 + 0.058 / etot;
  c2 = 0.090;  // Fixed
  float err3 = c1 + c2 * ddz;
  if (ddz > 0.3) err3 = c1 + c2 * 0.3;  // Const at ddz>0.3

  // For e4
  float e4p = 0.25 - 0.668 * rr + 0.460 * rr * rr;
  c11 = 0.171958 + 0.0142421 * logE - 0.00214827 * logE * logE;
  //  c11 = 0.171085 + 0.0156215*logE - -0.0025809*logE*logE;
  float err4 = 0.102 - 1.43 * c11 * rr + c11 * rr * rr;  // Min is set to x=1.43/2.
  err4 *= 1.1;

  chi2 = 0.;
  chi2 += (e1p - e1t) * (e1p - e1t) / err1 / err1;
  chi2 += (e2p - e2t) * (e2p - e2t) / err2 / err2;
  chi2 += (e3p - e3t) * (e3p - e3t) / err3 / err3;
  chi2 += (e4p - e4t) * (e4p - e4t) / err4 / err4;
  ndf = 4;

  //  chi2 /= 1.1;
  float prob = TMath::Prob(chi2, ndf);

  return prob;
}
*/

void BEmcRecCEMC::CorrectShowerDepth(float E, float xA, float yA, float zA, float& xC, float& yC, float& zC)
{
  /*
  xC = xA;
  yC = yA;
  zC = zA;
  */

  float logE = log(0.1);
  if (E > 0.1) logE = std::log(E);

  // Rotate by phi (towers are tilted by a fixed angle in phi by ~9 deg?)
  // Just tuned from sim data
  float phi = 0.002 - 0.001 * logE;
  xC = xA * std::cos(phi) - yA * std::sin(phi);
  yC = xA * std::sin(phi) + yA * std::cos(phi);

  // Correction in z
  // Just tuned for sim data ... don't fully understand why it works like that
  float rA = std::sqrt(xA * xA + yA * yA);
  //  float theta_twr = GetTowerTheta(xA,yA,zA);
  float theta_twr;
  if (std::fabs(zA) <= 15)
    theta_twr = 0;
  else if (zA > 15)
    theta_twr = std::atan2(zA - 15, rA);
  else
    theta_twr = std::atan2(zA + 15, rA);

  float theta_tr = std::atan2(zA - fVz, rA);
  float L = -1.3 + 0.7 * logE;  // Shower CG in long. direction
  float dz = L * std::sin(theta_tr - theta_twr) / std::cos(theta_twr);

  dz -= fVz * 0.10;

  zC = zA - dz;

  return;
}

void BEmcRecCEMC::CorrectEnergy(float Energy, float /*x*/, float /*y*/,
                                float& Ecorr)
{
  // Corrects the EM Shower Energy for attenuation in fibers and
  // long energy leakage
  //
  // (x,y) - impact position (cm) in Sector frame
  /*
  float sinT;
  float att, leak, corr;
  const float leakPar = 0.0033; // parameter from fit
  const float attPar = 120; // Attenuation in module (cm)
  const float X0 = 2; // radiation length (cm)

  *Ecorr = Energy;
  if( Energy < 0.01 ) return;

  GetImpactAngle(x, y, &sinT); // sinT not used so far
  leak = 2-sqrt(1+leakPar*log(1+Energy)*log(1+Energy));
  att = exp(log(Energy)*X0/attPar);
  corr = leak*att;
  *Ecorr = Energy/corr;
  */
  Ecorr = Energy;
}

void BEmcRecCEMC::CorrectECore(float Ecore, float /*x*/, float /*y*/, float& Ecorr)
{
  // Corrects the EM Shower Core Energy for attenuation in fibers,
  // long energy leakage and angle dependance
  //
  // (x,y) - impact position (cm) in Sector frame

  const float c0 = 0.950;  // For no threshold
  Ecorr = Ecore / c0;
}

void BEmcRecCEMC::CorrectPosition(float Energy, float x, float y,
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

  if (Energy < 0.01) return;
  /*
  float xA, yA, zA;
  Tower2Global(Energy, x, y, xA, yA, zA);
  zA -= fVz;
  float sinTx = xA / sqrt(xA * xA + zA * zA);
  float sinTy = yA / sqrt(yA * yA + zA * zA);
  */
  float sinTx = 0;
  float sinTy = 0;

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
  bx = 0.15 + t * sin2Tx;
  by = 0.15 + t * sin2Ty;

  x0 = x + xZero;
  ix0 = EmcCluster::lowint(x0 + 0.5);

  if (EmcCluster::ABS(x0 - ix0) <= 0.5)
  {
    x0 = (ix0 - xZero) + bx * asinh(2. * (x0 - ix0) * sinh(0.5 / bx));
  }
  else
  {
    x0 = x;
    std::cout << "????? Something wrong in BEmcRecCEMC::CorrectPosition: x = "
              << x << " dx = " << x0 - ix0 << std::endl;
  }

  // Correct for phi bias within module of 8 towers
  int ix8 = int(x + 0.5) / 8;
  float x8 = x + 0.5 - ix8 * 8 - 4;  // from -4 to +4
  float dx = 0.10 * x8 / 4.;
  if (std::fabs(x8) > 3.3) dx = 0;  // Don't correct near the module edge
  //  dx = 0;

  xc = x0 - dx;
  while (xc < -0.5) xc += float(fNx);
  while (xc >= fNx - 0.5) xc -= float(fNx);

  y0 = y + yZero;
  iy0 = EmcCluster::lowint(y0 + 0.5);

  if (EmcCluster::ABS(y0 - iy0) <= 0.5)
  {
    y0 = (iy0 - yZero) + by * asinh(2. * (y0 - iy0) * sinh(0.5 / by));
  }
  else
  {
    y0 = y;
    std::cout << "????? Something wrong in BEmcRecCEMC::CorrectPosition: y = "
              << y << "dy = " << y0 - iy0 << std::endl;
  }
  yc = y0;
}

/*
void BEmcRecCEMC::CorrectPosition(float Energy, float x, float y,
                                  float* pxc, float* pyc)
{
  // Corrects the Shower Center of Gravity for the systematic error due to
  // the limited tower size and angle shift
  //
  // Everything here is in cell units.
  // (x,y) - CG position, (*pxc,*pyc) - corrected position

  float xShift, yShift, xZero, yZero, bx, by;
  float t, x0, y0;
  int ix0, iy0;
  //  int signx, signy;

  const float Xrad = 0.3;  // !!!!! Need to put correct value
  const float Remc = 90.;  // EMCal inner radius. !!!!! Should be obtained from geometry container

  *pxc = x;
  *pyc = y;
  //  return;

  SetProfileParameters(0, Energy, x, y);
  // if( fSinTx >= 0 ) signx =  1;
  // else 	   signx = -1;
  // if( fSinTy >= 0 ) signy =  1;
  // else 	   signy = -1;
  t = 5.0 + 1.0 * log(Energy);         // In Rad Length units
  t *= (Xrad / Remc / GetModSizex());  // !!!!!
  xShift = t * fSinTx;
  yShift = t * fSinTy;
  // xZero=xShift-(0.417*EmcCluster::ABS(fSinTx)+1.500*fSinTx*fSinTx)*signx;
  // yZero=yShift-(0.417*EmcCluster::ABS(fSinTy)+1.500*fSinTy*fSinTy)*signy;
  xZero = xShift;                  // ...Somehow this works better !!!!!
  yZero = yShift;                  // ...Somehow this works better !!!!!
  t = 0.98 + 0.98 * sqrt(Energy);  // !!!!! Still from PHENIX
  bx = 0.15 + t * fSinTx * fSinTx;
  by = 0.15 + t * fSinTy * fSinTy;

  x0 = x;
  x0 = x0 - xShift + xZero;
  ix0 = EmcCluster::lowint(x0 + 0.5);
  if (EmcCluster::ABS(x0 - ix0) <= 0.5)
  {
    x0 = (ix0 - xZero) + bx * asinh(2. * (x0 - ix0) * sinh(0.5 / bx));
    *pxc = x0;
  }
  else
  {
    *pxc = x - xShift;
    std::cout << "????? Something wrong in CorrectPosition: x = "
         << x << " dx = " << x0 - ix0 << std::endl;
  }

  y0 = y;
  y0 = y0 - yShift + yZero;
  iy0 = EmcCluster::lowint(y0 + 0.5);
  if (EmcCluster::ABS(y0 - iy0) <= 0.5)
  {
    y0 = (iy0 - yZero) + by * asinh(2. * (y0 - iy0) * sinh(0.5 / by));
    *pyc = y0;
  }
  else
  {
    *pyc = y - yShift;
    std::cout << "????? Something wrong in CorrectPosition: y = "
         << y << " dy = " << y0 - iy << std::endl;
  }
}
*/
