#include "BEmcRecCEMC.h"
#include "BEmcCluster.h"

#include <cmath>

void BEmcRecCEMC::Tower2Global(float E, float xC, float yC,
                               float& xA, float& yA, float& zA)
{
  xA = 0;
  yA = 0;
  zA = 0;
}

void BEmcRecCEMC::CorrectEnergy(float Energy, float x, float y,
                                float* Ecorr)
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
}

void BEmcRecCEMC::CorrectECore(float Ecore, float x, float y, float* Ecorr)
{
  // Corrects the EM Shower Core Energy for attenuation in fibers,
  // long energy leakage and angle dependance
  //
  // (x,y) - impact position (cm) in Sector frame

  const float c0 = 0.950;  // For no threshold
  *Ecorr = Ecore / c0;
}

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
    printf("????? Something wrong in CorrectPosition: x=%f  dx=%f\n", x, x0 - ix0);
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
    printf("????? Something wrong in CorrectPosition: y=%f  dy=%f\n", y, y0 - iy0);
  }
}
