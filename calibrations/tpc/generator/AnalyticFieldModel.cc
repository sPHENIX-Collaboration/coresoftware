
#include "AnalyticFieldModel.h"
#include "TFormula.h"
#include "TVector3.h"

#include <cstdio>

AnalyticFieldModel::AnalyticFieldModel(float _ifc_radius, float _ofc_radius, float _z_max, float scalefactor)
{
  double ifc_radius = _ifc_radius;
  double ofc_radius = _ofc_radius;
  double tpc_halfz = _z_max;

  double sum = ifc_radius + ofc_radius;   //338 in ALICE, [3] in args
  double prod = ifc_radius * ofc_radius;  //21250.75 in ALICE [4] in args
  double diff = ofc_radius - ifc_radius;

  double a = ofc_radius * ofc_radius;
  a *= (diff);
  a *= (diff);
  //a =  a/1000.0;
  a = 1 / a;
  a *= scalefactor;
  double b = 0.5;
  double c = 1.0 / (((tpc_halfz) / 2.0) * ((tpc_halfz) / 2.0));
  double d = sum;
  double e = prod;

  vTestFunction1 = new TFormula("f1", "[0]*(x^4 - [3] *x^3 + [4] * x^2)*cos([1]* y)^2*exp(-1* [2] * z^2)");
  rhoTestFunction1 = new TFormula("ff1", "[0]*(((16.0 * x^2 - 9.0 * [3] * x + 4.0*[4]) *cos([1] * y)^2 * exp(-1 *[2]*z^2)) - ((x^2 -  [3] * x + [4]) * 2 * [1]^2 * cos(2 * [1] * y) * exp(-1 *[2]*z^2)) + ((x^4 -  [3] * x^3 + [4] * x^2) * cos([1] * y)^2 * (4*[2]^2*z^2 - 2 * [2]) * exp(-1 *[2]*z^2)))");

  erTestFunction1 = new TFormula("er", " [0]*(4*x^3 - 3 * [3] *x^2 + 2 * [4] * x)*cos([1]* y)^2*exp(-1* [2] * z^2)");
  ePhiTestFunction1 = new TFormula("ePhi",
                                   "  [0]*(x^3 - [3] *x^2 +  [4] * x)* -1  * [1] * sin(2 * [1]* y)*exp(-1* [2] * z^2)");
  ezTestFunction1 = new TFormula("ez",
                                 " [0]*(x^4 - [3] *x^3 + [4] * x^2)*cos([1]* y)^2*-1*2*[2]*z*exp(-1* [2] * z^2)");

  intErDzTestFunction1 = new TFormula("intErDz",
                                      " [0]*(4*x^3 - 3 * [3] *x^2 + 2 * [4] * x)*cos([1]* y)^2*((sqrt(pi)*TMath::Erf(sqrt([2]) * z))/(2 * sqrt([2]))) ");
  intEPhiDzTestFunction1 = new TFormula("intEPhiDz",
                                        "[0]* (x^3 - [3] *x^2 +  [4] * x)* -1  * [1] * sin(2 * [1]* y)*((sqrt(pi)*TMath::Erf(sqrt([2]) * z))/(2 * sqrt([2])))");
  intEzDzTestFunction1 = new TFormula("intEzDz",
                                      "[0]* (x^4 - [3] *x^3 + [4] * x^2)*cos([1]* y)^2*exp(-1* [2] * z^2)");

  printf("Setting Analytic Formula, variables:\n");
  printf("ifc=%f\tofc=%f\tdelz=%f\ndiff=%f\tscale=%f\n", ifc_radius, ofc_radius, tpc_halfz, diff, scalefactor);
  printf("a=%E\nb=%E\nc=%E\nd=%f\ne=%f\n", a, b, c, d, e);

  vTestFunction1->SetParameters(a, b, c, d, e);
  rhoTestFunction1->SetParameters(a, b, c, d, e);
  //printf("rho value at  (rmid,1,zmid)=%f\n",

  erTestFunction1->SetParameters(-a, b, c, d, e);
  ePhiTestFunction1->SetParameters(-a, b, c, d, e);
  ezTestFunction1->SetParameters(-a, b, c, d, e);
  intErDzTestFunction1->SetParameters(-a, b, c, d, e);
  intEPhiDzTestFunction1->SetParameters(-a, b, c, d, e);
  intEzDzTestFunction1->SetParameters(-a, b, c, d, e);
  return;
}

TVector3 AnalyticFieldModel::E(TVector3 pos)
{  //field as a function of position
  //in rhat phihat zhat coordinates: (at phi=0, phi is the +Y position, Perp is the +X direction and Z is Z)
  TVector3 ret(erTestFunction1->Eval(pos.Perp(), pos.Phi(), pos.Z()),
               ePhiTestFunction1->Eval(pos.Perp(), pos.Phi(), pos.Z()),
               ezTestFunction1->Eval(pos.Perp(), pos.Phi(), pos.Z()));
  //now rotate this to the position we evaluated it at, to match the global coordinate system.
  ret.RotateZ(pos.Phi());
  return ret;
}

double AnalyticFieldModel::Rho(TVector3 pos)
{                                             //charge density as a function of position
  const double alice_chargescale = 8.85e-14;  //their rho has charge density in units of C/cm^3 /eps0.  This is eps0 in (V*cm)/C units so that I can multiple by the volume in cm^3 to get Q in C.
  //at phi=0, phi is the +Y position, Perp is the +X direction and Z is Z.
  return alice_chargescale * rhoTestFunction1->Eval(pos.Perp(), pos.Phi(), pos.Z());
}

TVector3 AnalyticFieldModel::Eint(float zfinal, TVector3 pos)
{  //field integral from 'pos' to z-position zfinal.
  //in rhat phihat zhat coordinates: (at phi=0, phi is the +Y position, Perp is the +X direction and Z is Z)
  TVector3 eintI, eintF;
  eintI.SetXYZ(intErDzTestFunction1->Eval(pos.Perp(), pos.Phi(), pos.Z()),
               intEPhiDzTestFunction1->Eval(pos.Perp(), pos.Phi(), pos.Z()),
               intEzDzTestFunction1->Eval(pos.Perp(), pos.Phi(), pos.Z()));
  eintF.SetXYZ(intErDzTestFunction1->Eval(pos.Perp(), pos.Phi(), zfinal),
               intEPhiDzTestFunction1->Eval(pos.Perp(), pos.Phi(), zfinal),
               intEzDzTestFunction1->Eval(pos.Perp(), pos.Phi(), zfinal));

  TVector3 ret = eintF - eintI;
  // printf("Integrating z=%E to z=%E, delz=%E.  Before rotation, field integrals: (xyz)  (%E,%E,%E) to (%E,%E,%E), diff=(%E,%E,%E)\n",pos.Z(),zfinal,zfinal-pos.Z(),eintI.X(),eintI.Y(),eintI.Z(),eintF.X(),eintF.Y(),eintF.Z(),ret.X(),ret.Y(),ret.Z());
  //now rotate this to the position we evaluated it at, to match the global coordinate system.
  ret.RotateZ(pos.Phi());
  return ret;
}
