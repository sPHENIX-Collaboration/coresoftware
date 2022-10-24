#ifndef ANALYTICFIELDMODEL_H
#define ANALYTICFIELDMODEL_H

#include <TVector3.h>

#include <cmath>

class TFormula;

class AnalyticFieldModel
{
 public:
  double ifc_radius = NAN;
  double ofc_radius = NAN;
  double tpc_halfz = NAN;

  AnalyticFieldModel(float _ifc_radius, float _ofc_radius, float _z_max, float scalefactor = 1);
  //! delete copy ctor and assignment opertor (cppcheck)
  explicit AnalyticFieldModel(const AnalyticFieldModel &) = delete;
  AnalyticFieldModel &operator=(const AnalyticFieldModel &) = delete;

  TVector3 E(TVector3 pos);                     //field as a function of position
  double Rho(TVector3 pos);                     //charge density as a function of position
  TVector3 Eint(float zfinal, TVector3 start);  //field integral from start point to position zfinal.

 private:
  TFormula *vTestFunction1 = nullptr;
  TFormula *rhoTestFunction1 = nullptr;
  TFormula *erTestFunction1 = nullptr;
  TFormula *ePhiTestFunction1 = nullptr;
  TFormula *ezTestFunction1 = nullptr;
  TFormula *intErDzTestFunction1 = nullptr;

  TFormula *intEPhiDzTestFunction1 = nullptr;

  TFormula *intEzDzTestFunction1 = nullptr;
};

#endif
