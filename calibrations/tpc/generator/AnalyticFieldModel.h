//#pragma once
#ifndef ANALYTICFIELDMODEL_H
#define ANALYTICFIELDMODEL_H
#include "TVector3.h"
#include "TFormula.h"

class AnalyticFieldModel{
public:
  double ifc_radius;
  double ofc_radius;
  double tpc_halfz;

  


 
  
  AnalyticFieldModel(float _ifc_radius, float _ofc_radius, float _z_max, float scalefactor=1);

  TVector3 E(TVector3 pos);//field as a function of position
  double Rho(TVector3 pos);//charge density as a function of position
  TVector3 Eint(float zfinal, TVector3 start);//field integral from start point to position zfinal.

 private:
  
  TFormula *vTestFunction1;
  TFormula *rhoTestFunction1;
  TFormula *erTestFunction1;
  TFormula *ePhiTestFunction1;
  TFormula *ezTestFunction1;
  TFormula *intErDzTestFunction1;
			
  TFormula *intEPhiDzTestFunction1;
			
  TFormula *intEzDzTestFunction1;
			

};

#endif
