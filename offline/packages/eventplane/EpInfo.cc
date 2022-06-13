#include "EpInfo.h"
#include "TMath.h"
#include <iostream>

EpInfo::EpInfo(){
  Reset(); 
}

void EpInfo::Reset()
{
  for (int iorder=0; iorder<_EpOrderMax; iorder++){
    for (int xy=0; xy<2; xy++){
      QrawOneSide[iorder][xy] = 0.0;
    }
  }

  for (int iorder=0; iorder<_EpOrderMax; iorder++){
    PsiRaw[iorder] = -999.0;
  }
}

// ===================== Access to Q-vectors ==========================

//------------------------------ Raw Q --------------------------------
TVector2 EpInfo::RawQ(int order){
  if (ArgumentOutOfBounds(order)){
    TVector2 crap(-999,-999);
    return crap;
  }
  TVector2 q(QrawOneSide[order-1][0],QrawOneSide[order-1][1]);
  return q;
}

// --------------------- Wheel sum-of-weights, raw ----------------------
double EpInfo::SWRaw(int order){
  if (ArgumentOutOfBounds(order)) return -999;
  return WheelSumWeightsRaw[order-1];
}


// ===================== Access to Event-plane angles ====================

//------------------------- raw EP angles --------------------------------
double EpInfo::RawPsi(int order){
  if (ArgumentOutOfBounds(order)) return -999;
  return Range(PsiRaw[order-1],order);
}
//-----------------------------------------------------------------------



//----- Simple method to put angles in a convenient range: (0,2pi/n) ----
double EpInfo::Range(double psi,int order){
  if (ArgumentOutOfBounds(order)) return -999;
  double wrap = (2.0*M_PI)/(double)order;
  if (psi<0.0) psi += (1.0+(int)(fabs(psi)/wrap))*wrap;
  else{ if (psi>wrap) psi -= ((int)(psi/wrap))*wrap;}
  return psi;
}
  
//--------- protection against argument out-of-bounds -------
bool EpInfo::ArgumentOutOfBounds(int order){
  if ((order<1)||(order>_EpOrderMax)){
    std::cout << "\n *** Invalid order requested ***\n";
    std::cout << "  order must be between 1 (for first-order EP) and " << _EpOrderMax
	      << ".    To change the upuper limit, edit StEpdUtil/EpInfo.h\n";
    std::cout << "  I will now return you an invalid result.  Have a nice day\n";
    return true;
  }
  return false;
}

