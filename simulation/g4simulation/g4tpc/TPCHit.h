// Points to one hit in TPC fiducial volume
// wraps any existing data structure
// Author: Carlos Perez
#ifndef __TPCHIT_H__
#define __TPCHIT_H__

#include <cmath>
#include "g4main/PHG4Hit.h" // base clase of data container

class TPCHit {
 public:
  TPCHit();
  virtual ~TPCHit();
  void Assign( PHG4Hit *hit ) {fHit = hit;}
  virtual void CopyFrom( TPCHit *hit ) {Assign( hit->GetDatum() );}
  void Clear() {fHit=NULL;}

  float GetX() {return fHit->get_avg_x();}
  float GetY() {return fHit->get_avg_y();}
  float GetZ() {return fHit->get_avg_z();}
  float GetR();
  float GetPhi();
  float GetT() {return fHit->get_avg_t();}
  float GetL() {return std::sqrt( std::pow(fHit->get_x(0)-fHit->get_x(1),2) +
				  std::pow(fHit->get_y(0)-fHit->get_y(1),2) +
				  std::pow(fHit->get_z(0)-fHit->get_z(1),2) );}
  float GetDEnergy() {return fHit->get_eion()*1e6;} // in keV
  int GetTrack() {return fHit->get_trkid();}

  PHG4Hit* GetDatum() {return fHit;}

 protected:
  PHG4Hit *fHit; //! not owned

  ClassDef(TPCHit,1);
};

#endif
