// TPC HIT class
// Stores  one hit in TPC fiducial volume
// Author: Carlos Perez
#ifndef __TPCHIT_H__
#define __TPCHIT_H__

#include "vHit.h"

class TPCHit : public vHit {
 public:
  TPCHit();
  virtual ~TPCHit();
  Float_t GetDEnergy() {return fDEnergy;}
  void SetDEnergy(Float_t v) {fDEnergy=v;}

  Float_t GetR() {return X(0);}
  Float_t GetPhi() {return X(1);}
  Float_t GetZ() {return X(2);}
  Float_t GetL() {return X(3);}
  void SetR(Float_t v) {SetX(0,v);}
  void SetPhi(Float_t v) {SetX(1,v);}
  void SetZ(Float_t v) {SetX(2,v);}
  void SetL(Float_t v) {SetX(3,v);}

  void AddR(Float_t d) {AddX(0,d);}
  void AddPhi(Float_t d) {AddX(1,d);}
  void AddZ(Float_t d) {AddX(2,d);}

  virtual void CopyFrom(TPCHit *th);

 protected:
  Float_t fDEnergy; // Ionizing energy from step in keV
};

#endif
