#ifndef __TPCCLOUD_H__
#define __TPCCLOUD_H__

#include <TMath.h>
#include <TPCbase/TPCHit.h>

class TPCCloud : public TPCHit {
 public:
  TPCCloud();
  virtual ~TPCCloud();
  void SetMST(float v) {fMST=v;}
  void SetMSL0(float v) {fMSL0=v;}
  void SetMSL1(float v) {fMSL1=v;}
  float GetRMST() {return TMath::Sqrt(fMST);}
  float GetRMSL0() {return TMath::Sqrt(fMSL0);}
  float GetRMSL1() {return TMath::Sqrt(fMSL1);}
  void AddMST(float v) {fMST+=v;}
  void AddMSL0(float v) {fMSL0+=v;}
  void AddMSL1(float v) {fMSL1+=v;}
  void SetElectrons(float v) {SetDEnergy(v);}
  float GetElectrons() {return GetDEnergy();}
  virtual void CopyFrom(TPCHit *hit);

 protected:
  float fMST; // mean square in tranverse direction
  float fMSL0; // mean square in longudinal direction
  float fMSL1; // mean square in longudinal direction
};

#endif
