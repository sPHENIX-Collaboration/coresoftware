// Base CLUSTER class
// Stores 3D point related reconstructed cluster for central arm
// Origin: Carlos Perez
#ifndef __vCLUSTER_H__
#define __vCLUSTER_H__

#include <phool/PHObject.h>

class vCluster : public PHObject {
 public:
  vCluster();
  virtual ~vCluster();
  Float_t GetX(Int_t i) const {return fX[i];}
  Float_t GetCovariance(int i) {return fCov[i];}
  Float_t GetSignal() {return fSgn;}
  UChar_t GetSize(int i) {return fSize[i];}

  void SetX(Int_t i, Float_t v) {fX[i]=v;}
  void SetCovariance(Int_t i, Float_t v) {fCov[i]=v;}
  void SetSignal(Float_t v) {fSgn=v;}
  void SetSize(Int_t i, UChar_t v) {fSize[i]=v;}

  virtual void CopyFrom(vCluster *vc);

 protected:
  Float_t fSgn; // signal associated
  Float_t fX[3]; // X0 X1 X2
  Float_t fCov[6]; // X0X0 X1X1 X2X2 X0X1 X0X2 X1X2
  UChar_t fSize[3]; // S0 S1 S2
};

#endif
