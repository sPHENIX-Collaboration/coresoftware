// Base HIT class
// Stores 3D point related SimHit attached to track in simulation
// Author: Carlos Perez
#ifndef __vHIT_H__
#define __vHIT_H__

#include <phool/PHObject.h>

class vHit : public PHObject {
 public:
  vHit();
  virtual ~vHit();
  Int_t GetTrack() const {return fTrack;}
  void SetTrack(Int_t track) {fTrack=track;}

  Float_t X(Int_t i) const {return fX[i];}
  void SetX(Int_t i, Float_t v) {fX[i]=v;}
  void AddX(Int_t i, Float_t d) {fX[i]+=d;}

  virtual void CopyFrom(vHit *vh);

 protected:
  Int_t fTrack; // geant track id
  Float_t fX[4]; // X0 X1 X2 L

  ClassDef(vHit,1);
};

#endif
