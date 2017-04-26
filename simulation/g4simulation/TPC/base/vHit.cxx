// Base HIT class
// Stores 3D point related SimHit attached to track in simulation
// Author: Carlos Perez
#include "vHit.h"

//==========
vHit::vHit():
 fTrack(-1) {
  for(int i=0; i!=4; ++i) fX[i] = 0.0;
}
//==========
vHit::~vHit() {
}
//==========
void vHit::CopyFrom(vHit *vh) {
  fTrack = vh->fTrack;
  for(int i=0; i!=4; ++i) fX[i] = vh->fX[i];
}
