#ifndef __SIMPLETRACK3D__
#define __SIMPLETRACK3D__

#include <vector>
#include "SimpleHit3D.h"

class SimpleTrack3D {
 public:
  SimpleTrack3D(float phi0 = 0., float d0 = 0., float kappa0 = 0.,
                float dzdl0 = 0., float z00 = 0., unsigned int idx = 0)
      : phi(phi0), d(d0), kappa(kappa0), dzdl(dzdl0), z0(z00), index(idx) {}
  ~SimpleTrack3D() {}

  std::vector<SimpleHit3D> hits;
  float phi, d, kappa, dzdl, z0;
  unsigned int index;
};

#endif // __SIMPLETRACK3D__
