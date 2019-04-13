#ifndef PHG4ParticleGun_H__
#define PHG4ParticleGun_H__

#include "PHG4ParticleGeneratorBase.h"

class PHG4Particle;

class PHG4ParticleGun : public PHG4ParticleGeneratorBase
{
 public:
  PHG4ParticleGun(const std::string &name = "PGUN");
  virtual ~PHG4ParticleGun();
  int process_event(PHCompositeNode *topNode);

 protected:
};

#endif
