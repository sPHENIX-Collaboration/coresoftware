#include "PHG4ParticleGun.h"

#include "PHG4InEvent.h"
#include "PHG4Particlev1.h"

#include <phool/getClass.h>

#include <vector>  // for vector, vector<>::const_iterator

class PHCompositeNode;
class PHG4Particle;

PHG4ParticleGun::PHG4ParticleGun(const std::string &name)
  : PHG4ParticleGeneratorBase(name)
{
  return;
}

int PHG4ParticleGun::process_event(PHCompositeNode *topNode)
{
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  ReuseExistingVertex(topNode);  // checks if we should reuse existing vertex
  int vtxindex = ineve->AddVtx(get_vtx_x(), get_vtx_y(), get_vtx_z(), get_t0());
  std::vector<PHG4Particle *>::const_iterator iter;
  for (iter = particlelist_begin(); iter != particlelist_end(); ++iter)
  {
    PHG4Particle *particle = new PHG4Particlev1(*iter);
    SetParticleId(particle, ineve);
    ineve->AddParticle(vtxindex, particle);
  }
  if (Verbosity() > 0)
  {
    ineve->identify();
  }
  return 0;
}
