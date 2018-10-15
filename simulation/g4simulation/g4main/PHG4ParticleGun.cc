#include "PHG4ParticleGun.h"
#include "PHG4Particlev1.h"

#include "PHG4InEvent.h"

#include <phool/getClass.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>

#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4ParticleDefinition.hh>

using namespace std;

PHG4ParticleGun::PHG4ParticleGun(const string &name): 
  PHG4ParticleGeneratorBase(name)
{
  return;
}

PHG4ParticleGun::~PHG4ParticleGun()
{
  return;
}

int
PHG4ParticleGun::process_event(PHCompositeNode *topNode)
{
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode,"PHG4INEVENT");
  ReuseExistingVertex(topNode); // checks if we should reuse existing vertex
  int vtxindex = ineve->AddVtx(vtx_x,vtx_y,vtx_z,t0);
  vector<PHG4Particle *>::const_iterator iter;
  for (iter = particlelist.begin(); iter != particlelist.end(); ++iter)
    {
      PHG4Particle *particle = new PHG4Particlev1(*iter);
      SetParticleId(particle,ineve);
      ineve->AddParticle(vtxindex, particle);
    }
  if (Verbosity() > 0)
    {
      ineve->identify();
    }
  return 0;
}

