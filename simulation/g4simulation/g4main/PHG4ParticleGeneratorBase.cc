#include "PHG4ParticleGeneratorBase.h"

#include "PHG4Particle.h"  // for PHG4Particle
#include "PHG4Particlev1.h"

#include "PHG4InEvent.h"
#include "PHG4TruthInfoContainer.h"
#include "PHG4VtxPoint.h"

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>      // for PHDataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>

#include <HepMC/SimpleVector.h>  // for FourVector

#include <gsl/gsl_rng.h>

#include <cassert>
#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, basic_ostream
#include <iterator>  // for operator!=, reverse_iterator
#include <map>       // for map<>::const_iterator, map
#include <utility>   // for pair

using namespace std;

PHG4ParticleGeneratorBase::PHG4ParticleGeneratorBase(const string &name)
  : SubsysReco(name)
  , reuse_existing_vertex(0)
  , vtx_x(0)
  , vtx_y(0)
  , vtx_z(0)
  , t0(0)
{
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  gsl_rng_set(m_RandomGenerator, seed);
  return;
}

PHG4ParticleGeneratorBase::~PHG4ParticleGeneratorBase()
{
  while (particlelist.begin() != particlelist.end())
  {
    delete particlelist.back();
    particlelist.pop_back();
  }
  gsl_rng_free(m_RandomGenerator);
  return;
}

int PHG4ParticleGeneratorBase::get_pdgcode(const std::string &name) const
{
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName = name;
  G4ParticleDefinition *particledef = particleTable->FindParticle(particleName);
  if (particledef)
  {
    return particledef->GetPDGEncoding();
  }
  return 0;
}

string
PHG4ParticleGeneratorBase::get_pdgname(const int pdgcode) const
{
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *particledef = particleTable->FindParticle(pdgcode);
  if (particledef)
  {
    return particledef->GetParticleName();
  }
  // if we cannot find the particle definition we'll make it ia geantino
  return "geantino";
}

double
PHG4ParticleGeneratorBase::get_mass(const int pdgcode) const
{
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *particledef = particleTable->FindParticle(get_pdgname(pdgcode));
  if (particledef)
  {
    return particledef->GetPDGMass() / GeV;
  }
  return 0;
}

void PHG4ParticleGeneratorBase::set_name(const std::string &particle)
{
  CheckAndCreateParticleVector();
  particlelist[0]->set_name(particle);
  particlelist[0]->set_pid(get_pdgcode(particle));
  return;
}

void PHG4ParticleGeneratorBase::set_pid(const int pid)
{
  CheckAndCreateParticleVector();
  particlelist[0]->set_pid(pid);
}

void PHG4ParticleGeneratorBase::set_mom(const double x, const double y, const double z)
{
  CheckAndCreateParticleVector();
  particlelist[0]->set_px(x);
  particlelist[0]->set_py(y);
  particlelist[0]->set_pz(z);
  return;
}

void PHG4ParticleGeneratorBase::set_vtx(const double x, const double y, const double z)
{
  vtx_x = x;
  vtx_y = y;
  vtx_z = z;
  return;
}

int PHG4ParticleGeneratorBase::InitRun(PHCompositeNode *topNode)
{
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!ineve)
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode;
    dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

    ineve = new PHG4InEvent();
    PHDataNode<PHObject> *newNode = new PHDataNode<PHObject>(ineve, "PHG4INEVENT", "PHObject");
    dstNode->addNode(newNode);
  }
  return 0;
}

int PHG4ParticleGeneratorBase::process_event(PHCompositeNode *topNode)
{
  cout << PHWHERE << " " << Name() << " using empty process_event" << endl;
  return 0;
}

void PHG4ParticleGeneratorBase::PrintParticles(const std::string &what) const
{
  vector<PHG4Particle *>::const_iterator iter;
  int i = 0;
  for (iter = particlelist.begin(); iter != particlelist.end(); ++iter)
  {
    cout << "particle " << i << endl;
    (*iter)->identify();
    i++;
  }
}

void PHG4ParticleGeneratorBase::AddParticle(const std::string &particle, const double x, const double y, const double z)
{
  PHG4Particle *part = new PHG4Particlev1(particle, get_pdgcode(particle), x, y, z);
  particlelist.push_back(part);
}

void PHG4ParticleGeneratorBase::AddParticle(const int pid, const double x, const double y, const double z)
{
  PHG4Particle *particle = new PHG4Particlev1();
  particle->set_pid(pid);
  particle->set_px(x);
  particle->set_py(y);
  particle->set_pz(z);
  particlelist.push_back(particle);
}

void PHG4ParticleGeneratorBase::CheckAndCreateParticleVector()
{
  if (!particlelist.size())
  {
    PHG4Particle *part = new PHG4Particlev1();
    particlelist.push_back(part);
  }
  return;
}

void PHG4ParticleGeneratorBase::SetParticleId(PHG4Particle *particle, PHG4InEvent *ineve)
{
  if ((particle->get_name()).size() == 0)  // no size -> empty name string
  {
    particle->set_name(get_pdgname(particle->get_pid()));
  }
  if (particle->get_pid() == 0)
  {
    particle->set_pid(get_pdgcode(particle->get_name()));
  }
  if (m_EmbedFlag)
  {
    ineve->AddEmbeddedParticle(particle, m_EmbedFlag);
  }
  return;
}

void PHG4ParticleGeneratorBase::set_seed(const unsigned int iseed)
{
  seed = iseed;
  cout << Name() << " random seed: " << seed << endl;
  gsl_rng_set(m_RandomGenerator, seed);
}

int PHG4ParticleGeneratorBase::ReuseExistingVertex(PHCompositeNode *topNode)
{
  if (!reuse_existing_vertex)
  {
    return 0;
  }

  PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

  if (genevtmap)
  {
    // use the highest priority HepMC subevent's vertex, ie that with highest embedding ID

    PHHepMCGenEventMap::ConstReverseIter iter =
        genevtmap->rbegin();

    if (iter != genevtmap->rend())
    {
      const PHHepMCGenEvent *hepmc_evt = iter->second;

      assert(hepmc_evt);

      const HepMC::FourVector &vtx = hepmc_evt->get_collision_vertex();

      set_vtx(vtx.x(), vtx.y(), vtx.z());

      if (Verbosity() > 0)
      {
        cout << "PHG4ParticleGeneratorBase::ReuseExistingVertex - reuse PHHepMCGenEventMap vertex "
             << vtx.x() << ", " << vtx.y() << ", " << vtx.z() << " cm. Source event:"
             << endl;
        hepmc_evt->identify();
      }

      return 1;
    }
  }

  // try PHG4INEVENT
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");

  if (ineve->GetNVtx() > 0)
  {
    std::pair<std::map<int, PHG4VtxPoint *>::const_iterator,
              std::map<int, PHG4VtxPoint *>::const_iterator>
        range = ineve->GetVertices();
    std::map<int, PHG4VtxPoint *>::const_iterator iter = range.first;
    PHG4VtxPoint *vtx = iter->second;

    if (!vtx)
    {
      cout << PHWHERE << "::Error - PHG4SimpleEventGenerator expects an existing vertex in PHG4InEvent, but none exists" << endl;
      exit(1);
    }
    if (Verbosity() > 0)
    {
      cout << PHWHERE << "::Info - use this primary vertex from PHG4InEvent:" << endl;
      vtx->identify();
    }

    set_vtx(vtx->get_x(), vtx->get_y(), vtx->get_z());
    return 1;

  }  // if (_ineve->GetNVtx() > 0) {
  // try the next option for getting primary vertex

  PHG4TruthInfoContainer *truthInfoList =  //
      findNode::getClass<PHG4TruthInfoContainer>(topNode,
                                                 "G4TruthInfo");
  if (truthInfoList)
  {
    // embed to vertexes as set by primary vertex from the truth container, e.g. when embedding into DST

    PHG4VtxPoint *vtx = truthInfoList->GetPrimaryVtx(1);

    if (!vtx)
    {
      truthInfoList->identify();
      cout << PHWHERE << "::Error - PHG4SimpleEventGenerator expects an existing truth vertex in PHG4TruthInfoContainer, but none exists"
           << endl;
      exit(1);
    }

    if (Verbosity() > 0)
    {
      cout << PHWHERE << "::Info - use this primary vertex from PHG4TruthInfoContainer:" << endl;
      vtx->identify();
    }

    set_vtx(vtx->get_x(), vtx->get_y(), vtx->get_z());
    return 1;
  }

  // I am out of options.....

  cout << PHWHERE << "::Error - PHG4SimpleEventGenerator expects an existing truth vertex, but none exists"
       << endl;
  exit(1);

  return 0;
}
