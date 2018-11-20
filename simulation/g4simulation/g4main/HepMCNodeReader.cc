#include "HepMCNodeReader.h"
#include "PHG4InEvent.h"
#include "PHG4Particlev1.h"
#include "PHG4Particlev2.h"

#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/G4ParticleTable.hh>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <HepMC/GenEvent.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cassert>
#include <list>

using namespace std;

//// All length Units are in cm, no conversion to G4 internal units since
//// this is filled into our objects (PHG4VtxPoint and PHG4Particle)
//
//// pythia vtx time seems to be in mm/c
//const double mm_over_c_to_sec = 0.1 / GSL_CONST_CGS_SPEED_OF_LIGHT;
//// pythia vtx time seems to be in mm/c
//const double mm_over_c_to_nanosecond = mm_over_c_to_sec * 1e9;
/// \class  IsStateFinal

/// this predicate returns true if the input has no decay vertex
class IsStateFinal
{
 public:
  /// returns true if the GenParticle does not decay
  bool operator()(const HepMC::GenParticle *p)
  {
    if (!p->end_vertex() && p->status() == 1) return 1;
    return 0;
  }
};

static IsStateFinal isfinal;

HepMCNodeReader::HepMCNodeReader(const std::string &name)
  : SubsysReco(name)
  , use_seed(0)
  , seed(0)
  , vertex_pos_x(0.0)
  , vertex_pos_y(0.0)
  , vertex_pos_z(0.0)
  , vertex_t0(0.0)
  , width_vx(0.0)
  , width_vy(0.0)
  , width_vz(0.0)
{
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  return;
}

HepMCNodeReader::~HepMCNodeReader()
{
  gsl_rng_free(RandomGenerator);
}

int HepMCNodeReader::Init(PHCompositeNode *topNode)
{
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!ineve)
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode;
    dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

    ineve = new PHG4InEvent();
    PHDataNode<PHObject> *newNode =
        new PHDataNode<PHObject>(ineve, "PHG4INEVENT", "PHObject");
    dstNode->addNode(newNode);
  }
  unsigned int phseed = PHRandomSeed();  // fixed seed is handled in this funtcion, need to call it to preserve random numbder order even if we override it
  if (use_seed)
  {
    phseed = seed;
    cout << Name() << " override random seed: " << phseed << endl;
  }
  gsl_rng_set(RandomGenerator, phseed);
  return 0;
}

int HepMCNodeReader::process_event(PHCompositeNode *topNode)
{
  // For pile-up simulation: define GenEventMap
  PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

  if (!genevtmap)
  {
    static bool once = true;

    if (once and Verbosity())
    {
      once = false;

      cout << "HepMCNodeReader::process_event - No PHHepMCGenEventMap node. Do not perform HepMC->Geant4 input" << endl;
    }

    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!ineve)
  {
    cout << PHWHERE << "no PHG4INEVENT node" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  recoConsts *rc = recoConsts::instance();

  float worldsizex = rc->get_FloatFlag("WorldSizex");
  float worldsizey = rc->get_FloatFlag("WorldSizey");
  float worldsizez = rc->get_FloatFlag("WorldSizez");
  string worldshape = rc->get_CharFlag("WorldShape");

  enum
  {
    ShapeG4Tubs = 0,
    ShapeG4Box = 1
  };

  int ishape;
  if (worldshape == "G4Tubs")
  {
    ishape = ShapeG4Tubs;
  }
  else if (worldshape == "G4Box")
  {
    ishape = ShapeG4Box;
  }
  else
  {
    cout << PHWHERE << " unknown world shape " << worldshape << endl;
    exit(1);
  }

  // For pile-up simulation: loop over PHHepMC event map
  // insert highest embedding ID event first, whose vertex maybe resued in  PHG4ParticleGeneratorBase::ReuseExistingVertex()
  int vtxindex = -1;
  for (PHHepMCGenEventMap::ReverseIter iter = genevtmap->rbegin(); iter != genevtmap->rend(); ++iter)
  {
    PHHepMCGenEvent *genevt = iter->second;
    assert(genevt);

    if (genevt->is_simulated())
    {
      if (Verbosity())
      {
        cout << "HepMCNodeReader::process_event - this event is already simulated. Move on: ";
        genevt->identify();
      }

      continue;
    }

    HepMC::GenEvent *evt = genevt->getEvent();
    if (!evt)
    {
      cout << PHWHERE << " no evt pointer under HEPMC Node found:";
      genevt->identify();
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    genevt->is_simulated(true);

    const int embed_flag = genevt->get_embedding_id();

    double xshift = vertex_pos_x + genevt->get_collision_vertex().x();
    double yshift = vertex_pos_y + genevt->get_collision_vertex().y();
    double zshift = vertex_pos_z + genevt->get_collision_vertex().z();
    double tshift = vertex_t0 + genevt->get_collision_vertex().t();

    if (width_vx > 0.0)
      xshift += smeargauss(width_vx);
    else if (width_vx < 0.0)
      xshift += smearflat(width_vx);

    if (width_vy > 0.0)
      yshift += smeargauss(width_vy);
    else if (width_vy < 0.0)
      yshift += smearflat(width_vy);

    if (width_vz > 0.0)
      zshift += smeargauss(width_vz);
    else if (width_vz < 0.0)
      zshift += smearflat(width_vz);

    std::list<HepMC::GenParticle *> finalstateparticles;
    std::list<HepMC::GenParticle *>::const_iterator fiter;

    // units in G4 interface are GeV and CM as in PHENIX convention
    const double mom_factor = HepMC::Units::conversion_factor(evt->momentum_unit(), HepMC::Units::GEV);
    const double length_factor = HepMC::Units::conversion_factor(evt->length_unit(), HepMC::Units::CM);
    const double time_factor = HepMC::Units::conversion_factor(evt->length_unit(), HepMC::Units::CM) / GSL_CONST_CGS_SPEED_OF_LIGHT * 1e9;  // from length_unit()/c to ns

    for (HepMC::GenEvent::vertex_iterator v = evt->vertices_begin();
         v != evt->vertices_end();
         ++v)
    {
      finalstateparticles.clear();
      for (HepMC::GenVertex::particle_iterator p =
               (*v)->particles_begin(HepMC::children);
           p != (*v)->particles_end(HepMC::children); ++p)
      {
        if (isfinal(*p))
        {
          finalstateparticles.push_back(*p);
        }
      }

      if (!finalstateparticles.empty())
      {
        double xpos = (*v)->position().x() * length_factor + xshift;
        double ypos = (*v)->position().y() * length_factor + yshift;
        double zpos = (*v)->position().z() * length_factor + zshift;
        double time = (*v)->position().t() * time_factor + tshift;

        if (Verbosity() > 1)
        {
          cout << "Vertex : " << endl;
          (*v)->print();
          cout << "id: " << (*v)->barcode() << endl;
          cout << "x: " << xpos << endl;
          cout << "y: " << ypos << endl;
          cout << "z: " << zpos << endl;
          cout << "t: " << time << endl;
          cout << "Particles" << endl;
        }

        if (ishape == ShapeG4Tubs)
        {
          if (sqrt(xpos * xpos + ypos * ypos) > worldsizey / 2 ||
              fabs(zpos) > worldsizez / 2)
          {
            cout << "vertex x/y/z" << xpos << "/" << ypos << "/" << zpos
                 << " outside world volume radius/z (+-) " << worldsizex / 2
                 << "/" << worldsizez / 2 << ", dropping it and its particles"
                 << endl;
            continue;
          }
        }
        else if (ishape == ShapeG4Box)
        {
          if (fabs(xpos) > worldsizex / 2 || fabs(ypos) > worldsizey / 2 ||
              fabs(zpos) > worldsizez / 2)
          {
            cout << "Vertex x/y/z " << xpos << "/" << ypos << "/" << zpos
                 << " outside world volume x/y/z (+-) " << worldsizex / 2 << "/"
                 << worldsizey / 2 << "/" << worldsizez / 2
                 << ", dropping it and its particles" << endl;
            continue;
          }
        }
        else
        {
          cout << PHWHERE << " shape " << ishape << " not implemented. exiting"
               << endl;
          exit(1);
        }

        // For pile-up simulation: vertex position
        vtxindex = ineve->AddVtx(xpos, ypos, zpos, time);
        int trackid = -1;
        for (fiter = finalstateparticles.begin();
             fiter != finalstateparticles.end();
             ++fiter)
        {
          ++trackid;

          if (Verbosity() > 1) (*fiter)->print();

          PHG4Particle *particle = new PHG4Particlev1();
          particle->set_pid((*fiter)->pdg_id());
          particle->set_px((*fiter)->momentum().px() * mom_factor);
          particle->set_py((*fiter)->momentum().py() * mom_factor);
          particle->set_pz((*fiter)->momentum().pz() * mom_factor);
          particle->set_barcode((*fiter)->barcode());

          ineve->AddParticle(vtxindex, particle);

          if (embed_flag != 0) ineve->AddEmbeddedParticle(particle, embed_flag);
        }
      }  //      if (!finalstateparticles.empty())

    }  //    for (HepMC::GenEvent::vertex_iterator v = evt->vertices_begin();

  }  // For pile-up simulation: loop end for PHHepMC event map

  if (Verbosity() > 0) ineve->identify();

  return Fun4AllReturnCodes::EVENT_OK;
}

double HepMCNodeReader::smeargauss(const double width)
{
  if (width == 0) return 0;
  return gsl_ran_gaussian(RandomGenerator, width);
}

double HepMCNodeReader::smearflat(const double width)
{
  if (width == 0) return 0;
  return 2.0 * width * (gsl_rng_uniform_pos(RandomGenerator) - 0.5);
}

void HepMCNodeReader::VertexPosition(const double v_x, const double v_y,
                                     const double v_z)
{
  cout << "HepMCNodeReader::VertexPosition - WARNING - this function is depreciated. "
       << "HepMCNodeReader::VertexPosition() move all HEPMC subevents to a new vertex location. "
       << "This also leads to a different vertex is used for HepMC subevent in Geant4 than that recorded in the HepMCEvent Node."
       << "Recommendation: the vertex shifts are better controlled for individually HEPMC subevents in Fun4AllHepMCInputManagers and event generators."
       << endl;

  vertex_pos_x = v_x;
  vertex_pos_y = v_y;
  vertex_pos_z = v_z;
  return;
}

void HepMCNodeReader::SmearVertex(const double s_x, const double s_y,
                                  const double s_z)
{
  cout << "HepMCNodeReader::SmearVertex - WARNING - this function is depreciated. "
       << "HepMCNodeReader::SmearVertex() smear each HEPMC subevents to a new vertex location. "
       << "This also leads to a different vertex is used for HepMC subevent in Geant4 than that recorded in the HepMCEvent Node."
       << "Recommendation: the vertex smears are better controlled for individually HEPMC subevents in Fun4AllHepMCInputManagers and event generators."
       << endl;

  width_vx = s_x;
  width_vy = s_y;
  width_vz = s_z;
  return;
}

void HepMCNodeReader::Embed(const int i)
{
  cout << "HepMCNodeReader::Embed - WARNING - this function is depreciated. "
       << "Embedding IDs are controlled for individually HEPMC subevents in Fun4AllHepMCInputManagers and event generators."
       << endl;

  return;
}
