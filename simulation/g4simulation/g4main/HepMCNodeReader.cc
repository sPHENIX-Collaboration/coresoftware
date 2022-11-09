#include "HepMCNodeReader.h"
#include "PHG4InEvent.h"
#include "PHG4Particle.h"
#include "PHG4Particlev1.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#pragma GCC diagnostic pop
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>
#include <HepMC/IteratorRange.h>
#include <HepMC/SimpleVector.h>
#include <HepMC/Units.h>

#include <CLHEP/Vector/LorentzRotation.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <list>
#include <utility>

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
  , is_pythia(false)
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
  string worldshape = rc->get_StringFlag("WorldShape");

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

    if (Verbosity())
    {
      cout << __PRETTY_FUNCTION__ << " : L" << __LINE__ << " Found PHHepMCGenEvent:" << endl;
      genevt->identify();
    }

    const auto collisionVertex = genevt->get_collision_vertex();

    HepMC::GenEvent *evt = genevt->getEvent();
    if (!evt)
    {
      cout << PHWHERE << " no evt pointer under HEPMC Node found:";
      genevt->identify();
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    if (Verbosity())
    {
      cout << __PRETTY_FUNCTION__ << " : L" << __LINE__ << " Found HepMC::GenEvent:" << endl;
      evt->print();
    }

    genevt->is_simulated(true);

    const int embed_flag = genevt->get_embedding_id();

    double xshift = vertex_pos_x + genevt->get_collision_vertex().x();
    double yshift = vertex_pos_y + genevt->get_collision_vertex().y();
    double zshift = vertex_pos_z + genevt->get_collision_vertex().z();
    double tshift = vertex_t0 + genevt->get_collision_vertex().t();

    const CLHEP::HepLorentzRotation lortentz_rotation(genevt->get_LorentzRotation_EvtGen2Lab());

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
      if (Verbosity() > 1)
      {
        cout << __PRETTY_FUNCTION__ << " : L" << __LINE__ << " Found vertex:" << endl;
        (*v)->print();
      }

      finalstateparticles.clear();
      for (HepMC::GenVertex::particle_iterator p =
               (*v)->particles_begin(HepMC::children);
           p != (*v)->particles_end(HepMC::children); ++p)
      {
        if (Verbosity() > 1)
        {
          cout << __PRETTY_FUNCTION__ << " : L" << __LINE__ << " Found particle:" << endl;
          (*p)->print();
          cout << "end vertex " << (*p)->end_vertex() << endl;
        }
        if (isfinal(*p))
        {
          if (Verbosity() > 1)
          {
            cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
            cout << "\tparticle passed " << endl;
          }
          finalstateparticles.push_back(*p);
        }
        else
        {
          if (Verbosity() > 1)
          {
            cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
            cout << "\tparticle failed " << endl;
          }
        }
      }

      if (!finalstateparticles.empty())
      {
        CLHEP::HepLorentzVector lv_vertex((*v)->position().x(),
                                          (*v)->position().y(),
                                          (*v)->position().z(),
                                          (*v)->position().t());
	if(is_pythia)
	  {
	    lv_vertex.setX(collisionVertex.x());
	    lv_vertex.setY(collisionVertex.y());
	    lv_vertex.setZ(collisionVertex.z());
	    lv_vertex.setT(collisionVertex.t());
	    if (Verbosity() > 1)
	      {
		std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ 
			  << std::endl;
		std::cout << "\t vertex reset to collision vertex: " 
			  << lv_vertex << std::endl;
	      }
	  }

        // event gen frame to lab frame
        lv_vertex = lortentz_rotation(lv_vertex);

        double xpos = lv_vertex.x() * length_factor + xshift;
        double ypos = lv_vertex.y() * length_factor + yshift;
        double zpos = lv_vertex.z() * length_factor + zshift;
        double time = lv_vertex.t() * time_factor + tshift;

        if (Verbosity() > 1)
        {
          cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
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
        for (fiter = finalstateparticles.begin();
             fiter != finalstateparticles.end();
             ++fiter)
        {
          if (Verbosity() > 1)
          {
            cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
            (*fiter)->print();
          }

          CLHEP::HepLorentzVector lv_momentum((*fiter)->momentum().px(),
                                              (*fiter)->momentum().py(),
                                              (*fiter)->momentum().pz(),
                                              (*fiter)->momentum().e());

          // event gen frame to lab frame
          lv_momentum = lortentz_rotation(lv_momentum);

          PHG4Particle *particle = new PHG4Particlev1();
          particle->set_pid((*fiter)->pdg_id());
          particle->set_px(lv_momentum.x() * mom_factor);
          particle->set_py(lv_momentum.y() * mom_factor);
          particle->set_pz(lv_momentum.z() * mom_factor);
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

void HepMCNodeReader::Embed(const int)
{
  cout << "HepMCNodeReader::Embed - WARNING - this function is depreciated. "
       << "Embedding IDs are controlled for individually HEPMC subevents in Fun4AllHepMCInputManagers and event generators."
       << endl;

  return;
}
