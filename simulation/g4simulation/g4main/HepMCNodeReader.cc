#include "HepMCNodeReader.h"
#include "PHG4InEvent.h"
#include "PHG4Particlev1.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phhepmc/PHHepMCGenEvent.h>

#include <phool/getClass.h>
#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <HepMC/GenEvent.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <list>

using namespace std;

// All length Units are in cm, no conversion to G4 internal units since
// this is filled into our objects (PHG4VtxPoint and PHG4Particle)

const double mm_over_c_to_sec = 0.1/GSL_CONST_CGS_SPEED_OF_LIGHT; // pythia vtx time seems to be in mm/c
/// \class  IsStateFinal

/// this predicate returns true if the input has no decay vertex
class IsStateFinal {
public:
    /// returns true if the GenParticle does not decay
    bool operator()( const HepMC::GenParticle* p ) { 
	if ( !p->end_vertex() && p->status()==1 ) return 1;
	return 0;
    }
};

static IsStateFinal isfinal;

HepMCNodeReader::HepMCNodeReader(const std::string &name):
  SubsysReco(name),
  _embed_flag(0),
  vertex_pos_x(0),
  vertex_pos_y(0),
  vertex_pos_z(0),
  width_vx(0),
  width_vy(0),
  width_vz(0)
{
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned int seed = PHRandomSeed(); // fixed seed is handled in this funtcion
  gsl_rng_set(RandomGenerator,seed);
  return;
}

HepMCNodeReader::~HepMCNodeReader()
{
  gsl_rng_free (RandomGenerator);
}

int
HepMCNodeReader::Init(PHCompositeNode *topNode)
{
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode,"PHG4INEVENT");
  if (!ineve)
    {
      PHNodeIterator iter( topNode );
      PHCompositeNode *dstNode;
      dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

      ineve = new PHG4InEvent();
      PHDataNode<PHObject> *newNode = new PHDataNode<PHObject>(ineve, "PHG4INEVENT", "PHObject");
      dstNode->addNode(newNode);
    }
  return 0;
}

int
HepMCNodeReader::process_event(PHCompositeNode *topNode)
{
  recoConsts *rc = recoConsts::instance();
  float worldsizex = rc->get_FloatFlag("WorldSizex");
  float worldsizey = rc->get_FloatFlag("WorldSizey");
  float worldsizez = rc->get_FloatFlag("WorldSizez");
  string worldshape = rc->get_CharFlag("WorldShape");
  enum {ShapeG4Tubs = 0, ShapeG4Box = 1};
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
  double xshift = vertex_pos_x; + smeargauss(width_vx);
  double yshift = vertex_pos_y; + smeargauss(width_vy);
  double zshift = vertex_pos_z; + smeargauss(width_vz);

  if (width_vx > 0.0) xshift += smeargauss(width_vx);
  else                xshift += smearflat(width_vx);

  if (width_vy > 0.0) yshift += smeargauss(width_vy);
  else                yshift += smearflat(width_vy);
  
  if (width_vz > 0.0) zshift += smeargauss(width_vz);
  else                zshift += smearflat(width_vz);
  
  PHHepMCGenEvent *genevt = findNode::getClass<PHHepMCGenEvent>(topNode,"PHHepMCGenEvent");
  
  HepMC::GenEvent *evt = genevt->getEvent();
  if (!evt)
    {
      cout << PHWHERE << " no evt pointer under HEPMC Node found" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode,"PHG4INEVENT");
  if (!ineve)
    {
      cout << PHWHERE << "no PHG4INEVENT node" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  std::list<HepMC::GenParticle*> finalstateparticles;
  std::list<HepMC::GenParticle*>::const_iterator  fiter;
// units in G4 interface are GeV and CM
  const double mom_factor = HepMC::Units::conversion_factor( evt->momentum_unit(), (HepMC::Units::MomentumUnit) (genevt->get_momentumunit())); 
  const double length_factor = HepMC::Units::conversion_factor( evt->length_unit(), (HepMC::Units::LengthUnit) (genevt->get_lengthunit()) );
  for ( HepMC::GenEvent::vertex_iterator v = evt->vertices_begin();
	v != evt->vertices_end(); ++v )
    {
      finalstateparticles.clear();
      for (HepMC::GenVertex::particle_iterator p = (*v)->particles_begin(HepMC::children); p != (*v)->particles_end(HepMC::children); ++p)
	{
	  if (isfinal(*p))
	    {
	      finalstateparticles.push_back(*p);
	    }
	}
      if (!finalstateparticles.empty())
	{
	  double xpos = (*v)->position().x()*length_factor + xshift;
	  double ypos = (*v)->position().y()*length_factor + yshift;
	  double zpos = (*v)->position().z()*length_factor + zshift;
 	  if (verbosity > 1)
	    {
	      cout << "Vertex : " << endl;
	      (*v)->print();
	      cout << "id: " << (*v)->barcode() << endl;
	      cout << "x: " << xpos << endl;
	      cout << "y: " << ypos << endl;
	      cout << "z: " << zpos << endl;
	      cout << "t: " << (*v)->position().t()*mm_over_c_to_sec << endl;
	      cout << "Particles" << endl;
	    }

      if (ishape == ShapeG4Tubs)
        {
          if (sqrt(xpos*xpos + ypos*ypos) > worldsizey / 2 || fabs(zpos) > worldsizez / 2)
            {
              cout << "vertex x/y/z" << xpos << "/" << ypos << "/" << zpos
                   << " outside world volume radius/z (+-) " << worldsizex / 2 << "/"
                   << worldsizez / 2
                   << ", dropping it and its particles" << endl;
              continue;
            }
        }
      else if (ishape == ShapeG4Box)
        {
	  if (fabs(xpos) > worldsizex/2 || fabs(ypos) > worldsizey/2 || fabs(zpos) > worldsizez/2)
	    {
	      cout << "Vertex x/y/z " << xpos << "/" << ypos << "/" << zpos
                   << " outside world volume x/y/z (+-) " << worldsizex/2 << "/"
                   << worldsizey/2 << "/" << worldsizez/2
                   << ", dropping it and its particles" << endl;
	      continue;
	    }
        }
      else
        {
          cout << PHWHERE << " shape " << ishape << " not implemented. exiting" << endl;
          exit(1);
        }

	  ineve->AddVtxHepMC((*v)->barcode(), xpos, ypos, zpos, (*v)->position().t()*mm_over_c_to_sec);
	  for (fiter = finalstateparticles.begin(); fiter != finalstateparticles.end(); ++fiter)
	    {
	      if (verbosity > 1)
		{
		  (*fiter)->print();
		}
	      PHG4Particle *particle = new PHG4Particlev1();
	      particle->set_pid((*fiter)->pdg_id());
	      particle->set_px((*fiter)->momentum().px()*mom_factor);
	      particle->set_py((*fiter)->momentum().py()*mom_factor);
	      particle->set_pz((*fiter)->momentum().pz()*mom_factor);
	      ineve->AddParticle((*v)->barcode(), particle);
	      if (_embed_flag != 0) ineve->AddEmbeddedParticle(particle,_embed_flag);
	    }
	}
    }
  if (verbosity > 0)
    {
      ineve->identify();
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

double
HepMCNodeReader::smeargauss(const double width)
{
  if (width == 0)
    {
      return 0;
    }
  return gsl_ran_gaussian(RandomGenerator,width);
}

double
HepMCNodeReader::smearflat(const double width)
{
  if (width == 0)
    {
      return 0;
    }
  return 2.0*width*(gsl_rng_uniform_pos(RandomGenerator) - 0.5);
}

void
HepMCNodeReader::VertexPosition(const double v_x, const double v_y, const double v_z)
{
  vertex_pos_x = v_x;
  vertex_pos_y = v_y;
  vertex_pos_z = v_z;
  return;
}

void
HepMCNodeReader::SmearVertex(const double s_x, const double s_y, const double s_z)
{
  width_vx = s_x;
  width_vy = s_y;
  width_vz = s_z;
  return;
}
