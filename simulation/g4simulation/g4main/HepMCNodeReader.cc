#include "HepMCNodeReader.h"
#include "PHG4InEvent.h"
#include "PHG4Particlev1.h"

#include <phhepmc/PHHepMCGenEvent.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>
#include <fun4all/recoConsts.h>

#include <HepMC/GenEvent.h>
#include <gsl/gsl_const.h>

#include <list>

using namespace std;
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
  SubsysReco(name)
{}

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
      if (finalstateparticles.size())
	{
	  double xpos = (*v)->position().x()*length_factor;
	  double ypos = (*v)->position().y()*length_factor;
	  double zpos = (*v)->position().z()*length_factor;
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
	  for (fiter = finalstateparticles.begin(); fiter != finalstateparticles.end(); fiter++)
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
	    }
	}
    }
  if (verbosity > 0)
    {
      ineve->identify();
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

