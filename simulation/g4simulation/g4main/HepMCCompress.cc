#include "HepMCCompress.h"
#include "PHG4InEvent.h"
#include "PHG4Particle.h"
#include <vararray/VariableArrayContainer.h>
#include <vararray/VariableArray.h>
#include <vararray/VariableArrayIds.h>
#include <half/half.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <HepMC/GenEvent.h>
#include <gsl/gsl_const.h>

#include <list>

using namespace std;
const double mm_over_c_to_sec = 0.1 / GSL_CONST_CGS_SPEED_OF_LIGHT; // pythia vtx time seems to be in mm/c
/// \class  IsStateFinal

/// this predicate returns true if the input has no decay vertex
class IsStateFinal
{
public:
  /// returns true if the GenParticle does not decay
  bool operator()( const HepMC::GenParticle* p )
  {
    if ( !p->end_vertex() && p->status() == 1 ) return 1;
    return 0;
  }
};

static IsStateFinal isfinal;

HepMCCompress::HepMCCompress(const std::string &name):
  SubsysReco(name)
{}

int
HepMCCompress::Init(PHCompositeNode *topNode)
{
  VariableArrayContainer *vararraycontainer = findNode::getClass<VariableArrayContainer>(topNode, "HEPMC_VarArray");
  if (!vararraycontainer)
    {
      PHNodeIterator iter( topNode );
      PHCompositeNode *dstNode;
      dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

      vararraycontainer = new VariableArrayContainer();
      VariableArray *vararray = new VariableArray(varids::G4VTXV1);
      vararraycontainer->AddVarArray(vararray);
      vararray = new VariableArray(varids::G4PARTICLEV1);
      vararraycontainer->AddVarArray(vararray);

      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(vararraycontainer, "HEPMC_VarArray", "PHObject");
      dstNode->addNode(newNode);
    }
  return 0;
}

int
HepMCCompress::process_event(PHCompositeNode *topNode)
{
  HepMC::GenEvent *evt = findNode::getClass<HepMC::GenEvent>(topNode, "HEPMC");
  if (!evt)
    {
      cout << PHWHERE << " no evt pointer under HEPMC Node found" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  VariableArrayContainer *vararraycontainer = findNode::getClass<VariableArrayContainer>(topNode, "HEPMC_VarArray");
  if (!vararraycontainer)
    {
      cout << PHWHERE << "no PHG4INEVENT node" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  vector<short> shepmcvtxvec;

  std::list<HepMC::GenParticle*> finalstateparticles;
  // units in G4 interface are GeV and CM
  //  const double mom_factor = HepMC::Units::conversion_factor( evt->momentum_unit(), HepMC::Units::GEV );
  const double length_factor = HepMC::Units::conversion_factor( evt->length_unit(), HepMC::Units::CM );
  for ( HepMC::GenEvent::vertex_iterator v = evt->vertices_begin();
        v != evt->vertices_end(); ++v )
    {
      finalstateparticles.clear();
      for (HepMC::GenVertex::particle_iterator p = (*v)->particles_begin(HepMC::children); p != (*v)->particles_end(HepMC::children); ++p)
        {
          if (isfinal(*p))
            {
              if (!select_pid.empty())
                {
                  if (select_pid.find((*p)->pdg_id()) != select_pid.end())
                    {
                      finalstateparticles.push_back(*p);
                    }
                  continue;
                }
              if (!exclude_pid.empty())
                {
                  if (exclude_pid.find((*p)->pdg_id()) != exclude_pid.end())
                    {
                      continue;
                    }
                }
              finalstateparticles.push_back(*p);
            }
        }
      if (!finalstateparticles.empty())
        {
          //  	  cout << "Vertex : " << endl;
          //  	    (*v)->print();
          //  	    cout << "id: " << (*v)->barcode() << endl;
          //  	    cout << "x: " << (*v)->position().x() << endl;
          //  	    cout << "y: " << (*v)->position().y() << endl;
          //  	    cout << "z: " << (*v)->position().z() << endl;
          // 	    cout << "t: " << (*v)->position().t() << endl;
          // 	  cout << "Particles" << endl;
          shepmcvtxvec.push_back((*v)->barcode());
          shepmcvtxvec.push_back(FloatToInt((*v)->position().x()*length_factor));
          shepmcvtxvec.push_back(FloatToInt((*v)->position().y()*length_factor));
          shepmcvtxvec.push_back(FloatToInt((*v)->position().z()*length_factor));

          // 	  for (fiter = finalstateparticles.begin(); fiter != finalstateparticles.end(); fiter++)
          // 	    {
          // 	      //	      (*fiter)->print();
          // 	      PHG4Particle *particle = new PHG4Particle();
          // 	      particle->set_pid((*fiter)->pdg_id());
          // 	      particle->set_px((*fiter)->momentum().px()*mom_factor);
          // 	      particle->set_py((*fiter)->momentum().py()*mom_factor);
          // 	      particle->set_pz((*fiter)->momentum().pz()*mom_factor);
          // 	      ineve->AddParticle((*v)->barcode(), particle);
          // 	    }
          // 	}
        }
    }
  //  ineve->identify();
  return Fun4AllReturnCodes::EVENT_OK;
}

short int
HepMCCompress::FloatToInt(const float rval) const
{
  half ftoi(rval);
  return ftoi.bits();
}
