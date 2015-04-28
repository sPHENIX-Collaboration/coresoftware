#include "PHG4SimpleEventGenerator.h"

#include "PHG4Particlev2.h"
#include "PHG4InEvent.h"
#include "PHG4VtxPoint.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>

#include <TRandom3.h>

#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4ParticleDefinition.hh>

#include <cstdlib>
#include <cmath>

using namespace std;

PHG4SimpleEventGenerator::PHG4SimpleEventGenerator(const string &name): 
  SubsysReco(name),
  _seed(0),
  _rand(NULL),
  _particle_codes(),
  _particle_names(),
  _reuse_existing_vertex(false),
  _vertex_func_x(Uniform),
  _vertex_func_y(Uniform),
  _vertex_func_z(Uniform),
  _vertex_x(0.0),
  _vertex_y(0.0),
  _vertex_z(0.0),
  _vertex_width_x(0.0),
  _vertex_width_y(0.0),
  _vertex_width_z(0.0),
  _vertex_offset_x(0.0),
  _vertex_offset_y(0.0),
  _vertex_offset_z(0.0),
  _vertex_size_func_r(Uniform),
  _vertex_size_mean(0.0),
  _vertex_size_width(0.0),
  _eta_min(-1.25),
  _eta_max(1.25),
  _phi_min(-M_PI),
  _phi_max(M_PI),
  _pt_min(0.0),
  _pt_max(10.0),
  _embedflag(0),
  _ineve(NULL) {
  return;
}

PHG4SimpleEventGenerator::~PHG4SimpleEventGenerator() {
  delete _rand;
  return;
}

void PHG4SimpleEventGenerator::set_seed(int seed) {
  _seed = seed;
  _rand = new TRandom3(seed);
  return;
}

void PHG4SimpleEventGenerator::add_particles(std::string name, unsigned int num) {
  _particle_names.push_back(std::make_pair(name,num));
  return;
}

void PHG4SimpleEventGenerator::add_particles(int pid, unsigned int num) {
  _particle_codes.push_back(std::make_pair(pid,num));
  return;
}

void PHG4SimpleEventGenerator::set_eta_range(double min, double max) {
  _eta_min = min;
  _eta_max = max;
  return;
}

void PHG4SimpleEventGenerator::set_phi_range(double min, double max) {
  _phi_min = min;
  _phi_max = max;
  return;
}

void PHG4SimpleEventGenerator::set_pt_range(double min, double max) {
  _pt_min = min;
  _pt_max = max;
  return;
}

void PHG4SimpleEventGenerator::set_vertex_distribution_function(FUNCTION x, FUNCTION y, FUNCTION z) {
  _vertex_func_x = x;
  _vertex_func_y = y;
  _vertex_func_z = z;
  return;
}

void PHG4SimpleEventGenerator::set_vertex_distribution_mean(double x, double y, double z) {
  _vertex_x = x;
  _vertex_y = y;
  _vertex_z = z;
  return;
}

void PHG4SimpleEventGenerator::set_vertex_distribution_width(double x, double y, double z) {
  _vertex_width_x = x;
  _vertex_width_y = y;
  _vertex_width_z = z;
  return;
}

void PHG4SimpleEventGenerator::set_existing_vertex_offset_vector(double x, double y, double z) {
  _vertex_offset_x = x;
  _vertex_offset_y = y;
  _vertex_offset_z = z;
  return;
}

void PHG4SimpleEventGenerator::set_vertex_size_function(FUNCTION r) {
  _vertex_size_func_r = r;
  return;
}

void PHG4SimpleEventGenerator::set_vertex_size_parameters(double mean, double width) {
  _vertex_size_mean = mean;
  _vertex_size_width = width;
  return;
}

int PHG4SimpleEventGenerator::InitRun(PHCompositeNode *topNode) {

  if ((_vertex_func_x != Uniform)&&(_vertex_func_x != Gaus)) {
    cout << PHWHERE << "::Error - unknown vertex distribution function requested" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if ((_vertex_func_y != Uniform)&&(_vertex_func_y != Gaus)) {
    cout << PHWHERE << "::Error - unknown vertex distribution function requested" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if ((_vertex_func_z != Uniform)&&(_vertex_func_z != Gaus)) {
    cout << PHWHERE << "::Error - unknown vertex distribution function requested" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!_ineve) {
    PHNodeIterator iter( topNode );
    PHCompositeNode *dstNode;
    dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
      
    _ineve = new PHG4InEvent();
    PHDataNode<PHObject> *newNode = new PHDataNode<PHObject>(_ineve, "PHG4INEVENT", "PHObject");
    dstNode->addNode(newNode);
  }

  if (verbosity >= 0) {
    cout << "================ PHG4SimpleEventGenerator::InitRun() ======================" << endl;
    cout << " Random seed = " << _seed << endl;
    cout << " Particles:" << endl;
    for (unsigned int i=0; i<_particle_codes.size(); ++i) {
      cout << "    " << _particle_codes[i].first << ", count = " << _particle_codes[i].second << endl;
    }
    for (unsigned int i=0; i<_particle_names.size(); ++i) {
      cout << "    " << _particle_names[i].first << ", count = " << _particle_names[i].second << endl;
    }
    if (_reuse_existing_vertex) {
      cout << " Vertex Distribution: Set to reuse a previously generated sim vertex" << endl;
      cout << " Vertex offset vector (x,y,z) = (" << _vertex_offset_x << ","<< _vertex_offset_y << ","<< _vertex_offset_z << ")" << endl;
    } else {
      cout << " Vertex Distribution Function (x,y,z) = ("; 
      if (_vertex_func_x == Uniform) cout << "Uniform,";
      else if (_vertex_func_x == Gaus) cout << "Gaus,";
      if (_vertex_func_y == Uniform) cout << "Uniform,";
      else if (_vertex_func_y == Gaus) cout << "Gaus,";
      if (_vertex_func_z == Uniform) cout << "Uniform";
      else if (_vertex_func_z == Gaus) cout << "Gaus";
      cout << ")" << endl;
      cout << " Vertex mean (x,y,z) = (" << _vertex_x << ","<< _vertex_y << ","<< _vertex_z << ")" << endl;
      cout << " Vertex width (x,y,z) = (" << _vertex_width_x << ","<< _vertex_width_y << ","<< _vertex_width_z << ")" << endl;
    }
    cout << " Vertex size function (r) = (";
    if (_vertex_size_func_r == Uniform) cout << "Uniform";
    else if (_vertex_size_func_r == Gaus) cout << "Gaus";
    cout << ")" << endl;
    cout << " Vertex size (mean) = (" << _vertex_size_mean << ")" << endl;
    cout << " Vertex size (width) = (" << _vertex_size_width << ")" << endl;
    cout << " Eta range = " << _eta_min << " - " << _eta_max << endl;
    cout << " Phi range = " << _phi_min << " - " << _phi_max << endl;
    cout << " pT range = " << _pt_min << " - " << _pt_max << endl;
    cout << "===========================================================================" << endl;
  }

  // the definition table should be filled now, so convert codes into names
  for (unsigned int i=0;i<_particle_codes.size();++i) {
    int pdgcode = _particle_codes[i].first;
    unsigned int count = _particle_codes[i].second;
    string pdgname = get_pdgname(pdgcode);
    _particle_names.push_back(make_pair(pdgname,count));
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SimpleEventGenerator::process_event(PHCompositeNode *topNode) {

  double vertex_x = 0.0;
  double vertex_y = 0.0;
  double vertex_z = 0.0;

  if (!_reuse_existing_vertex) {
    // generate a new vertex point

    if (_vertex_func_x == Uniform) {
      vertex_x = _rand->Uniform(_vertex_x-_vertex_width_x,_vertex_x+_vertex_width_x);
    } else if (_vertex_func_x == Gaus) {
      vertex_x = _rand->Gaus(_vertex_x,_vertex_width_x);
    }

    if (_vertex_func_y == Uniform) {
      vertex_y = _rand->Uniform(_vertex_y-_vertex_width_y,_vertex_y+_vertex_width_y);
    } else if (_vertex_func_y == Gaus) {
      vertex_y = _rand->Gaus(_vertex_y,_vertex_width_y);
    }

    if (_vertex_func_z == Uniform) {
      vertex_z = _rand->Uniform(_vertex_z-_vertex_width_z,_vertex_z+_vertex_width_z);
    } else if (_vertex_func_z == Gaus) {
      vertex_z = _rand->Gaus(_vertex_z,_vertex_width_z);
    }
  } else {
    
    if (_ineve->GetNVtx() == 0) {
      cout << PHWHERE << "::Error - PHG4SimpleEventGenerator expects an existing truth vertex, but none exists" << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    // use the first vertex in the list
    std::pair< std::map<int, PHG4VtxPoint *>::const_iterator, 
	       std::map<int, PHG4VtxPoint *>::const_iterator > 
      range = _ineve->GetVertices();
    std::map<int, PHG4VtxPoint* >::const_iterator iter = range.first;
    PHG4VtxPoint* vtx = iter->second;

    if (!vtx) {
      cout << PHWHERE << "::Error - PHG4SimpleEventGenerator expects an existing truth vertex, but none exists" << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    vertex_x = vtx->get_x();
    vertex_y = vtx->get_y();
    vertex_z = vtx->get_z();

    vertex_x += _vertex_offset_x;
    vertex_y += _vertex_offset_y;
    vertex_z += _vertex_offset_z;
  }

  int vtxindex = -1;
  int trackid = -1;
  for (unsigned int i=0; i<_particle_names.size(); ++i) {

    string pdgname = _particle_names[i].first;
    int pdgcode = get_pdgcode(pdgname);
    unsigned int nparticles = _particle_names[i].second;
    
    for (unsigned int j=0; j<nparticles; ++j) {

      if ((_vertex_size_width > 0.0)||(_vertex_size_mean != 0.0)) {

	double r = 0.0;
	if (_vertex_size_func_r == Uniform) {
	  r = _rand->Uniform(_vertex_size_mean-_vertex_size_width,_vertex_size_mean+_vertex_size_width);
	} else if (_vertex_size_func_r == Gaus) {
	  r = fabs(_rand->Gaus(_vertex_size_mean,_vertex_size_width));
	}

	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	_rand->Sphere(x,y,z,r);

	vtxindex = _ineve->AddVtx(vertex_x+x,vertex_y+y,vertex_z+z,0.0);
      } else if ((i==0)&&(j==0)) {
	vtxindex = _ineve->AddVtx(vertex_x,vertex_y,vertex_z,0.0);
      }

      ++trackid;
   
      double eta = _rand->Uniform(_eta_min,_eta_max);
      double phi = _rand->Uniform(_phi_min,_phi_max);
      double pt  = _rand->Uniform(_pt_min,_pt_max);   
      double px = pt*cos(phi);
      double py = pt*sin(phi);
      double pz = pt*sinh(eta);
      double m = get_mass(pdgcode);
      double e = sqrt(px*px+py*py+pz*pz+m*m);

      PHG4Particle *particle = new PHG4Particlev2();
      particle->set_track_id(trackid);
      particle->set_vtx_id(vtxindex);
      particle->set_parent_id(0);
      particle->set_name(pdgname);
      particle->set_pid(pdgcode);
      particle->set_px(px);
      particle->set_py(py);
      particle->set_pz(pz);
      particle->set_e(e);

      _ineve->AddParticle(vtxindex, particle);
      if (_embedflag != 0) _ineve->AddEmbeddedParticle(particle);
    }
  }

  if (verbosity > 0) {
    cout << "====================== PHG4SimpleEventGenerator::process_event() =====================" << endl;
    _ineve->identify();
    cout << "======================================================================================" << endl;
  } 

  return Fun4AllReturnCodes::EVENT_OK;
}

// call only during execution!
int PHG4SimpleEventGenerator::get_pdgcode(std::string name) {
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName = name;
  G4ParticleDefinition* particledef = particleTable->FindParticle(particleName);
  if (!particledef) return 0;
  return particledef->GetPDGEncoding();
}

// call only during execution!
std::string PHG4SimpleEventGenerator::get_pdgname(int pdgcode) {
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particledef = particleTable->FindParticle(pdgcode);
  if (!particledef) return 0;
  return particledef->GetParticleName();
}

// call only during execution!
double PHG4SimpleEventGenerator::get_mass(int pdgcode) {
  if (pdgcode == 0) return 0.0;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle_definition = particleTable->FindParticle(get_pdgname(pdgcode));
  if (!particle_definition) return 0.0;
  return particle_definition->GetPDGMass();
}
