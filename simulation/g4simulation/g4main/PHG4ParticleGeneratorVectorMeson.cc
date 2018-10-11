#include "PHG4ParticleGeneratorVectorMeson.h"
#include "PHG4Particlev1.h"

#include "PHG4InEvent.h"

#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHRandomSeed.h>

#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4ParticleDefinition.hh>

#include <TLorentzVector.h>
#include <TF1.h>
#include <TRandom3.h>

#include <gsl/gsl_randist.h>

using namespace std;

PHG4ParticleGeneratorVectorMeson::PHG4ParticleGeneratorVectorMeson(const string &name): 
  PHG4ParticleGeneratorBase(name),
  decay1_names(),
  decay2_names(),
  decay_vtx_offset_x(),
  decay_vtx_offset_y(),
  decay_vtx_offset_z(),
  _vertex_func_x(Uniform),
  _vertex_func_y(Uniform),
  _vertex_func_z(Uniform),
  _t0(0.0),
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
  read_vtx_from_hepmc(true),
  y_min(0.),
  y_max(0.),
  eta_min(-1.0),
  eta_max(1.0),
  mom_min(0.0),
  mom_max(10.0),
  pt_min(4.),
  pt_max(4.),
  mass(9.46),
  width(54.02e-6),
  m1(0.511e-3),
  m2(0.511e-3),
  _histrand_init(0),
  decay1("e+"),
  decay2("e-"),
  fsin(NULL),
  frap(NULL),
  fpt(NULL),
  trand(NULL),
  ineve(NULL)
{

  // From PDG:
  // Upsilon 1S has mass 9.4603, width 54.02 keV
  // Upsilon 2S has mass 10.0233, width 31.98 keV
  // Upsilon 3S has mass 10.3552, width 20.32 keV
  return;
}

void 
PHG4ParticleGeneratorVectorMeson::add_decay_particles(const std::string &name1, const std::string &name2, const unsigned int decay_id)
{
  decay1_names.insert(std::pair<unsigned int, std::string>(decay_id, name1));
  decay2_names.insert(std::pair<unsigned int, std::string>(decay_id, name2));
  decay_vtx_offset_x.insert(std::pair<unsigned int, double>(decay_id, 0.));
  decay_vtx_offset_y.insert(std::pair<unsigned int, double>(decay_id, 0.));
  decay_vtx_offset_z.insert(std::pair<unsigned int, double>(decay_id, 0.));
  return;
}

void
PHG4ParticleGeneratorVectorMeson::set_decay_vertex_offset(double dx, double dy, double dz, const unsigned int decay_id)
{
  decay_vtx_offset_x.find(decay_id)->second = dx;
  decay_vtx_offset_y.find(decay_id)->second = dy;
  decay_vtx_offset_z.find(decay_id)->second = dz;
  return;
}

void
PHG4ParticleGeneratorVectorMeson::set_eta_range(const double min, const double max)
{
  eta_min = min;
  eta_max = max;
  return;
}


void
PHG4ParticleGeneratorVectorMeson::set_rapidity_range(const double min, const double max)
{
  y_min = min;
  y_max = max;
  return;
}


void
PHG4ParticleGeneratorVectorMeson::set_mom_range(const double min, const double max)
{
  mom_min = min;
  mom_max = max;
  return;
}

void
PHG4ParticleGeneratorVectorMeson::set_pt_range(const double min, const double max)
{
  pt_min = min;
  pt_max = max;
  return;
}

void 
PHG4ParticleGeneratorVectorMeson::set_vertex_distribution_function(FUNCTION x, FUNCTION y, FUNCTION z) {
  _vertex_func_x = x;
  _vertex_func_y = y;
  _vertex_func_z = z;
  return;
}


void 
PHG4ParticleGeneratorVectorMeson::set_vertex_distribution_mean(const double x, const double y, const double z) {
  _vertex_x = x;
  _vertex_y = y;
  _vertex_z = z;
  return;
}


void 
PHG4ParticleGeneratorVectorMeson::set_vertex_distribution_width(const double x, const double y, const double z) {
  _vertex_width_x = x;
  _vertex_width_y = y;
  _vertex_width_z = z;
  return;
}

void 
PHG4ParticleGeneratorVectorMeson::set_existing_vertex_offset_vector(const double x, const double y, const double z) {
  _vertex_offset_x = x;
  _vertex_offset_y = y;
  _vertex_offset_z = z;
  return;
}

void 
PHG4ParticleGeneratorVectorMeson::set_vertex_size_function(FUNCTION r) {
  _vertex_size_func_r = r;
  return;
}

void 
PHG4ParticleGeneratorVectorMeson::set_vertex_size_parameters(const double mean, const double width) {
  _vertex_size_mean = mean;
  _vertex_size_width = width;
  return;
}


void
PHG4ParticleGeneratorVectorMeson::set_mass(const double mass_in)
{
  mass = mass_in;
  return;
}

void
PHG4ParticleGeneratorVectorMeson::set_width(const double width_in)
{
  width = width_in;
  return;
}


void
PHG4ParticleGeneratorVectorMeson::set_decay_types(const std::string &name1, const std::string &name2)
{
  double mmuon = 105.64e-3;
  double melectron = 0.511e-3;

  decay1 = name1;
  if(decay1.compare("e+")==0 || decay1.compare("e-")==0)
    m1 = melectron;
  else if(decay1.compare("mu+")==0 || decay1.compare("mu-")==0)
    m1 = mmuon;
  else
    {
      cout << "Do not recognize the decay type " << decay1 << " will assume e+" << endl;
      decay1 = "e+";
      m1 = melectron;
    }


  decay2 = name2;
  if(decay2.compare("e+")==0 || decay2.compare("e-")==0)
    m2 = melectron;
  else if(decay2.compare("mu+")==0 || decay2.compare("mu-")==0)
    m2 = mmuon;
  else
    {
      cout << "Do not recognize the decay type " << decay2 << " will assume e-" << endl;
      decay2 = "e-";
      m2 = melectron;
    }

 return;
}


int
PHG4ParticleGeneratorVectorMeson::InitRun(PHCompositeNode *topNode)
{
  cout << "PHG4ParticleGeneratorVectorMeson::InitRun started." << endl;

  trand = new TRandom3();
  unsigned int iseed = PHRandomSeed(); // fixed seed handles in PHRandomSeed()
  cout << Name() << " random seed: " << iseed << endl;
  trand->SetSeed(iseed);
  if (_histrand_init)
    {
      iseed = PHRandomSeed();
      cout << Name() << " histrand random seed: " << iseed << endl;
      gRandom->SetSeed(iseed);
    }

  fsin = new TF1("fsin","sin(x)",0,M_PI);

  // From a fit to Pythia rapidity distribution for Upsilon(1S)
  frap = new TF1("frap","gaus(0)",y_min,y_max);
  frap->SetParameter(0,1.0);
  frap->SetParameter(1,0.0);
  frap->SetParameter(2,0.8749);

  // The dN/dPT  distribution is described by:
  fpt = new TF1("fpt","2.0*3.14159*x*[0]*pow((1 + x*x/(4*[1]) ),-[2])",pt_min,pt_max);
  fpt->SetParameter(0,72.1);
  fpt->SetParameter(1,26.516);
  fpt->SetParameter(2,10.6834);

  ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!ineve) {
    PHNodeIterator iter( topNode );
    PHCompositeNode *dstNode;
    dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

    ineve = new PHG4InEvent();
    PHDataNode<PHObject> *newNode = new PHDataNode<PHObject>(ineve, "PHG4INEVENT", "PHObject");
    dstNode->addNode(newNode);
  }


  cout << "PHG4ParticleGeneratorVectorMeson::InitRun endeded." << endl;
  return 0;
}



int
PHG4ParticleGeneratorVectorMeson::process_event(PHCompositeNode *topNode)
{
  if (!ineve) cout<<" G4InEvent not found "<<endl;

  // Generate a new set of vectors for the vector meson for each event
  // These are the momentum and direction vectors for the pre-decay vector meson

  // taken randomly from a fitted pT distribution to Pythia Upsilons

  double pt = 0.0;
  if(pt_max !=pt_min)
    {
      pt = fpt->GetRandom();
    }
  else
    {
      pt = pt_min;
    }
  // taken randomly from a fitted rapidity distribution to Pythia Upsilons

  double y = 0.0;
  if(y_max != y_min)
    {
      y  = frap->GetRandom();
    }
  else
    {
      y = y_min;
    }
  // 0 and 2*M_PI identical, so use gsl_rng_uniform which excludes 1.0
  double phi  = (2.0*M_PI)*gsl_rng_uniform(RandomGenerator);

  // The mass of the meson is taken from a Breit-Wigner lineshape

  double mnow = trand->BreitWigner(mass,width);

  // Get the pseudorapidity, eta, from the rapidity, mass and pt

  double mt = sqrt(pow(mnow,2) + pow(pt,2));
  double eta = TMath::ASinH(TMath::SinH(y)*mt/pt);

  // Put it in a TLorentzVector

  TLorentzVector vm;
  vm.SetPtEtaPhiM(pt,eta,phi,mnow);

  int vtxindex = -9;

  if (!ReuseExistingVertex(topNode))
    {
      // If not reusing existing vertex Randomly generate vertex position in z

      //                   mean   width
      vtx_x = smearvtx(_vertex_x,_vertex_width_x,_vertex_func_x);
      vtx_y = smearvtx(_vertex_y,_vertex_width_y,_vertex_func_y);
      vtx_z = smearvtx(_vertex_z,_vertex_width_z,_vertex_func_z);
    }
        vtx_x += _vertex_offset_x;
        vtx_y += _vertex_offset_y;
        vtx_z += _vertex_offset_z;



  for (std::map<unsigned int, std::string>::iterator it=decay1_names.begin(); it !=decay1_names.end() ; ++it){

    unsigned int decay_id = it->first;
    std::string decay1_name = it->second;
    std::string decay2_name;
    std::map<unsigned int, std::string>::iterator jt = decay2_names.find(decay_id);
    std::map<unsigned int, double>::iterator xt = decay_vtx_offset_x.find(decay_id);
    std::map<unsigned int, double>::iterator yt = decay_vtx_offset_y.find(decay_id);
    std::map<unsigned int, double>::iterator zt = decay_vtx_offset_z.find(decay_id);
    

    if (jt != decay2_names.end() && xt != decay_vtx_offset_x.end() && yt != decay_vtx_offset_y.end() && zt != decay_vtx_offset_z.end()){
    decay2_name = jt->second;
    set_decay_types(decay1_name, decay2_name);
    set_existing_vertex_offset_vector(xt->second, yt->second, zt->second);

    } else{
      cout << PHWHERE << "Decay particles && vertex info can't be found !!" << endl;
      exit(1);
    }

	// 3D Randomized vertex
      if ((_vertex_size_width > 0.0)||(_vertex_size_mean != 0.0)) {
	_vertex_size_mean = sqrt(pow(vtx_x,2)+pow(vtx_y,2)+pow(vtx_z,2));
        double r = smearvtx(_vertex_size_mean,_vertex_size_width,_vertex_size_func_r);
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        gsl_ran_dir_3d(RandomGenerator,&x,&y,&z);
        x *= r;
        y *= r;
        z *= r;
	vtxindex = ineve->AddVtx(vtx_x+x,vtx_y+y,vtx_z+z,t0);
      }
      else if (decay_id==0)
      {
	vtxindex = ineve->AddVtx(vtx_x,vtx_y,vtx_z,t0);
      }
	

    // Now decay it 
    // Get the decay energy and momentum in the frame of the vector meson - this correctly handles decay particles of any mass.

    double E1 = (pow(mnow,2) - pow(m2,2) + pow(m1,2)) / (2.0 * mnow);
    double p1 = sqrt( ( pow(mnow,2) - pow(m1+m2,2) )*( pow(mnow,2) - pow(m1-m2,2) ) ) / (2.0 * mnow);
  
    // In the frame of the vector meson, get a random theta and phi angle for particle 1
    // Assume angular distribution in the frame of the decaying meson that is uniform in phi and goes as sin(theta) in theta 
    // particle 2 has particle 1 momentum reflected through the origin

    double th1 = fsin->GetRandom();
    // 0 and 2*M_PI identical, so use gsl_rng_uniform which excludes 1.0
    double phi1 = 2.0*M_PI*gsl_rng_uniform(RandomGenerator);

    // Put particle 1 into a TLorentzVector

    double px1 = p1*sin(th1)*cos(phi1);
    double py1 = p1*sin(th1)*sin(phi1);
    double pz1 = p1*cos(th1);
    TLorentzVector v1;
    v1.SetPxPyPzE(px1,py1,pz1,E1);

    // now boost the decay product v1 into the lab using a vector consisting of the beta values of the vector meson
    // where p/E is v/c if we use GeV/c for p and GeV for E

    double betax = vm.Px()/vm.E();
    double betay = vm.Py()/vm.E();
    double betaz = vm.Pz()/vm.E();
    v1.Boost(betax,betay,betaz);

    // The second decay product's lab vector is the difference between the original meson and the boosted decay product 1

    TLorentzVector v2 = vm - v1;

    // Add the boosted decay particles to the particle list for the event

    AddParticle(decay1_name,v1.Px(),v1.Py(),v1.Pz());
    AddParticle(decay2_name,v2.Px(),v2.Py(),v2.Pz());

    // Now output the list of boosted decay particles to the node tree

    vector<PHG4Particle *>::const_iterator iter;
    for (iter = particlelist.begin(); iter != particlelist.end(); ++iter)
    {
      PHG4Particle *particle = new PHG4Particlev1(*iter);
      SetParticleId(particle,ineve);
      ineve->AddParticle(vtxindex, particle);
      if(embedflag!=0) { ineve->AddEmbeddedParticle(particle,embedflag); }
    }
    // List what has been put into ineve for this event

    if(Verbosity() > 0)
      {  
      ineve->identify();

      // Print some check output 
      cout << endl << "Output some sanity check info from PHG4ParticleGeneratorVectorMeson:" << endl;

      cout << "  Vertex for this event (X,Y,Z) is (" << vtx_x << ", " << vtx_y << ", " << vtx_z << ")" << endl;
      // Print the decay particle kinematics

      cout << "  Decay particle 1:"
	   << " px " << v1.Px()
	   << " py " << v1.Py()
	   << " pz " << v1.Pz()
	   << " eta " << v1.PseudoRapidity()
	   << " phi " << v1.Phi()
	   << " theta " << v1.Theta()
	   << " pT " << v1.Pt()
	   << " mass " << v1.M()
	   << " E " << v1.E()
	   << endl;

      cout << "  Decay particle 2:"
	   << " px " << v2.Px()
	   << " py " << v2.Py()
	   << " pz " << v2.Pz()
	   << " eta " << v2.PseudoRapidity()
	   << " phi " << v2.Phi()
	   << " theta " << v2.Theta()
	   << " pT " << v2.Pt()
	   << " mass " << v2.M()
	   << " E " << v2.E()
	   << endl; 

      // Print the input vector meson kinematics
      cout << "  Vector meson input kinematics:     mass " << vm.M()
	   << " px " << vm.Px()
  	   << " py " << vm.Py()
	   << " pz " << vm.Pz()
	   << " eta " << vm.PseudoRapidity()
	   << " y " << vm.Rapidity()
	   << " pt " << vm.Pt()
	   << " E " << vm.E()
	   << endl;

      // Now, as a check, reconstruct the mass from the particle 1 and 2 kinematics

      TLorentzVector vreco = v1 + v2;

      cout << "  Reco'd vector meson kinematics:    mass " << vreco.M()
	   << " px " << vreco.Px()
	   << " py " << vreco.Py()
	   << " pz " << vreco.Pz()
	   << " eta " << vreco.PseudoRapidity()
	   << " y " << vreco.Rapidity()
	   << " pt " << vreco.Pt()
	   << " E " << vreco.E()
	   << endl;


      }
  } // decay particles

  // Reset particlelist for the next event
  while(particlelist.begin() != particlelist.end())
    {
      delete particlelist.back();
      particlelist.pop_back();
    }

  return 0;
}

double
PHG4ParticleGeneratorVectorMeson::smearvtx(const double position, const double width, FUNCTION dist) const
{
  double res = position;
  if (dist == Uniform)
    {
      res = (position-width) + 2*gsl_rng_uniform_pos(RandomGenerator)*width;
    }
  else if (dist == Gaus)
    {
      res = position + gsl_ran_gaussian(RandomGenerator,width);
    }
  return res;
}

