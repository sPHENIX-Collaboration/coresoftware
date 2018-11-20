#include "PHG4ParticleGeneratorD0.h"
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
#include <gsl/gsl_const.h>

using namespace std;

PHG4ParticleGeneratorD0::PHG4ParticleGeneratorD0(const string &name): 
  PHG4ParticleGeneratorBase(name),
  vtx_zmin(-10.),
  vtx_zmax(10),
  y_min(0.),
  y_max(0.),
  eta_min(-1.0),
  eta_max(1.0),
  mom_min(0.0),
  mom_max(10.0),
  pt_min(4.),
  pt_max(4.),
  mass(1.86486),
  m1(0.493677),
  m2(0.13957018),
  fsin(NULL),
  frap(NULL),
  fpt(NULL)
{
  return;
}

void
PHG4ParticleGeneratorD0::set_eta_range(const double min, const double max)
{
  eta_min = min;
  eta_max = max;
  return;
}


void
PHG4ParticleGeneratorD0::set_rapidity_range(const double min, const double max)
{
  y_min = min;
  y_max = max;
  return;
}


void
PHG4ParticleGeneratorD0::set_mom_range(const double min, const double max)
{
  mom_min = min;
  mom_max = max;
  return;
}

void
PHG4ParticleGeneratorD0::set_pt_range(const double min, const double max)
{
  pt_min = min;
  pt_max = max;
  return;
}

void
PHG4ParticleGeneratorD0::set_vtx_zrange(const double zmin, const double zmax)
{
  vtx_zmin = zmin;
  vtx_zmax = zmax;

  return;
}

void
PHG4ParticleGeneratorD0::set_mass(const double mass_in)
{
  mass = mass_in;
  return;
}

int
PHG4ParticleGeneratorD0::InitRun(PHCompositeNode *topNode)
{
  unsigned int iseed = PHRandomSeed(); // fixed seed handled in PHRandomSeed()
  cout << Name() << " random seed: " << iseed << endl;
  gRandom->SetSeed(iseed);

  fsin = new TF1("fsin","sin(x)",0,M_PI);

  // From a fit to Pythia rapidity distribution for Upsilon(1S)
  frap = new TF1("frap","gaus(0)",y_min,y_max);
  frap->SetParameter(0,1.0);
  frap->SetParameter(1,0.0);
  frap->SetParameter(2,0.8749);


  // The dN/dPT  distribution is described by:
  fpt = new TF1("fpt","[0]*x*pow((1.0 + x*x/[1]/[1]),[2])",pt_min,pt_max);
  fpt->SetParameter(0,1.16930e+04);
  fpt->SetParameter(1,3.03486e+00);
  fpt->SetParameter(2,-5.42819e+00);


  return 0;
}

int
PHG4ParticleGeneratorD0::process_event(PHCompositeNode *topNode)
{
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode,"PHG4INEVENT");

  double tau = 410.1e-15; // D0 life time in seconds
  double cc = GSL_CONST_CGSM_SPEED_OF_LIGHT; // speed of light cm/sec
  double ctau = tau*cc*1.0e+04;  // ctau in micrometers

  // If not reusing existing vertex Randomly generate vertex position in z 
  if (! ReuseExistingVertex(topNode))
    {

      if (vtx_zmax != vtx_zmin)
	{
	  vtx_z = (vtx_zmax - vtx_zmin) * gsl_rng_uniform_pos(RandomGenerator) + vtx_zmin;
	}
      else
	{
	  vtx_z = vtx_zmin;
	}
    }
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

  double mnow = mass;

  // Get the pseudorapidity, eta, from the rapidity, mass and pt

  double mt = sqrt(pow(mnow,2) + pow(pt,2));
  double eta = TMath::ASinH(TMath::SinH(y)*mt/pt);

  // Put it in a TLorentzVector

  TLorentzVector vd0;
  vd0.SetPtEtaPhiM(pt,eta,phi,mnow);

  double beta = vd0.Beta();
  double gamma = vd0.Gamma();
  double lifetime = gsl_ran_exponential(RandomGenerator,410.1e-03)* 1.0e-12; // life time in seconds
  double lifepath = lifetime*gamma*beta*cc;  // path in cm
  if (Verbosity() > 0)
    {
      cout << "D0 px,py,pz: " << vd0.Px() << " " << vd0.Py() << " " << vd0.Pz() << " " << beta << " " << gamma << endl;
      cout << "   ctau = " << ctau << " " << lifetime << " " << lifepath << " " << lifepath*1.0e+04 << endl;
    }  

  vtx_x = vd0.Px()/vd0.P()*lifepath;
  vtx_y = vd0.Py()/vd0.P()*lifepath;
  vtx_z = vtx_z + vd0.Pz()/vd0.P()*lifepath;
  t0 = lifetime;
  int vtxindex = ineve->AddVtx(vtx_x,vtx_y,vtx_z,t0);
  if (Verbosity() > 0)
    {
      cout << "  XY vertex: " << sqrt(vtx_x*vtx_x+vtx_y*vtx_y) << " " << sqrt(vtx_x*vtx_x+vtx_y*vtx_y)*1.0e+04 << endl;
    }

  // Now decay it 
  // Get the decay energy and momentum in the frame of the D0 - this correctly handles decay particles of any mass.

  double E1 = (pow(mnow,2) - pow(m2,2) + pow(m1,2)) / (2.0 * mnow);
  double p1 = sqrt( ( pow(mnow,2) - pow(m1+m2,2) )*( pow(mnow,2) - pow(m1-m2,2) ) ) / (2.0 * mnow);
  
  // In the frame of the D0, get a random theta and phi angle for particle 1
  // Assume angular distribution in the frame of the decaying D0 that is uniform in phi and goes as sin(theta) in theta 
  // particle 2 has particle 1 momentum reflected through the origin

  double th1 = fsin->GetRandom();
  double phi1 = 2.0*M_PI*gsl_rng_uniform_pos(RandomGenerator);

  // Put particle 1 into a TLorentzVector

  double px1 = p1*sin(th1)*cos(phi1);
  double py1 = p1*sin(th1)*sin(phi1);
  double pz1 = p1*cos(th1);
  TLorentzVector v1;
  v1.SetPxPyPzE(px1,py1,pz1,E1);

  // now boost the decay product v1 into the lab using a vector consisting of the beta values of the vector meson
  // where p/E is v/c if we use GeV/c for p and GeV for E

  double betax = vd0.Px()/vd0.E();
  double betay = vd0.Py()/vd0.E();
  double betaz = vd0.Pz()/vd0.E();
  v1.Boost(betax,betay,betaz);

  // The second decay product's lab vector is the difference between the original D0 and the boosted decay product 1

  TLorentzVector v2 = vd0 - v1;

  // Add the boosted decay particles to the particle list for the event

  AddParticle(-321,v1.Px(),v1.Py(),v1.Pz()); // K-
  AddParticle(211,v2.Px(),v2.Py(),v2.Pz());  // pi+

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
      cout << endl << "Output some sanity check info from PHG4ParticleGeneratorD0:" << endl;

      cout << "Event vertex: (" << vtx_x << ", " << vtx_y << ", " << vtx_z << ")" << endl;
      
      cout << "Kaon : " << v1.Pt() << " " << v1.PseudoRapidity() << " " << v1.M() << endl;
      cout << "pion : " << v2.Pt() << " " << v2.PseudoRapidity() << " " << v2.M() << endl;
      cout << "D0 : " << vd0.Pt() << " " << vd0.PseudoRapidity() << " " << vd0.M() << endl;
      TLorentzVector vreco = v1 + v2;
      cout << "reconstructed D0 : " << vreco.Pt() << " " << vreco.PseudoRapidity() << " " << vreco.M() << endl;

/*
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
      cout << " D0 input kinematics:     mass " << vd0.M()
	   << " px " << vd0.Px()
  	   << " py " << vd0.Py()
	   << " pz " << vd0.Pz()
	   << " eta " << vd0.PseudoRapidity()
	   << " y " << vd0.Rapidity()
	   << " pt " << vd0.Pt()
	   << " E " << vd0.E()
	   << endl;

      // Now, as a check, reconstruct the mass from the particle 1 and 2 kinematics

      TLorentzVector vreco = v1 + v2;

      cout << "  Reconstructed D0 kinematics:    mass " << vreco.M()
	   << " px " << vreco.Px()
	   << " py " << vreco.Py()
	   << " pz " << vreco.Pz()
	   << " eta " << vreco.PseudoRapidity()
	   << " y " << vreco.Rapidity()
	   << " pt " << vreco.Pt()
	   << " E " << vreco.E()
	   << endl;
*/

    }

  // Reset particlelist for the next event
  while(particlelist.begin() != particlelist.end())
    {
      delete particlelist.back();
      particlelist.pop_back();
    }

  return 0;
}

