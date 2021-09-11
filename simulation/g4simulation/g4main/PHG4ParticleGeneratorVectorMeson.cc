#include "PHG4ParticleGeneratorVectorMeson.h"

#include "PHG4InEvent.h"
#include "PHG4Particle.h"  // for PHG4Particle
#include "PHG4Particlev1.h"

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>      // for PHDataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TF1.h>
#include <TLorentzVector.h>
#include <TRandom.h>  // for TRandom
#include <TRandom3.h>
#include <TSystem.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform, gsl_rng_uniform_pos

#include <cmath>     // for sin, sqrt, cos, M_PI
#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, basic_ostream, basic_...
#include <utility>   // for pair
#include <vector>    // for vector, vector<>::const_iterator

PHG4ParticleGeneratorVectorMeson::PHG4ParticleGeneratorVectorMeson(const std::string &name)
  : PHG4ParticleGeneratorBase(name)
{
  set_upsilon_1s();  // make mass and width of 1S default
  return;
}

PHG4ParticleGeneratorVectorMeson::
    ~PHG4ParticleGeneratorVectorMeson()
{
  delete trand;
}

void PHG4ParticleGeneratorVectorMeson::add_decay_particles(const std::string &name1, const unsigned int decay_id)
{
  if (name1.compare("e") == 0)
  {
    add_decay_particles("e+", "e-", decay_id);
  }
  else if (name1.compare("mu") == 0)
  {
    add_decay_particles("mu+", "mu-", decay_id);
  }
  else
  {
    std::cout << "Invalid decay " << name1 << ", valid is e or mu" << std::endl;
    gSystem->Exit(1);
  }
  return;
}

void PHG4ParticleGeneratorVectorMeson::add_decay_particles(const std::string &name1, const std::string &name2, const unsigned int decay_id)
{
  // check for valid select ion (e+,e- or mu+,mu-)
  if ((name1.compare("e-") == 0 && name2.compare("e+") == 0) ||
      (name1.compare("e+") == 0 && name2.compare("e-") == 0) ||
      (name1.compare("mu+") == 0 && name2.compare("mu-") == 0) ||
      (name1.compare("mu-") == 0 && name2.compare("mu+") == 0))
  {
    decay1_names.insert(std::pair<unsigned int, std::string>(decay_id, name1));
    decay2_names.insert(std::pair<unsigned int, std::string>(decay_id, name2));
    decay_vtx_offset_x.insert(std::pair<unsigned int, double>(decay_id, 0.));
    decay_vtx_offset_y.insert(std::pair<unsigned int, double>(decay_id, 0.));
    decay_vtx_offset_z.insert(std::pair<unsigned int, double>(decay_id, 0.));
    return;
  }
  std::cout << "invalid decay channel Y --> " << name1 << " + " << name2 << std::endl;
  gSystem->Exit(1);
}

void PHG4ParticleGeneratorVectorMeson::set_decay_vertex_offset(double dx, double dy, double dz, const unsigned int decay_id)
{
  decay_vtx_offset_x.find(decay_id)->second = dx;
  decay_vtx_offset_y.find(decay_id)->second = dy;
  decay_vtx_offset_z.find(decay_id)->second = dz;
  return;
}

void PHG4ParticleGeneratorVectorMeson::set_eta_range(const double min, const double max)
{
  eta_min = min;
  eta_max = max;
  return;
}

void PHG4ParticleGeneratorVectorMeson::set_rapidity_range(const double min, const double max)
{
  y_min = min;
  y_max = max;
  return;
}

void PHG4ParticleGeneratorVectorMeson::set_mom_range(const double min, const double max)
{
  mom_min = min;
  mom_max = max;
  return;
}

void PHG4ParticleGeneratorVectorMeson::set_pt_range(const double min, const double max)
{
  pt_min = min;
  pt_max = max;
  return;
}

void PHG4ParticleGeneratorVectorMeson::set_vertex_distribution_function(FUNCTION x, FUNCTION y, FUNCTION z)
{
  _vertex_func_x = x;
  _vertex_func_y = y;
  _vertex_func_z = z;
  return;
}

void PHG4ParticleGeneratorVectorMeson::set_vertex_distribution_mean(const double x, const double y, const double z)
{
  _vertex_x = x;
  _vertex_y = y;
  _vertex_z = z;
  return;
}

void PHG4ParticleGeneratorVectorMeson::set_vertex_distribution_width(const double x, const double y, const double z)
{
  _vertex_width_x = x;
  _vertex_width_y = y;
  _vertex_width_z = z;
  return;
}

void PHG4ParticleGeneratorVectorMeson::set_existing_vertex_offset_vector(const double x, const double y, const double z)
{
  _vertex_offset_x = x;
  _vertex_offset_y = y;
  _vertex_offset_z = z;
  return;
}

void PHG4ParticleGeneratorVectorMeson::set_vertex_size_function(FUNCTION r)
{
  _vertex_size_func_r = r;
  return;
}

void PHG4ParticleGeneratorVectorMeson::set_vertex_size_parameters(const double mean, const double width)
{
  _vertex_size_mean = mean;
  _vertex_size_width = width;
  return;
}

void PHG4ParticleGeneratorVectorMeson::set_decay_types(const std::string &name1, const std::string &name2)
{
  //http://pdg.lbl.gov/2020/listings/rpp2020-list-muon.pdf
  static const double mmuon = 105.6583745e-3;       //+-0.0000024e-3
                                                    // http://pdg.lbl.gov/2020/listings/rpp2020-list-electron.pdf
  static const double melectron = 0.5109989461e-3;  //+-0.0000000031e-3

  decay1 = name1;
  if (decay1.compare("e+") == 0 || decay1.compare("e-") == 0)
  {
    m1 = melectron;
  }
  else if (decay1.compare("mu+") == 0 || decay1.compare("mu-") == 0)
  {
    m1 = mmuon;
  }
  else
  {
    std::cout << "Do not recognize the decay type " << decay1 << std::endl;
    gSystem->Exit(1);
  }

  decay2 = name2;
  if (decay2.compare("e+") == 0 || decay2.compare("e-") == 0)
  {
    m2 = melectron;
  }
  else if (decay2.compare("mu+") == 0 || decay2.compare("mu-") == 0)
  {
    m2 = mmuon;
  }
  else
  {
    std::cout << "Do not recognize the decay type " << decay2 << std::endl;
    gSystem->Exit(1);
  }

  return;
}

int PHG4ParticleGeneratorVectorMeson::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4ParticleGeneratorVectorMeson::InitRun started." << std::endl;
  }
  trand = new TRandom3();
  unsigned int iseed = PHRandomSeed();  // fixed seed handles in PHRandomSeed()
  trand->SetSeed(iseed);
  if (_histrand_init)
  {
    iseed = PHRandomSeed();
    gRandom->SetSeed(iseed);
  }

  fsin = new TF1("fsin", "sin(x)", 0, M_PI);

  // From a fit to Pythia rapidity distribution for Upsilon(1S)
  frap = new TF1("frap", "gaus(0)", y_min, y_max);
  frap->SetParameter(0, 1.0);
  frap->SetParameter(1, 0.0);
  frap->SetParameter(2, 0.8749);

  // The dN/dPT  distribution is described by:
  fpt = new TF1("fpt", "2.0*3.14159*x*[0]*pow((1 + x*x/(4*[1]) ),-[2])", pt_min, pt_max);
  fpt->SetParameter(0, 72.1);
  fpt->SetParameter(1, 26.516);
  fpt->SetParameter(2, 10.6834);

  ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!ineve)
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode;
    dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

    ineve = new PHG4InEvent();
    PHDataNode<PHObject> *newNode = new PHDataNode<PHObject>(ineve, "PHG4INEVENT", "PHObject");
    dstNode->addNode(newNode);
  }

  if (Verbosity() > 0)
  {
    std::cout << "PHG4ParticleGeneratorVectorMeson::InitRun endeded." << std::endl;
  }
  return 0;
}

int PHG4ParticleGeneratorVectorMeson::process_event(PHCompositeNode *topNode)
{
  if (!ineve) std::cout << " G4InEvent not found " << std::endl;

  // Generate a new set of vectors for the vector meson for each event
  // These are the momentum and direction vectors for the pre-decay vector meson

  // taken randomly from a fitted pT distribution to Pythia Upsilons

  double pt = 0.0;
  if (pt_max != pt_min)
  {
    pt = fpt->GetRandom();
  }
  else
  {
    pt = pt_min;
  }
  // taken randomly from a fitted rapidity distribution to Pythia Upsilons

  double y = 0.0;
  if (y_max != y_min)
  {
    y = frap->GetRandom();
  }
  else
  {
    y = y_min;
  }
  // 0 and 2*M_PI identical, so use gsl_rng_uniform which excludes 1.0
  double phi = (2.0 * M_PI) * gsl_rng_uniform(RandomGenerator());

  // The mass of the meson is taken from a Breit-Wigner lineshape

  double mnow = trand->BreitWigner(mass, m_Width);

  // Get the pseudorapidity, eta, from the rapidity, mass and pt

  double mt = sqrt(mnow * mnow + pt * pt);
  double eta = asinh(sinh(y) * mt / pt);

  // Put it in a TLorentzVector

  TLorentzVector vm;
  vm.SetPtEtaPhiM(pt, eta, phi, mnow);

  int vtxindex = -9;

  if (!ReuseExistingVertex(topNode))
  {
    // If not reusing existing vertex Randomly generate vertex position in z

    //                   mean   width
    set_vtx(smearvtx(_vertex_x, _vertex_width_x, _vertex_func_x),
            smearvtx(_vertex_y, _vertex_width_y, _vertex_func_y),
            smearvtx(_vertex_z, _vertex_width_z, _vertex_func_z));
  }
  set_vtx(get_vtx_x() + _vertex_offset_x,
          get_vtx_y() + _vertex_offset_y,
          get_vtx_z() + _vertex_offset_z);

  for (std::map<unsigned int, std::string>::iterator it = decay1_names.begin(); it != decay1_names.end(); ++it)
  {
    unsigned int decay_id = it->first;
    std::string decay1_name = it->second;
    std::string decay2_name;
    std::map<unsigned int, std::string>::iterator jt = decay2_names.find(decay_id);
    std::map<unsigned int, double>::iterator xt = decay_vtx_offset_x.find(decay_id);
    std::map<unsigned int, double>::iterator yt = decay_vtx_offset_y.find(decay_id);
    std::map<unsigned int, double>::iterator zt = decay_vtx_offset_z.find(decay_id);

    if (jt != decay2_names.end() && xt != decay_vtx_offset_x.end() && yt != decay_vtx_offset_y.end() && zt != decay_vtx_offset_z.end())
    {
      decay2_name = jt->second;
      set_decay_types(decay1_name, decay2_name);
      set_existing_vertex_offset_vector(xt->second, yt->second, zt->second);
    }
    else
    {
      std::cout << PHWHERE << "Decay particles && vertex info can't be found !!" << std::endl;
      exit(1);
    }

    // 3D Randomized vertex
    if ((_vertex_size_width > 0.0) || (_vertex_size_mean != 0.0))
    {
      _vertex_size_mean = sqrt(get_vtx_x() * get_vtx_x() +
                               get_vtx_y() * get_vtx_y() +
                               get_vtx_z() * get_vtx_z());
      double r = smearvtx(_vertex_size_mean, _vertex_size_width, _vertex_size_func_r);
      double x1 = 0.0;
      double y1 = 0.0;
      double z1 = 0.0;
      gsl_ran_dir_3d(RandomGenerator(), &x1, &y1, &z1);
      x1 *= r;
      y1 *= r;
      z1 *= r;
      vtxindex = ineve->AddVtx(get_vtx_x() + x1, get_vtx_y() + y1, get_vtx_z() + z1, get_t0());
    }
    else if (decay_id == 0)
    {
      vtxindex = ineve->AddVtx(get_vtx_x(), get_vtx_y(), get_vtx_z(), get_t0());
    }

    // Now decay it
    // Get the decay energy and momentum in the frame of the vector meson - this correctly handles decay particles of any mass.

    double E1 = (mnow * mnow - m2 * m2 + m1 * m1) / (2.0 * mnow);
    double p1 = sqrt((mnow * mnow - (m1 + m2) * (m1 + m2)) * (mnow * mnow - (m1 - m2) * (m1 - m2))) / (2.0 * mnow);

    // In the frame of the vector meson, get a random theta and phi angle for particle 1
    // Assume angular distribution in the frame of the decaying meson that is uniform in phi and goes as sin(theta) in theta
    // particle 2 has particle 1 momentum reflected through the origin

    double th1 = fsin->GetRandom();
    // 0 and 2*M_PI identical, so use gsl_rng_uniform which excludes 1.0
    double phi1 = 2.0 * M_PI * gsl_rng_uniform(RandomGenerator());

    // Put particle 1 into a TLorentzVector

    double px1 = p1 * sin(th1) * cos(phi1);
    double py1 = p1 * sin(th1) * sin(phi1);
    double pz1 = p1 * cos(th1);
    TLorentzVector v1;
    v1.SetPxPyPzE(px1, py1, pz1, E1);

    // now boost the decay product v1 into the lab using a vector consisting of the beta values of the vector meson
    // where p/E is v/c if we use GeV/c for p and GeV for E

    double betax = vm.Px() / vm.E();
    double betay = vm.Py() / vm.E();
    double betaz = vm.Pz() / vm.E();
    v1.Boost(betax, betay, betaz);

    // The second decay product's lab vector is the difference between the original meson and the boosted decay product 1

    TLorentzVector v2 = vm - v1;

    // Add the boosted decay particles to the particle list for the event

    AddParticle(decay1_name, v1.Px(), v1.Py(), v1.Pz());
    AddParticle(decay2_name, v2.Px(), v2.Py(), v2.Pz());

    // Now output the list of boosted decay particles to the node tree

    for (std::vector<PHG4Particle *>::const_iterator iter = particlelist_begin(); iter != particlelist_end(); ++iter)
    {
      PHG4Particle *particle = new PHG4Particlev1(*iter);
      SetParticleId(particle, ineve);
      ineve->AddParticle(vtxindex, particle);
      if (EmbedFlag() != 0)
      {
        ineve->AddEmbeddedParticle(particle, EmbedFlag());
      }
    }
    // List what has been put into ineve for this event

    if (Verbosity() > 0)
    {
      ineve->identify();

      // Print some check output
      std::cout << std::endl
                << "Output some sanity check info from PHG4ParticleGeneratorVectorMeson:" << std::endl;

      std::cout << "  Vertex for this event (X,Y,Z) is (" << get_vtx_x() << ", " << get_vtx_y() << ", " << get_vtx_z() << ")" << std::endl;
      // Print the decay particle kinematics

      std::cout << "  Decay particle 1:"
                << " px " << v1.Px()
                << " py " << v1.Py()
                << " pz " << v1.Pz()
                << " eta " << v1.PseudoRapidity()
                << " phi " << v1.Phi()
                << " theta " << v1.Theta()
                << " pT " << v1.Pt()
                << " mass " << v1.M()
                << " E " << v1.E()
                << std::endl;

      std::cout << "  Decay particle 2:"
                << " px " << v2.Px()
                << " py " << v2.Py()
                << " pz " << v2.Pz()
                << " eta " << v2.PseudoRapidity()
                << " phi " << v2.Phi()
                << " theta " << v2.Theta()
                << " pT " << v2.Pt()
                << " mass " << v2.M()
                << " E " << v2.E()
                << std::endl;

      // Print the input vector meson kinematics
      std::cout << "  Vector meson input kinematics:     mass " << vm.M()
                << " px " << vm.Px()
                << " py " << vm.Py()
                << " pz " << vm.Pz()
                << " eta " << vm.PseudoRapidity()
                << " y " << vm.Rapidity()
                << " pt " << vm.Pt()
                << " E " << vm.E()
                << std::endl;

      // Now, as a check, reconstruct the mass from the particle 1 and 2 kinematics

      TLorentzVector vreco = v1 + v2;

      std::cout << "  Reco'd vector meson kinematics:    mass " << vreco.M()
                << " px " << vreco.Px()
                << " py " << vreco.Py()
                << " pz " << vreco.Pz()
                << " eta " << vreco.PseudoRapidity()
                << " y " << vreco.Rapidity()
                << " pt " << vreco.Pt()
                << " E " << vreco.E()
                << std::endl;
    }
  }  // decay particles

  ResetParticleList();

  return 0;
}

double
PHG4ParticleGeneratorVectorMeson::smearvtx(const double position, const double /*width*/, FUNCTION dist) const
{
  double res = position;
  if (dist == Uniform)
  {
    res = (position - m_Width) + 2 * gsl_rng_uniform_pos(RandomGenerator()) * m_Width;
  }
  else if (dist == Gaus)
  {
    res = position + gsl_ran_gaussian(RandomGenerator(), m_Width);
  }
  return res;
}

void PHG4ParticleGeneratorVectorMeson::set_upsilon_1s()
{
  // http://pdg.lbl.gov/2020/listings/rpp2020-list-upsilon-1S.pdf
  set_mass(9.4603);     //+- 0.00026
  set_width(54.02e-6);  // +- 1.25e-6
}

void PHG4ParticleGeneratorVectorMeson::set_upsilon_2s()
{
  // http://pdg.lbl.gov/2020/listings/rpp2020-list-upsilon-2S.pdf
  set_mass(10.02326);   // +- 0.00031
  set_width(31.98e-6);  // +- 2.63e-6
}

void PHG4ParticleGeneratorVectorMeson::set_upsilon_3s()
{
  //http://pdg.lbl.gov/2020/listings/rpp2020-list-upsilon-3S.pdf
  set_mass(10.3552);    // +- 0.0005
  set_width(20.32e-6);  // +- 1.85e-6
}
