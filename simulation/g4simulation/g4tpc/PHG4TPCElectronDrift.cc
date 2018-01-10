#include "PHG4TPCElectronDrift.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <TSystem.h>

#include <Geant4/G4SystemOfUnits.hh>

#include <gsl/gsl_randist.h>

#include <iostream>

using namespace std;

PHG4TPCElectronDrift::PHG4TPCElectronDrift(const std::string& name):
  SubsysReco(name),
  PHG4ParameterInterface(name),
  diffusion_trans(NAN),
  diffusion_long(NAN),
  drift_velocity(NAN),
  electrons_per_gev(NAN)
{  
  InitializeParameters();
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  set_seed(PHRandomSeed()); // fixed seed is handled in this funtcion
  return;
}

int PHG4TPCElectronDrift::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      exit(1);
    }
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));
//  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "PAR" ));
  string paramnodename = "G4CELLPARAM_" + detector;
  hitnodename = "G4HIT_" + detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
    {
      cout << Name() << " Could not locate G4HIT node " << hitnodename << endl;
      topNode->print();
      gSystem->Exit(1);
      exit(1);
    }
   UpdateParametersWithMacro();
  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode =  dynamic_cast<PHCompositeNode*>(runIter.findFirst("PHCompositeNode",detector));
  if (! RunDetNode)
    {
      RunDetNode = new PHCompositeNode(detector);
      runNode->addNode(RunDetNode);
    }
  SaveToNodeTree(RunDetNode,paramnodename);
  diffusion_long = get_double_param("diffusion_long");
  diffusion_trans = get_double_param("diffusion_trans");
  drift_velocity = get_double_param("drift_velocity");
  electrons_per_gev = get_double_param("electrons_per_gev");
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TPCElectronDrift::process_event(PHCompositeNode *topNode)
{
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
    {
      cout << "Could not locate g4 hit node " << hitnodename << endl;
      gSystem->Exit(1);
    }
  PHG4HitContainer::ConstIterator hiter;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  double tpc_length = 211.;
// apply binning hardcoded 360 in phi, 40 in r, rmin = 30cm, rmax=75
  double rbins = 40;
  double rbinwidth = (75.-30.)/rbins;
  double phibinwidth = 2*M_PI/360.;
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      double eion = hiter->second->get_eion();
      unsigned int n_electrons = gsl_ran_poisson(RandomGenerator,eion*electrons_per_gev);
      double dx = (hiter->second->get_x(1) - hiter->second->get_x(0))/n_electrons;
      double dy = (hiter->second->get_y(1) - hiter->second->get_y(0))/n_electrons;
      double dz = (hiter->second->get_z(1) - hiter->second->get_z(0))/n_electrons;
      double dt = (hiter->second->get_t(1) - hiter->second->get_t(0))/n_electrons;
      cout << "layer: " << hiter->second->get_layer() << ", dz: " << dz << endl;
      double x_start = hiter->second->get_x(0) + dx/2.;
      double y_start = hiter->second->get_y(0) + dy/2.;
      double z_start = hiter->second->get_z(0) + dz/2.;
      double t_start = hiter->second->get_t(0) + dt/2.;
      cout << "g4hit created electrons: " << n_electrons 
	   << " from " << eion*1000000 << " keV" << endl;
      for (unsigned int i=0; i<n_electrons; i++)
      {
	// cout << "drift for x: " << x_start
	//       << ", y: " << y_start
	//       << ",z: " << z_start << endl;
	x_start += dx;
	y_start += dy;
	z_start += dz;
	double t_start = (tpc_length/2. - fabs(z_start))/drift_velocity;
// now the drift
	double z_sigma = sqrt(diffusion_long*diffusion_long*(tpc_length/2. - fabs(z_start)));

	double ranphi = gsl_ran_flat(RandomGenerator,-M_PI,M_PI);
	double r_sigma =  sqrt(diffusion_trans*diffusion_trans*(tpc_length/2. - fabs(z_start)));
	double rantrans = gsl_ran_gaussian(RandomGenerator,r_sigma);
	double x_final = x_start + rantrans*cos(ranphi);
	double y_final = y_start + rantrans*sin(ranphi);
	double rantime = gsl_ran_gaussian(RandomGenerator,z_sigma)/drift_velocity;
	double t_final = t_start + rantime;
	double phi = atan2(y_final,x_final);
	double rad = sqrt(x_final*x_final + y_final*y_final);
	// cout << "x_s: " << x_start
	//      << ", x_f: " << x_final
	//      << ", y_s: " << y_start
	//      << ", y_f: " << y_final
	//      << ", t_s: " << t_start
	//      << ", t_f: " << t_final
	//      << ", phi: " << phi
	//      << ", rad: " << rad
	//      << endl;
	if (rad < 30 || rad >75)
	{
	  continue;
	}
	int phibin = (phi+M_PI)/phibinwidth;
	int radbin = (rad-30)/rbinwidth;
	cout << ", phi: " << phi
	     << ", phibin: " << phibin
	     << ", rad: " << rad
	     << ", radbin: " << radbin
	     << endl;
      }
//      gSystem->Exit(0);
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TPCElectronDrift::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TPCElectronDrift::set_seed(const unsigned int iseed)
{
  seed = iseed;
  gsl_rng_set(RandomGenerator,seed);
}

void PHG4TPCElectronDrift::SetDefaultParameters()
{
// Data on gasses @20 C and 760 Torr from the following source:
// http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
double Ne_dEdx = 1.56;    // keV/cm
double CF4_dEdx = 7.00;   // keV/cm
// double Ne_NPrimary = 12;    // Number/cm
// double CF4_NPrimary = 51;   // Number/cm
double Ne_NTotal = 43;     // Number/cm
double CF4_NTotal = 100;   // Number/cm
double TPC_NTot = 0.9*Ne_NTotal + 0.1*CF4_NTotal;
double TPC_dEdx = 0.90 * Ne_dEdx + 0.10 * CF4_dEdx;
double  TPC_ElectronsPerKeV = TPC_NTot / TPC_dEdx;
  set_default_double_param("diffusion_long",0.015); // cm/SQRT(cm)
  set_default_double_param("diffusion_trans",0.006); // cm/SQRT(cm)
  set_default_double_param("drift_velocity",8.0 / 1000.0); // cm/ns
  set_default_double_param("electrons_per_gev",TPC_ElectronsPerKeV*1000000.);
  return;
}
