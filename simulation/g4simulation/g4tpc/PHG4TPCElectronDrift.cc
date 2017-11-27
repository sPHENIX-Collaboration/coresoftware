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
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      double eion = hiter->second->get_eion();
      unsigned int n_electrons = gsl_ran_poisson(RandomGenerator,eion*electrons_per_gev);
      double dx = (hiter->second->get_x(1) - hiter->second->get_x(0))/n_electrons;
      double dy = (hiter->second->get_y(1) - hiter->second->get_y(0))/n_electrons;
      double dz = (hiter->second->get_z(1) - hiter->second->get_z(0))/n_electrons;
      cout << "layer: " << hiter->second->get_layer() << ", dz: " << dz << endl;
      double x_start = hiter->second->get_x(0) + dx/2.;
      double y_start = hiter->second->get_y(0) + dy/2.;
      double z_start = hiter->second->get_z(0) + dz/2.;
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
  set_default_double_param("electrons_per_gev",30.*1000000.);
  return;
}
