//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in TpcRS.h.
//
// TpcRS(const std::string &name = "TpcRS")
// everything is keyed to TpcRS, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// TpcRS::~TpcRS()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int TpcRS::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int TpcRS::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int TpcRS::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT;
//   proceed but do not save this event in output (needs output manager setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT;
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int TpcRS::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int TpcRS::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int TpcRS::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int TpcRS::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void TpcRS::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

#include "TpcRS.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TSystem.h>

//____________________________________________________________________________..
TpcRS::TpcRS(const std::string &name)
  : SubsysReco(name)
{
  std::cout << "TpcRS::TpcRS(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
TpcRS::~TpcRS()
{
  delete simulator;
  std::cout << "TpcRS::~TpcRS() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int TpcRS::Init(PHCompositeNode *topNode)
{
  std::cout << "TpcRS::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRS::InitRun(PHCompositeNode *topNode)
{
  std::cout << "TpcRS::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRS::process_event(PHCompositeNode *topNode)
{
  std::cout << "TpcRS::process_event(PHCompositeNode *topNode)" << std::endl;
  PHG4TruthInfoContainer *TruthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!TruthInfo)
  {
    std::cout << PHWHERE << " ERROR: Can't find G4TruthInfo" << std::endl;
    gSystem->Exit(-1);
  }

  PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  if (!hits)
  {
    std::cout << "no G4HITS_TPC node" << std::endl;
    gSystem->Exit(-1);
  }
  PHG4HitContainer::ConstRange hit_range = hits->getHits();
  for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
  {
    PHG4Particle *particle = TruthInfo->GetParticle(hit_iter->second->get_trkid());
    tpcrs::SimulatedHit simu_hit;
    if (hit_iter->second->get_trkid() == 0)
    {
      std::cout << "trackid=0, x: " << hit_iter->second->get_avg_x() << std::endl;
    }
    simu_hit.track_id = hit_iter->second->get_trkid();
    simu_hit.particle_id = particle->get_pid();
    simu_hit.x = hit_iter->second->get_avg_x();
    simu_hit.y = hit_iter->second->get_avg_y();
    simu_hit.z = hit_iter->second->get_avg_z();
    simu_hit.px = 0.5 * (hit_iter->second->get_px(0) + hit_iter->second->get_px(1));
    simu_hit.py = 0.5 * (hit_iter->second->get_py(0) + hit_iter->second->get_py(1));
    simu_hit.pz = 0.5 * (hit_iter->second->get_pz(0) + hit_iter->second->get_pz(1));
    simu_hit.de = hit_iter->second->get_edep();
    simu_hit.ds = hit_iter->second->get_path_length();
    simu_hit.tof = hit_iter->second->get_avg_t();
    simu_hits.push_back(simu_hit);
  }
  simulator->Distort(begin(simu_hits), end(simu_hits), std::back_inserter(dist_hits));

  //  simulator->Digitize(begin(simu_hits), end(simu_hits), std::back_inserter(digi_hits));

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRS::ResetEvent(PHCompositeNode *topNode)
{
  std::cout << "TpcRS::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  simu_hits.clear();
  dist_hits.clear();
  digi_hits.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRS::EndRun(const int runnumber)
{
  std::cout << "TpcRS::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRS::End(PHCompositeNode *topNode)
{
  std::cout << "TpcRS::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRS::Reset(PHCompositeNode *topNode)
{
  std::cout << "TpcRS::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void TpcRS::Print(const std::string &what) const
{
  std::cout << "TpcRS::Print(const std::string &what) const Printing info for " << what << std::endl;
}

void TpcRS::SetupConfigurator(const std::string &filename)
{
  cfg = new tpcrs::Configurator("simple", filename);
  simulator = new tpcrs::Simulator(*cfg);
}
