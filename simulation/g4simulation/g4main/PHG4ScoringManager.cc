// $Id: $

/*!
 * \file PHG4ScoringManager.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4ScoringManager.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <Geant4/G4RunManager.hh>
#include <Geant4/G4ScoringManager.hh>
#include <Geant4/G4UImanager.hh>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <cassert>
#include <iostream>

using namespace std;

PHG4ScoringManager::PHG4ScoringManager()
  : SubsysReco("PHG4ScoringManager")
{
}

PHG4ScoringManager::~PHG4ScoringManager()
{
}

//_________________________________________________________________
int PHG4ScoringManager::Init(PHCompositeNode *topNode)
{
  return 0;
}

int PHG4ScoringManager::InitRun(PHCompositeNode *topNode)
{
  //1. check G4RunManager
  G4RunManager *runManager = G4RunManager::GetRunManager();
  if (runManager == nullptr)
  {
    cout << "PHG4ScoringManager::InitRun - fatal error: G4RunManager was not initialized yet. Please do include the Geant4 simulation in this Fun4All run." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  //2. Init scoring manager
  G4ScoringManager *scoringManager = G4ScoringManager::GetScoringManager();
  assert(scoringManager);

  //3. run scoring commands
  G4UImanager *UImanager = G4UImanager::GetUIpointer();
  assert(UImanager);

  for (const string &cmd : m_commands)
  {
    if (Verbosity() >= VERBOSITY_SOME)
    {
      cout << "PHG4ScoringManager::InitRun - execute Geatn4 command: " << cmd << endl;
    }
    UImanager->ApplyCommand(cmd.c_str());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_________________________________________________________________
void PHG4ScoringManager::G4Command(const string &cmd)
{
  m_commands.push_back(cmd);
  return;
}
//_________________________________________________________________

//_________________________________________________________________
int PHG4ScoringManager::process_event(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4ScoringManager::ResetEvent(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//_________________________________________________________________
int PHG4ScoringManager::End(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
