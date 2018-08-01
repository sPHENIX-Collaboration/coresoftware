// $Id: $

/*!
 * \file PHG4ScoringManager.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4ScoringManager.h"

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <Geant4/G4RunManager.hh>
#include <Geant4/G4ScoringManager.hh>
#include <Geant4/G4UImanager.hh>

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TH3.h>
#include <TString.h>

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
  if (Verbosity() >= VERBOSITY_SOME)
    cout << "PHG4ScoringManager::get_HistoManager - Making PHTFileServer " << m_outputFileName
         << endl;
  PHTFileServer::get().open(m_outputFileName, "RECREATE");

  Fun4AllHistoManager *hm = getHistoManager();
  assert(hm);
  TH1D *h = new TH1D("hNormalization",  //
                     "Normalization;Items;Summed quantity", 10, .5, 10.5);
  int i = 1;
  h->GetXaxis()->SetBinLabel(i++, "Event count");
  //  h->GetXaxis()->SetBinLabel(i++, "Collision count");
  //  h->GetXaxis()->SetBinLabel(i++, "G4Hit count");
  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  return Fun4AllReturnCodes::EVENT_OK;
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
  Fun4AllHistoManager *hm = getHistoManager();
  assert(hm);

  TH1D *h_norm = dynamic_cast<TH1D *>(hm->getHisto("hNormalization"));
  assert(h_norm);

  h_norm->Fill("Event count", 1);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4ScoringManager::ResetEvent(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//_________________________________________________________________
int PHG4ScoringManager::End(PHCompositeNode *)
{
  if (not m_outputFileName.empty())
  {
    PHTFileServer::get().cd(m_outputFileName);
    cout << "PHG4ScoringManager::End - save results to " << m_outputFileName << endl;

    Fun4AllHistoManager *hm = getHistoManager();
    assert(hm);
    for (unsigned int i = 0; i < hm->nHistos(); i++)
      hm->getHisto(i)->Write();
  }  //   if (not m_outputFileName.empty())

  return Fun4AllReturnCodes::EVENT_OK;
}

Fun4AllHistoManager *
PHG4ScoringManager::getHistoManager()
{
  static string histname("PHG4ScoringManager_HISTOS");
  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllHistoManager *hm = se->getHistoManager(histname);
  if (not hm)
  {
    if (Verbosity())
      cout
          << "PHG4ScoringManager::get_HistoManager - Making Fun4AllHistoManager " << histname
          << endl;
    hm = new Fun4AllHistoManager(histname);
    se->registerHistoManager(hm);
  }
  assert(hm);
  return hm;
}
