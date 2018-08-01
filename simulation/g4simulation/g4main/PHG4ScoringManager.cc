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
#include <Geant4/G4VPrimitiveScorer.hh>

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TString.h>

#include <boost/format.hpp>

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

  //4 init IOs
  if (Verbosity() >= VERBOSITY_SOME)
    cout << "PHG4ScoringManager::InitRun - Making PHTFileServer " << m_outputFileName
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

    makeScoringHistograms();

    Fun4AllHistoManager *hm = getHistoManager();
    assert(hm);
    for (unsigned int i = 0; i < hm->nHistos(); i++)
      hm->getHisto(i)->Write();

  }  //   if (not m_outputFileName.empty())

  return Fun4AllReturnCodes::EVENT_OK;
}

//! based on G4VScoreWriter::DumpAllQuantitiesToFile()
void PHG4ScoringManager::makeScoringHistograms()
{
  Fun4AllHistoManager *hm = getHistoManager();
  assert(hm);

  G4ScoringManager *scoringManager = G4ScoringManager::GetScoringManager();
  assert(scoringManager);

  for (unsigned int imesh = 0; imesh < scoringManager->GetNumberOfMesh(); ++imesh)
  {
    G4VScoringMesh *g4mesh = scoringManager->GetMesh(imesh);
    assert(g4mesh);

    const string meshName(g4mesh->GetWorldName().data());
    if (Verbosity())
    {
      cout << "PHG4ScoringManager::makeScoringHistograms - processing mesh " << meshName << endl;
    }

    // descriptors
    G4int nMeshSegments[3];  // number of segments of the mesh
    g4mesh->GetNumberOfSegments(nMeshSegments);

    G4String divisionAxisNames[3];
    g4mesh->GetDivisionAxisNames(divisionAxisNames);

    MeshScoreMap fSMap = g4mesh->GetScoreMap();
    MeshScoreMap::const_iterator msMapItr = fSMap.begin();
    for (; msMapItr != fSMap.end(); msMapItr++)
    {
      G4String psname = msMapItr->first;
      std::map<G4int, G4double *> &score = *(msMapItr->second->GetMap());

      G4double unitValue = g4mesh->GetPSUnitValue(psname);
      G4String unit = g4mesh->GetPSUnit(psname);

      const string hname = boost::str(boost::format("hScore_%1%_%2%") % meshName.data() % psname.data());
      const string htitle = boost::str(boost::format("Mesh %1%, Primitive scorer %2%: score [%3%]") % meshName.c_str() % psname.data() % unit.data());

      if (Verbosity())
        cout << "PHG4ScoringManager::makeScoringHistograms - processing mesh " << meshName
             << "  scorer " << psname
             << "  with axis: "
             << "# i" << divisionAxisNames[0]
             << ", i" << divisionAxisNames[1]
             << ", i" << divisionAxisNames[2]
             << ", value "
             << "[unit: " << unit << "]."
             << " Saving to histogram " << hname << " : " << htitle
             << endl;

      //book histogram
      TH3 *h = new TH3D(hname.c_str(),   //
                        htitle.c_str(),  //
                        nMeshSegments[0],
                        -.5, nMeshSegments[0] - .5,
                        nMeshSegments[1],
                        -.5, nMeshSegments[1] - .5,
                        nMeshSegments[2],
                        -.5, nMeshSegments[2] - .5);
      hm->registerHisto(h);

      h->GetXaxis()->SetTitle(divisionAxisNames[0].data() + TString(" index"));
      h->GetYaxis()->SetTitle(divisionAxisNames[1].data() + TString(" index"));
      h->GetZaxis()->SetTitle(divisionAxisNames[2].data() + TString(" index"));

      // write quantity
      for (int x = 0; x < nMeshSegments[0]; x++)
      {
        for (int y = 0; y < nMeshSegments[1]; y++)
        {
          for (int z = 0; z < nMeshSegments[2]; z++)
          {
            const int idx = x * nMeshSegments[1] * nMeshSegments[2] + y * nMeshSegments[2] + z;

            std::map<G4int, G4double *>::iterator value = score.find(idx);

            if (value != score.end())
            {
              h->SetBinContent(x+1, y+1, z+1, *(value->second) / unitValue);
            }

          }  //           for (int z = 0; z < fNMeshSegments[2]; z++)
        }
      }  //       for (int x = 0; x < nMeshSegments[0]; x++)

    }  //     for (; msMapItr != fSMap.end(); msMapItr++)

  }  //   for (int imesh = 0; imesh < scoringManager->GetNumberOfMesh(); ++imesh)
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
