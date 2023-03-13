// $Id: $

/*!
 * \file PHG4ScoringManager.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4ScoringManager.h"

#include "PHG4InEvent.h"
#include "PHG4Particle.h"

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <fun4all/Fun4AllBase.h>  // for Fun4AllBase::VERBO...
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <phool/getClass.h>

#include <HepMC/SimpleVector.h>  // for FourVector

#include <TAxis.h>  // for TAxis
#include <TDatabasePDG.h>
#include <TH1.h>
#include <TH3.h>           // for TH3, TH3D
#include <TNamed.h>        // for TNamed
#include <TParticlePDG.h>  // for TParticlePDG
#include <TVector3.h>

#include <Geant4/G4RunManager.hh>
#include <Geant4/G4ScoringManager.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4THitsMap.hh>     // for G4THitsMap
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Types.hh>        // for G4int, G4double
#include <Geant4/G4UImanager.hh>
#include <Geant4/G4VScoringMesh.hh>  // for G4VScoringMesh
#include <Geant4/G4Version.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/format.hpp>
#pragma GCC diagnostic pop

#include <cassert>
#include <cmath>  // for fabs, M_PI
#include <iostream>
#include <limits>   // for numeric_limits
#include <map>      // for _Rb_tree_const_ite...
#include <utility>  // for pair

using namespace std;

PHG4ScoringManager::PHG4ScoringManager()
  : SubsysReco("PHG4ScoringManager")
{
}

int PHG4ScoringManager::InitRun(PHCompositeNode */*topNode*/)
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
  h->GetXaxis()->SetBinLabel(i++, "Collision count");
  h->GetXaxis()->SetBinLabel(i++, "Event count accepted");
  h->GetXaxis()->SetBinLabel(i++, "Collision count accepted");
  //  h->GetXaxis()->SetBinLabel(i++, "G4Hit count");
  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  hm->registerHisto(new TH1D("hNChEta",  //
                             "Charged particle #eta distribution;#eta;Count",
                             1000, -5, 5));

  hm->registerHisto(new TH1D("hVertexZ",  //
                             "Vertex z distribution;z [cm];Count",
                             10000, m_vertexHistRange.first, m_vertexHistRange.second));

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

  PHHepMCGenEventMap *geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!geneventmap)
  {
    static bool once = true;
    if (once)
    {
      once = false;
      cout << "PHG4ScoringManager::process_event - - missing node PHHepMCGenEventMap. Skipping HepMC stat." << std::endl;
    }
  }
  else
  {
    h_norm->Fill("Collision count", geneventmap->size());

    TH1D *hVertexZ = dynamic_cast<TH1D *>(hm->getHisto("hVertexZ"));
    assert(hVertexZ);

    for (const auto genevntpair : geneventmap->get_map())
    {
      const PHHepMCGenEvent *genevnt = genevntpair.second;
      assert(genevnt);

      if (genevnt->get_collision_vertex().z() < m_vertexAcceptanceRange.first
          or genevnt->get_collision_vertex().z() > m_vertexAcceptanceRange.second)
      {
        if (Verbosity() >= 2)
        {
          cout <<__PRETTY_FUNCTION__<<": get vertex "<<genevnt->get_collision_vertex().z()
              <<" which is outside range "<<m_vertexAcceptanceRange.first <<" to "<<m_vertexAcceptanceRange.second<<" cm:";
          genevnt->identify();
        }

        return Fun4AllReturnCodes::ABORTEVENT;
      }

      hVertexZ->Fill(genevnt->get_collision_vertex().z());
    }

    h_norm->Fill("Collision count accepted", geneventmap->size());
  }

  TH1D *hNChEta = dynamic_cast<TH1D *>(hm->getHisto("hNChEta"));
  assert(hNChEta);

  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!ineve)
  {
    cout << "PHG4ScoringManager::process_event - Error - "
         << "unable to find DST node "
         << "PHG4INEVENT" << endl;
  }
  else
  {
    const auto primary_range = ineve->GetParticles();
    for (auto particle_iter = primary_range.first;
         particle_iter != primary_range.second;
         ++particle_iter)
    {
      const PHG4Particle *p = particle_iter->second;
      assert(p);
      TParticlePDG *pdg_p = TDatabasePDG::Instance()->GetParticle(p->get_pid());
      assert(pdg_p);
      if (fabs(pdg_p->Charge()) > 0)
      {
        TVector3 pvec(p->get_px(), p->get_py(), p->get_pz());
        if (pvec.Perp2() > 0)
        {
          assert(hNChEta);
          hNChEta->Fill(pvec.PseudoRapidity());
        }
      }
    }  //          if (_load_all_particle) else
  }

  h_norm->Fill("Event count accepted", 1);

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
      cout << "PHG4ScoringManager::makeScoringHistograms - processing mesh " << meshName << ": " << endl;
      g4mesh->List();
    }

    // descriptors
    G4int nMeshSegments[3];  // number of segments of the mesh
    g4mesh->GetNumberOfSegments(nMeshSegments);

    G4String divisionAxisNames[3];
    g4mesh->GetDivisionAxisNames(divisionAxisNames);

    // process shape
    const G4ThreeVector meshSize = g4mesh->GetSize();
    const G4ThreeVector meshTranslate = g4mesh->GetTranslation();
#if G4VERSION_NUMBER >= 1060
    const G4VScoringMesh::MeshShape meshShape = g4mesh->GetShape();
#else
    const MeshShape meshShape = g4mesh->GetShape();
#endif
    //PHENIX units
    vector<double> meshBoundMin = {std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN()};
    //PHENIX units
    vector<double> meshBoundMax = {std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN()};
#if G4VERSION_NUMBER >= 1060
    if (meshShape == G4VScoringMesh::MeshShape::box)
#else
    if (meshShape == boxMesh)
#endif
    {
      meshBoundMin[0] = (-meshSize[0] + meshTranslate[0]) / cm;
      meshBoundMax[0] = (meshSize[0] + meshTranslate[0]) / cm;
      meshBoundMin[1] = (-meshSize[1] + meshTranslate[1]) / cm;
      meshBoundMax[1] = (meshSize[1] + meshTranslate[1]) / cm;
      meshBoundMin[2] = (-meshSize[2] + meshTranslate[2]) / cm;
      meshBoundMax[2] = (meshSize[2] + meshTranslate[2]) / cm;

      divisionAxisNames[0] += " [cm]";
      divisionAxisNames[1] += " [cm]";
      divisionAxisNames[2] += " [cm]";
    }
#if G4VERSION_NUMBER >= 1060
    else if (meshShape == G4VScoringMesh::MeshShape::cylinder)
#else
    else if (meshShape == cylinderMesh)
#endif
    {
      //      fDivisionAxisNames[0] = "Z";
      //      fDivisionAxisNames[1] = "PHI";
      //      fDivisionAxisNames[2] = "R";
      //      G4VSolid * tubsSolid = new G4Tubs(tubsName+"0", // name
      //                0.,           // R min
      //                fSize[0],     // R max
      //                fSize[1],     // Dz
      //                0.,           // starting phi
      //                                        twopi*rad);   // segment phi
      meshBoundMin[0] = (-meshSize[1] + meshTranslate[0]) / cm;
      meshBoundMax[0] = (meshSize[1] + meshTranslate[0]) / cm;
      meshBoundMin[1] = 0;
      meshBoundMax[1] = 2 * M_PI;
      meshBoundMin[2] = 0;
      meshBoundMax[2] = meshSize[0] / cm;

      divisionAxisNames[0] += " [cm]";
      divisionAxisNames[1] += " [rad]";
      divisionAxisNames[2] += " [cm]";
    }
    else
    {
      cout << "PHG4ScoringManager::makeScoringHistograms - Error - unsupported mesh shape " << (int) meshShape << ". Skipping this mesh!" << endl;
      g4mesh->List();
      continue;
    }

#if G4VERSION_NUMBER >= 1060
    G4VScoringMesh::MeshScoreMap fSMap = g4mesh->GetScoreMap();
    G4VScoringMesh::MeshScoreMap::const_iterator msMapItr = fSMap.begin();
#else
    MeshScoreMap fSMap = g4mesh->GetScoreMap();
    MeshScoreMap::const_iterator msMapItr = fSMap.begin();
#endif
    for (; msMapItr != fSMap.end(); ++msMapItr)
    {
      G4String psname = msMapItr->first;
#if G4VERSION_NUMBER >= 1040
      std::map<G4int, G4StatDouble *> &score = *(msMapItr->second->GetMap());
#else
      std::map<G4int, G4double *> &score = *(msMapItr->second->GetMap());
#endif
      G4double unitValue = g4mesh->GetPSUnitValue(psname);
      G4String unit = g4mesh->GetPSUnit(psname);

      const string hname = boost::str(boost::format("hScore_%1%_%2%") % meshName.data() % psname.data());
      const string htitle = boost::str(boost::format("Mesh %1%, Primitive scorer %2%: score [%3%]") % meshName.c_str() % psname.data() % unit.data());

      if (Verbosity())
      {
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
      }
      //book histogram
      TH3 *h = new TH3D(hname.c_str(),   //
                        htitle.c_str(),  //
                        nMeshSegments[0],
                        meshBoundMin[0], meshBoundMax[0],
                        nMeshSegments[1],
                        meshBoundMin[1], meshBoundMax[1],
                        nMeshSegments[2],
                        meshBoundMin[2], meshBoundMax[2]);
      hm->registerHisto(h);

      h->GetXaxis()->SetTitle(divisionAxisNames[0].data());
      h->GetYaxis()->SetTitle(divisionAxisNames[1].data());
      h->GetZaxis()->SetTitle(divisionAxisNames[2].data());

      // write quantity
      for (int x = 0; x < nMeshSegments[0]; x++)
      {
        for (int y = 0; y < nMeshSegments[1]; y++)
        {
          for (int z = 0; z < nMeshSegments[2]; z++)
          {
            const int idx = x * nMeshSegments[1] * nMeshSegments[2] + y * nMeshSegments[2] + z;

#if G4VERSION_NUMBER >= 1040
            std::map<G4int, G4StatDouble *>::iterator value = score.find(idx);
#else
            std::map<G4int, G4double *>::iterator value = score.find(idx);
#endif
            if (value != score.end())
            {
#if G4VERSION_NUMBER >= 1040
              h->SetBinContent(x + 1, y + 1, z + 1, (value->second->mean()) / unitValue);
#else
              h->SetBinContent(x + 1, y + 1, z + 1, *(value->second) / unitValue);
#endif
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
  if (!hm)
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
