/*
 * This file is part of KFParticle package
 * Copyright ( C ) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
 *
 * KFParticle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * ( at your option ) any later version.
 *
 * KFParticle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "KFParticle_sPHENIX.h"

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <TFile.h>

typedef std::pair<int, float> particle_pair;

KFParticle_Tools kfpTupleTools_Top;
KFParticle_particleList kfp_list;
std::map<std::string, particle_pair> particleList = kfp_list.getParticleList();

/// KFParticle constructor
KFParticle_sPHENIX::KFParticle_sPHENIX()
  : SubsysReco("KFPARTICLE")
  , m_verbosity(0)
  , m_has_intermediates_sPHENIX(false)
  , m_constrain_to_vertex_sPHENIX(false)
  , m_require_mva(false)
  , m_save_dst(0)
  , m_save_output(1)
  , m_outfile_name("outputData.root")
  , m_outfile(nullptr)
{
}

KFParticle_sPHENIX::KFParticle_sPHENIX(const std::string &name = "KFPARTICLE")
  : SubsysReco(name)
  , m_verbosity(0)
  , m_has_intermediates_sPHENIX(false)
  , m_constrain_to_vertex_sPHENIX(false)
  , m_require_mva(false)
  , m_save_dst(0)
  , m_save_output(1)
  , m_outfile_name("outputData.root")
  , m_outfile(nullptr)
{
}

int KFParticle_sPHENIX::Init(PHCompositeNode *topNode)
{
  if (m_save_output)
  {
    m_outfile = new TFile(m_outfile_name.c_str(), "RECREATE");
    if (m_verbosity > 0) std::cout << "Output nTuple: " << m_outfile_name << std::endl;
    initializeBranches();
  }

  if (m_save_dst) createParticleNode(topNode);

  if (m_require_mva)
  {
    TMVA::Reader *reader;
    std::vector<Float_t> MVA_parValues;
    tie(reader, MVA_parValues) = initMVA();
  }

  for (int i = 0; i < m_num_tracks; ++i)
    if (!particleList.count(m_daughter_name[i]))
    {
      std::cout << "Your track PID, " << m_daughter_name[i] << "is not in the particle list" << std::endl;
      std::cout << "Check KFParticle_particleList.cxx for a list of available particles" << std::endl;
      exit(0);
    }

  return 0;
}

int KFParticle_sPHENIX::process_event(PHCompositeNode *topNode)
{
  std::vector<KFParticle> mother, vertex;
  std::vector<std::vector<KFParticle>> daughters, intermediates;
  int nPVs, multiplicity;

  SvtxVertexMap *check_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vtx_map_node_name);
  if (check_vertexmap->size() == 0)
  {
    if (m_verbosity > 0) std::cout << "KFParticle: Event skipped as there are no vertices" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  SvtxTrackMap *check_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name);
  if (check_trackmap->size() == 0)
  {
    if (m_verbosity > 0) std::cout << "KFParticle: Event skipped as there are no tracks" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  createDecay(topNode, mother, vertex, daughters, intermediates, nPVs, multiplicity);

  if (!m_has_intermediates_sPHENIX) intermediates = daughters;
  if (!m_constrain_to_vertex_sPHENIX) vertex = mother;

  if (mother.size() != 0)
    for (unsigned int i = 0; i < mother.size(); ++i)
    {
      if (m_save_output) fillBranch(topNode, mother[i], vertex[i], daughters[i], intermediates[i], nPVs, multiplicity);
      if (m_save_dst) fillParticleNode(topNode, mother[i], daughters[i], intermediates[i]);

      if (m_verbosity > 0)
      {
        printParticles(mother[i], vertex[i], daughters[i], intermediates[i], nPVs, multiplicity);
        if (m_save_dst) printNode(topNode);
      }
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

int KFParticle_sPHENIX::End(PHCompositeNode *topNode)
{
  if (m_save_output)
  {
    m_outfile->Write();
    m_outfile->Close();
    delete m_outfile;
  }

  return 0;
}

void KFParticle_sPHENIX::printParticles(KFParticle motherParticle,
                                        KFParticle chosenVertex,
                                        std::vector<KFParticle> daughterParticles,
                                        std::vector<KFParticle> intermediateParticles,
                                        int numPVs, int numTracks)
{
  std::cout << "\n---------------KFParticle candidate information---------------" << std::endl;

  std::cout << "Mother information:" << std::endl;
  kfpTupleTools_Top.identify(motherParticle);

  if (m_has_intermediates_sPHENIX)
  {
    std::cout << "Intermediate state information:" << std::endl;
    for (unsigned int i = 0; i < intermediateParticles.size(); i++)
    {
      kfpTupleTools_Top.identify(intermediateParticles[i]);
    }
  }

  std::cout << "Final track information:" << std::endl;
  for (unsigned int i = 0; i < daughterParticles.size(); i++)
  {
    kfpTupleTools_Top.identify(daughterParticles[i]);
  }

  if (m_constrain_to_vertex_sPHENIX)
  {
    std::cout << "Primary vertex information:" << std::endl;
    std::cout << "(x,y,z) = (" << chosenVertex.GetX() << " +/- " << sqrt(chosenVertex.GetCovariance(0, 0)) << ", ";
    std::cout << chosenVertex.GetY() << " +/- " << sqrt(chosenVertex.GetCovariance(1, 1)) << ", ";
    std::cout << chosenVertex.GetZ() << " +/- " << sqrt(chosenVertex.GetCovariance(2, 2)) << ") cm\n"
         << std::endl;
  }

  std::cout << "The number of primary vertices is: " << numPVs << std::endl;
  std::cout << "The number of tracks in the event is: " << numTracks << std::endl;

  std::cout << "------------------------------------------------------------\n"
       << std::endl;
}
