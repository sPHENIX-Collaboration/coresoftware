/*
 * This file is part of KFParticle package
 * Copyright (C) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
 *
 * KFParticle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KFParticle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef KFPARTICLESPHENIX_KFPARTICLESPHENIX_H
#define KFPARTICLESPHENIX_KFPARTICLESPHENIX_H

#include "KFParticle_eventReconstruction.h"
#include "KFParticle_nTuple.h"
#include "KFParticle_DST.h"

//sPHENIX stuff
#include <fun4all/SubsysReco.h>

//KFParticle stuff
#include <KFParticle.h>

#include <map>
#include <string>

class TFile;

class KFParticle_sPHENIX : public SubsysReco, public KFParticle_nTuple, public KFParticle_DST, protected KFParticle_eventReconstruction
{
 public:
  KFParticle_sPHENIX();

  KFParticle_sPHENIX(const std::string &name);

  virtual ~KFParticle_sPHENIX(){}

  int Init(PHCompositeNode *topNode);

  int process_event(PHCompositeNode *topNode);

  int End(PHCompositeNode *topNode);

  /**
   * If verbosity is > 0, this will print out all candidate information:
   * masses, momenta and positions for mothers, intermediates and final state tracks,
   * PV position, number of vertices and number of tracks in the event (multiplicity)
   */
  void printParticles(KFParticle motherParticle,
                      KFParticle chosenVertex,
                      vector<KFParticle> daughterParticles,
                      vector<KFParticle> intermediateParticles,
                      int numPVs, int numTracks);

  ///Parameters for the user to vary
  static const int max_particles = 99;

  void setMotherName(const std::string mother_name)
  {
    m_mother_name = mother_name;
    m_mother_name_Tools = mother_name;
  }

  void useIntermediateName(bool use_intermediate_name) { m_use_intermediate_name = use_intermediate_name; }

  void hasIntermediateStates(bool has_intermediates)
  {
    m_has_intermediates = has_intermediates;
    m_has_intermediates_nTuple = has_intermediates;
    m_has_intermediates_sPHENIX = has_intermediates;
    m_has_intermediates_DST = has_intermediates;
  }

  void setNumberOfTracks(int num_tracks)
  {
    m_num_tracks = num_tracks;
    m_num_tracks_nTuple = num_tracks;
  }

  void setNumberTracksFromIntermeditateState(int num_tracks[max_particles])
  {
    for (int i = 0; i < max_particles; ++i) m_num_tracks_from_intermediate[i] = num_tracks[i];
  }

  void setNumberOfIntermediateStates(int n_intermediates)
  {
    m_num_intermediate_states = n_intermediates;
    m_num_intermediate_states_nTuple = n_intermediates;
  }

  void getChargeConjugate(bool get_charge_conjugate)
  {
    m_get_charge_conjugate_nTuple = get_charge_conjugate;
    m_get_charge_conjugate = get_charge_conjugate;
  }

  void setDaughters(std::pair<std::string, int> daughter_list[max_particles])
  {
    for (int i = 0; i < max_particles; ++i)
    {
      m_daughter_name_evt[i] = daughter_list[i].first;
      m_daughter_charge_evt[i] = daughter_list[i].second;
      m_daughter_name[i] = daughter_list[i].first;
      m_daughter_charge[i] = daughter_list[i].second;
    }
  }

  void setIntermediateStates(std::pair<std::string, int> intermediate_list[max_particles])
  {
    for (int i = 0; i < max_particles; ++i)
    {
      m_intermediate_name_ntuple[i] = intermediate_list[i].first;
      //m_intermediate_name_evt[i] = intermediate_list[i].first;
      //m_intermediate_charge_evt[i] = intermediate_list[i].second;
      m_intermediate_name[i] = intermediate_list[i].first;
      m_intermediate_charge[i] = intermediate_list[i].second;
    }
  }

  void setMinimumMass(float min_mass) { m_min_mass = min_mass; }

  void setMaximumMass(float max_mass) { m_max_mass = max_mass; }

  void setMinimumLifetime(float min_lifetime) { m_min_lifetime = min_lifetime; }

  void setMaximumLifetime(float max_lifetime) { m_max_lifetime = max_lifetime; }

  void setMinimumTrackPT(float pt) { m_track_pt = pt; }

  void setMaximumTrackPTchi2(float ptchi2) { m_track_ptchi2 = ptchi2; }

  void setMinimumTrackIPchi2(float ipchi2) { m_track_ipchi2 = ipchi2; }

  void setMaximumTrackchi2nDOF(float trackchi2ndof) { m_track_chi2ndof = trackchi2ndof; }

  void setMaximumDaughterDCA(float dca) { m_comb_DCA = dca; }

  void setMaximumVertexchi2nDOF(float vertexchi2nDOF) { m_vertex_chi2ndof = vertexchi2nDOF; }

  void setFlightDistancechi2(float fdchi2) { m_fdchi2 = fdchi2; }

  void setMinDIRA(float dira_min) { m_dira_min = dira_min; }

  void setMaxDIRA(float dira_max) { m_dira_max = dira_max; }

  void setMotherPT(float mother_pt) { m_mother_pt = mother_pt; }

  void setMotherIPchi2(float mother_ipchi2) { m_mother_ipchi2 = mother_ipchi2; }

  void constrainToPrimaryVertex(bool constrain_to_vertex)
  {
    m_constrain_to_vertex = constrain_to_vertex;
    m_constrain_to_vertex_nTuple = constrain_to_vertex;
    m_constrain_to_vertex_sPHENIX = constrain_to_vertex;
  }

  void constrainIntermediateMasses(bool constrain_int_mass) { m_constrain_int_mass = constrain_int_mass; }

  void setIntermediateMassRange(std::pair<float, float> intermediate_mass_range[max_particles])
  {
    for (int i = 0; i < max_particles; ++i) m_intermediate_mass_range[i] = intermediate_mass_range[i];
  }

  void setIntermediateMinPT(float intermediate_min_pt[max_particles])
  {
    for (int i = 0; i < max_particles; ++i) m_intermediate_min_pt[i] = intermediate_min_pt[i];
  }

  void useMVA(bool require_mva) { m_require_mva = require_mva; }

  void setNumMVAPars(unsigned int nPars) { m_nPars = nPars; }

  void setMVAVarList(std::string mva_variable_list[max_particles])
  {
    for (int i = 0; i < max_particles; ++i) m_mva_variable_list[i] = mva_variable_list[i];
  }

  void setMVAType(const std::string mva_type) { m_mva_type = mva_type; }

  void setMVAWeightsPath(const std::string mva_weights_path) { m_mva_path = mva_weights_path; }

  void setMVACutValue(float cut_value) { m_mva_cut_value = cut_value; }

  void Verbosity(int verbosity) { m_verbosity = verbosity; }

  void saveDST(bool save) { m_save_dst = save; }

  void saveTrackContainer(bool save) { m_write_track_container = save; }

  void saveParticleContainer(bool save) { m_write_particle_container = save; }

  void setContainerName(const std::string name) { m_container_name = name; }

  void saveOutput(bool save) { m_save_output = save; }

  void setOutputName(const std::string name) { m_outfile_name = name; }

  void doTruthMatching(bool truth) { m_truth_matching = truth; }

  void getDetectorInfo(bool detinfo) { m_detector_info = detinfo; }

  ///Use alternate vertex and track fitters
  void setVertexMapNodeName(std::string vtx_map_node_name) { m_vtx_map_node_name = vtx_map_node_name; }

  ///Use alternate vertex and track fitters
  void setTrackMapNodeName(std::string trk_map_node_name) { m_trk_map_node_name = trk_map_node_name; }

 private:
  bool m_verbosity;
  bool m_has_intermediates_sPHENIX;
  bool m_constrain_to_vertex_sPHENIX;
  bool m_require_mva;
  bool m_save_dst;
  bool m_save_output;
  std::string m_outfile_name;
  TFile *m_outfile;
};

#endif  //KFPARTICLESPHENIX_KFPARTICLESPHENIX_H
