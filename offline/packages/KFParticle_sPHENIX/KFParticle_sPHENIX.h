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

// There is something broken in this package. clang format puts
// #include "KFParticle_DST.h" first if there is no space but then
// loading KFParticle makes root barf. I have no time to track this down
// right now, it must be something in KFParticle_eventReconstruction.h and
// the include files it uses (some include guard misfiring?)

#include "KFParticle_DST.h"
#include "KFParticle_nTuple.h"

// sPHENIX stuff
#include <fun4all/SubsysReco.h>

// KFParticle stuff
#include <KFParticle.h>

#include <algorithm>  // for max
#include <memory>     // for allocator_traits<>::valu...
#include <string>
#include <utility>  // for pair
#include <vector>   // for vector

class PHCompositeNode;
class TFile;

class KFParticle_sPHENIX : public SubsysReco, public KFParticle_nTuple, public KFParticle_DST, protected KFParticle_eventReconstruction
{
 public:
  KFParticle_sPHENIX();

  explicit KFParticle_sPHENIX(const std::string &name);

  ~KFParticle_sPHENIX() override = default;

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  /**
   * If verbosity is > 0, this will print out all candidate information:
   * masses, momenta and positions for mothers, intermediates and final state tracks,
   * PV position, number of vertices and number of tracks in the event (multiplicity)
   */
  void printParticles(const KFParticle &motherParticle,
                      const KFParticle &chosenVertex,
                      const std::vector<KFParticle> &daughterParticles,
                      const std::vector<KFParticle> &intermediateParticles,
                      const int numPVs, const int numTracks);

  int parseDecayDescriptor();

  /// Parameters for the user to vary

  void setDecayDescriptor(const std::string &decayDescriptor) { m_decayDescriptor = decayDescriptor; }

  static const int max_particles = 99;

  void setMotherName(const std::string &mother_name)
  {
    m_mother_name = mother_name;
    m_mother_name_Tools = mother_name;
  }

  void hasIntermediateStates(bool has_intermediates = true)
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

  void setNumberTracksFromIntermeditateState(const std::vector<int> &num_tracks)
  {
    for (unsigned int i = 0; i < num_tracks.size(); ++i)
    {
      m_num_tracks_from_intermediate.push_back(num_tracks[i]);
      m_num_tracks_from_intermediate_nTuple.push_back(num_tracks[i]);
    }
  }

  void setNumberOfIntermediateStates(int n_intermediates)
  {
    m_num_intermediate_states = n_intermediates;
    m_num_intermediate_states_nTuple = n_intermediates;
  }

  void getChargeConjugate(bool get_charge_conjugate = true)
  {
    m_get_charge_conjugate_nTuple = get_charge_conjugate;
    m_get_charge_conjugate = get_charge_conjugate;
  }

  void setDaughters(std::vector<std::pair<std::string, int>> &daughter_list)
  {
    for (unsigned int i = 0; i < daughter_list.size(); ++i)
    {
      m_daughter_name.push_back(daughter_list[i].first);
      m_daughter_charge.push_back(daughter_list[i].second);
    }
  }

  void setIntermediateStates(const std::vector<std::pair<std::string, int>> &intermediate_list)
  {
    for (unsigned int i = 0; i < intermediate_list.size(); ++i)
    {
      m_intermediate_name_ntuple.push_back(intermediate_list[i].first);
      m_intermediate_name.push_back(intermediate_list[i].first);
      m_intermediate_charge.push_back(intermediate_list[i].second);
    }
  }

  void setMinimumMass(float min_mass) { m_min_mass = min_mass; }

  void setMaximumMass(float max_mass) { m_max_mass = max_mass; }

  void setDecayTimeRange_XY(float min_decayTime, float max_decayTime)
  {
    m_min_decayTime_xy = min_decayTime;
    m_max_decayTime_xy = max_decayTime;
  }

  void setDecayLengthRange_XY(float min_decayLength, float max_decayLength)
  {
    m_min_decayLength_xy = min_decayLength;
    m_max_decayLength_xy = max_decayLength;
  }

  void setDecayTimeRange(float min_decayTime, float max_decayTime)
  {
    m_min_decayTime = min_decayTime;
    m_max_decayTime = max_decayTime;
  }

  void setDecayLengthRange(float min_decayLength, float max_decayLength)
  {
    m_min_decayLength = min_decayLength;
    m_max_decayLength = max_decayLength;
  }

  void setMinDecayTimeSignificance(float min = 0) { m_mother_min_decay_time_significance = min; }

  void setMinDecayLengthSignificance(float min = 0) { m_mother_min_decay_length_significance = min; }

  void setMinDecayLengthSignificance_XY(float min = 0) { m_mother_min_decay_length_xy_significance = min; }

  void setMinimumTrackPT(float pt) { m_track_pt = pt; }

  void setMaximumTrackPTchi2(float ptchi2) { m_track_ptchi2 = ptchi2; }

  void setMinimumTrackIP_XY(float ip) { m_track_ip_xy = ip; }

  void setMinimumTrackIPchi2_XY(float ipchi2) { m_track_ipchi2_xy = ipchi2; }

  void setMinimumTrackIP(float ip) { m_track_ip = ip; }

  void setMinimumTrackIPchi2(float ipchi2) { m_track_ipchi2 = ipchi2; }

  void setMaximumTrackchi2nDOF(float trackchi2ndof) { m_track_chi2ndof = trackchi2ndof; }

  void setMinMVTXhits(int nHits) { m_nMVTXStates = nHits; } //Actually state counting but use this for backwards compatibility!

  void setMinINTThits(int nHits) { m_nINTTStates = nHits; } //Actually state counting but use this for backwards compatibility!

  void setMinTPChits(int nHits) { m_nTPCStates = nHits; } //Actually state counting but use this for backwards compatibility!

  void setMinTPOThits(int nHits) { m_nTPCStates = nHits; } //Actually state counting but use this for backwards compatibility!

  void setMaximumDaughterDCA_XY(float dca) { m_comb_DCA_xy = dca; }

  void setMaximumDaughterDCA(float dca) { m_comb_DCA = dca; }
 
  void setMinimumRadialSV(float min_rad_sv) { m_min_radial_SV = min_rad_sv; }

  void setMaximumVertexchi2nDOF(float vertexchi2nDOF) { m_vertex_chi2ndof = vertexchi2nDOF; }

  void setFlightDistancechi2(float fdchi2) { m_fdchi2 = fdchi2; }

  void setMinDIRA(float dira_min) { m_dira_min = dira_min; }

  void setMaxDIRA(float dira_max) { m_dira_max = dira_max; }

  void setMinDIRA_XY(float dira_min) { m_dira_xy_min = dira_min; }

  void setMaxDIRA_XY(float dira_max) { m_dira_xy_max = dira_max; }

  void setMotherPT(float mother_pt) { m_mother_pt = mother_pt; }

  void setMotherIP(float mother_ip) { m_mother_ip = mother_ip; }

  void setMotherIP_XY(float mother_ip) { m_mother_ip_xy = mother_ip; }

  void setMotherIPchi2(float mother_ipchi2) { m_mother_ipchi2 = mother_ipchi2; }

  void setMotherIPchi2_XY(float mother_ipchi2) { m_mother_ipchi2_xy = mother_ipchi2; }

  void setMaximumMotherVertexVolume(float vertexvol) { m_mother_vertex_volume = vertexvol; }

  void constrainToPrimaryVertex(bool constrain_to_vertex = true)
  {
    m_constrain_to_vertex = constrain_to_vertex;
    m_constrain_to_vertex_nTuple = constrain_to_vertex;
    m_constrain_to_vertex_sPHENIX = constrain_to_vertex;
  }

  void useMbdVertex(bool use = true)
  {
    m_use_mbd_vertex = use;
    m_use_mbd_vertex_truth = use;
  }

  void dontUseGlobalVertex(bool dont = true) { m_dont_use_global_vertex = m_dont_use_global_vertex_truth  = dont; }

  void useFakePrimaryVertex(bool use_fake = true)
  {
    m_use_fake_pv = use_fake;
    m_use_fake_pv_nTuple = use_fake;
  }

  void allowZeroMassTracks(bool allow = true) { m_allowZeroMassTracks = allow; }

  void extraolateTracksToSV(bool extrapolate = true)
  {
    m_extrapolateTracksToSV = extrapolate;
    m_extrapolateTracksToSV_nTuple = extrapolate;
  }

  void constrainIntermediateMasses(bool constrain_int_mass = true) { m_constrain_int_mass = constrain_int_mass; }

  void setIntermediateMassRange(const std::vector<std::pair<float, float>> &intermediate_mass_range)
  {
    for (unsigned int i = 0; i < intermediate_mass_range.size(); ++i) m_intermediate_mass_range.push_back(intermediate_mass_range[i]);
  }

  void setIntermediateMinPT(const std::vector<float> &intermediate_min_pt)
  {
    m_intermediate_min_pt = intermediate_min_pt;
  }

  void setIntermediateMinIP_XY(const std::vector<float> &intermediate_min_IP)
  {
    for (unsigned int i = 0; i < intermediate_min_IP.size(); ++i) m_intermediate_min_ip_xy.push_back(intermediate_min_IP[i]);
  }

  void setIntermediateIPRange_XY(const std::vector<std::pair<float, float> /*unused*/> &intermediate_IP_range)
  {
    for (unsigned int i = 0; i < intermediate_IP_range.size(); ++i)
    {
      m_intermediate_min_ip_xy.push_back(intermediate_IP_range[i].first);
      m_intermediate_max_ip_xy.push_back(intermediate_IP_range[i].second);
    }
  }

  void setIntermediateMinIP(const std::vector<float> &intermediate_min_IP)
  {
    for (unsigned int i = 0; i < intermediate_min_IP.size(); ++i) m_intermediate_min_ip.push_back(intermediate_min_IP[i]);
  }

  void setIntermediateIPRange(const std::vector<std::pair<float, float> /*unused*/> &intermediate_IP_range)
  {
    for (unsigned int i = 0; i < intermediate_IP_range.size(); ++i)
    {
      m_intermediate_min_ip.push_back(intermediate_IP_range[i].first);
      m_intermediate_max_ip.push_back(intermediate_IP_range[i].second);
    }
  }

  void setIntermediateMinIPchi2_XY(const std::vector<float> &intermediate_min_IPchi2)
  {
    for (unsigned int i = 0; i < intermediate_min_IPchi2.size(); ++i) m_intermediate_min_ipchi2_xy.push_back(intermediate_min_IPchi2[i]);
  }

  void setIntermediateIPchi2Range_XY(const std::vector<std::pair<float, float> /*unused*/> &intermediate_IPchi2_range)
  {
    for (unsigned int i = 0; i < intermediate_IPchi2_range.size(); ++i)
    {
      m_intermediate_min_ipchi2_xy.push_back(intermediate_IPchi2_range[i].first);
      m_intermediate_max_ipchi2_xy.push_back(intermediate_IPchi2_range[i].second);
    }
  }

  void setIntermediateMinIPchi2(const std::vector<float> &intermediate_min_IPchi2)
  {
    for (unsigned int i = 0; i < intermediate_min_IPchi2.size(); ++i) m_intermediate_min_ipchi2.push_back(intermediate_min_IPchi2[i]);
  }

  void setIntermediateIPchi2Range(const std::vector<std::pair<float, float> /*unused*/> &intermediate_IPchi2_range)
  {
    for (unsigned int i = 0; i < intermediate_IPchi2_range.size(); ++i)
    {
      m_intermediate_min_ipchi2.push_back(intermediate_IPchi2_range[i].first);
      m_intermediate_max_ipchi2.push_back(intermediate_IPchi2_range[i].second);
    }
  }

  void setIntermediateMinDIRA(const std::vector<float> &intermediate_min_DIRA)
  {
    for (unsigned int i = 0; i < intermediate_min_DIRA.size(); ++i) m_intermediate_min_dira.push_back(intermediate_min_DIRA[i]);
  }

  void setIntermediateMinFDchi2(const std::vector<float> &intermediate_min_FDchi2)
  {
    for (unsigned int i = 0; i < intermediate_min_FDchi2.size(); ++i) m_intermediate_min_fdchi2.push_back(intermediate_min_FDchi2[i]);
  }

  void setIntermediateMaxVertexVolume(const std::vector<float> &intermediate_max_vertexvol)
  {
    for (unsigned int i = 0; i < intermediate_max_vertexvol.size(); ++i) m_intermediate_vertex_volume.push_back(intermediate_max_vertexvol[i]);
  }

  void use2Dmatching(bool use_2D_matching_tools = true) { m_use_2D_matching_tools = use_2D_matching_tools; }

  void useMVA(bool require_mva = true) { m_require_mva = require_mva; }

  void setNumMVAPars(unsigned int nPars) { m_nPars = nPars; }

  void setMVAVarList(const std::vector<std::string> &mva_variable_list)
  {
    for (unsigned int i = 0; i < mva_variable_list.size(); ++i) m_mva_variable_list.push_back(mva_variable_list[i]);
  }

  void setMVAType(const std::string &mva_type) { m_mva_type = mva_type; }

  void setMVAWeightsPath(const std::string &mva_weights_path) { m_mva_path = mva_weights_path; }

  void setMVACutValue(float cut_value) { m_mva_cut_value = cut_value; }

  void saveDST(bool save = true) { m_save_dst = save; }

  void saveTrackContainer(bool save = true) { m_write_track_container = save; }

  void saveParticleContainer(bool save = true) { m_write_particle_container = save; }

  void setContainerName(const std::string &name) { m_container_name = name; }

  void saveOutput(bool save = true) { m_save_output = save; }

  void setOutputName(const std::string &name) { m_outfile_name = name; }

  void doTruthMatching(bool truth = true) { m_truth_matching = truth; }

  void getTriggerInfo(bool get = true) { m_get_trigger_info = get; }

  void getDetectorInfo(bool detinfo = true) { m_detector_info = detinfo; }

  void getCaloInfo(bool caloinfo = true) { m_calo_info = caloinfo; }

  void requireTrackEMCalMatch(bool require = true) { m_require_track_emcal_match = require; }

  void getAllPVInfo(bool pvinfo = true) { m_get_all_PVs = pvinfo; }

  void bunchCrossingZeroOnly(bool bcZeroOnly = true) { m_bunch_crossing_zero_only = bcZeroOnly; }

  void requireBunchCrossingMatch(bool require = true) { m_require_bunch_crossing_match = require; }

  void requireTrackVertexBunchCrossingMatch(bool require = true) { m_require_track_and_vertex_match = require; }

  void usePID(bool use = true){ m_use_PID = use; }
 
  void setPIDacceptFraction(float frac = 0.2){ m_dEdx_band_width = frac; }

  /// Use alternate vertex and track fitters
  void setVertexMapNodeName(const std::string &vtx_map_node_name) { m_vtx_map_node_name = m_vtx_map_node_name_nTuple = vtx_map_node_name; }

  /// Use alternate vertex and track fitters
  void setTrackMapNodeName(const std::string &trk_map_node_name) { m_trk_map_node_name = m_trk_map_node_name_nTuple = m_origin_track_map_node_name = trk_map_node_name; }

  void magFieldFile(const std::string &fname) { m_magField = fname; }

  void getField();

  void incrementCandidateCounter(){ candidateCounter += 1; }
  void setCandidateCounter(int countNum) { candidateCounter = countNum; }
  int getCandidateCounter() { return candidateCounter; }

 private:
  bool m_has_intermediates_sPHENIX;
  bool m_constrain_to_vertex_sPHENIX;
  bool m_require_mva;
  bool m_save_dst;
  bool m_save_output;
  int candidateCounter = 0;
  std::string m_outfile_name;
  TFile *m_outfile;
  std::string m_decayDescriptor;
  std::string m_magField = "FIELDMAP_TRACKING";
};

#endif  // KFPARTICLESPHENIX_KFPARTICLESPHENIX_H
