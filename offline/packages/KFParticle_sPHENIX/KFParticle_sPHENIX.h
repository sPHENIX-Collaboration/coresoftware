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
#include <iostream>
#include <memory>     // for allocator_traits<>::valu...
#include <string>
#include <utility>  // for pair
#include <vector>   // for vector
#include <stdint.h> //include for new member added
#include <cstdint>  // fixed width integer types used for BCO counters

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
                      int numPVs, int numTracks) const;

  int parseDecayDescriptor();

  /// Parameters for the user to vary

  void setDecayDescriptor(const std::string &decayDescriptor) { m_config.decayDescriptor = decayDescriptor; }

  static const int max_particles = 99;

  //DEPRECATED: Use setDecayDescriptor instead; still called internally by parseDecayDescriptor
  void setMotherName(const std::string &mother_name)
  {
    m_config.motherName = mother_name;
  }

  //DEPRECATED: Use setDecayDescriptor instead; still called internally by parseDecayDescriptor
  void hasIntermediateStates(bool has_intermediates = true)
  {
    m_config.hasIntermediates = has_intermediates;
  }

  //DEPRECATED: Use setDecayDescriptor instead; still called internally by parseDecayDescriptor
  void setNumberOfTracks(int num_tracks)
  {
    m_config.numTracks = num_tracks;
  }

  //DEPRECATED: Use setDecayDescriptor instead; still called internally by parseDecayDescriptor
  void setNumberTracksFromIntermeditateState(const std::vector<int> &num_tracks)
  {
    m_config.numTracksFromIntermediate.insert(m_config.numTracksFromIntermediate.end(), num_tracks.begin(), num_tracks.end());
  }

  //DEPRECATED: Use setDecayDescriptor instead; still called internally by parseDecayDescriptor
  void setNumberOfIntermediateStates(int n_intermediates)
  {
    m_config.numIntermediateStates = n_intermediates;
  }

  //DEPRECATED: Use setDecayDescriptor instead; still called internally by parseDecayDescriptor
  void getChargeConjugate(bool get_charge_conjugate = true)
  {
    m_config.chargeConjugate = get_charge_conjugate;
  }

  //DEPRECATED: Use setDecayDescriptor instead; still called internally by parseDecayDescriptor
  void setDaughters(std::vector<std::pair<std::string, int>> &daughter_list)
  {
    for (const auto &entry : daughter_list)
    {
      m_config.daughterNames.push_back(entry.first);
      m_config.daughterCharges.push_back(entry.second);
    }
  }

  //DEPRECATED: Use setDecayDescriptor instead; still called internally by parseDecayDescriptor
  void setIntermediateStates(const std::vector<std::pair<std::string, int>> &intermediate_list)
  {
    for (const auto &entry : intermediate_list)
    {
      m_config.intermediateNames.push_back(entry.first);
      m_config.intermediateCharges.push_back(entry.second);
    }
  }

  void setMinimumMass(float min_mass) { m_config.motherMinMass = min_mass; }

  void setMaximumMass(float max_mass) { m_config.motherMaxMass = max_mass; }

  void setDecayTimeRange_XY(float min_decayTime, float max_decayTime)
  {
    m_config.motherMinDecayTimeXY = min_decayTime;
    m_config.motherMaxDecayTimeXY = max_decayTime;
  }

  void setDecayLengthRange_XY(float min_decayLength, float max_decayLength)
  {
    m_config.motherMinDecayLengthXY = min_decayLength;
    m_config.motherMaxDecayLengthXY = max_decayLength;
  }

  void setDecayTimeRange(float min_decayTime, float max_decayTime)
  {
    m_config.motherMinDecayTime = min_decayTime;
    m_config.motherMaxDecayTime = max_decayTime;
  }

  void setDecayLengthRange(float min_decayLength, float max_decayLength)
  {
    m_config.motherMinDecayLength = min_decayLength;
    m_config.motherMaxDecayLength = max_decayLength;
  }

  void setMinDecayTimeSignificance(float min = 0) { m_config.motherMinDecayTimeSignificance = min; }

  void setMinDecayLengthSignificance(float min = 0) { m_config.motherMinDecayLengthSignificance = min; }

  void setMinDecayLengthSignificance_XY(float min = 0) { m_config.motherMinDecayLengthSignificanceXY = min; }

  void setMinimumTrackPT(float pt) { m_config.trackMinPt = pt; }

  void setMaximumTrackPT(float pt) { m_config.trackMaxPt = pt; }

  void setMaximumTrackPTchi2(float ptchi2) { m_config.trackMaxPtChi2 = ptchi2; }

  void setMinimumTrackPV_DCA_XY(float ip) { m_config.trackMinPVDCAXY = ip; }

  void setMinimumTrackPV_DCA_StdDev_XY(float ip_significance) { m_config.trackMinPVDCAStdDevXY = ip_significance; }

  void setMinimumTrackPV_DCA(float ip) { m_config.trackMinPVDCA = ip; }

  void setMinimumTrackPV_DCA_StdDev(float ip_significance) { m_config.trackMinPVDCAStdDev = ip_significance; }

  void setMaximumTrackchi2nDOF(float trackchi2ndof) { m_config.trackMaxChi2nDOF = trackchi2ndof; }

  void setMinMVTXhits(int nHits) { m_config.minMVTXStates = nHits; } //Actually state counting but use this for backwards compatibility!

  void setMinINTThits(int nHits) { m_config.minINTTStates = nHits; } //Actually state counting but use this for backwards compatibility!

  void setMinTPChits(int nHits) { m_config.minTPCStates = nHits; } //Actually state counting but use this for backwards compatibility!

  void setMinTPOThits(int nHits) { m_config.minTPOTStates = nHits; } //Actually state counting but use this for backwards compatibility!

  void useTPConlyTracks(bool use = true) { m_config.useTPCOnlyTracks = use; }

  void setMaximumDaughterDCA_XY(float dca) { m_config.maxDaughterDCAXY = dca; }

  void setMaximumDaughterDCA(float dca) { m_config.maxDaughterDCA = dca; }
 
  void setMinimumRadialSV(float min_rad_sv) { m_config.minRadialSV = min_rad_sv; }

  void setMaximumVertexchi2nDOF(float vertexchi2nDOF) { m_config.maxVertexChi2nDOF = vertexchi2nDOF; }

  void setFlightDistancechi2(float fdchi2) { m_config.motherMinFDChi2 = fdchi2; }

  void setMinDIRA(float dira_min) { m_config.motherMinDIRA = dira_min; }

  void setMaxDIRA(float dira_max) { m_config.motherMaxDIRA = dira_max; }

  void setMinDIRA_XY(float dira_min) { m_config.motherMinDIRAXY = dira_min; }

  void setMaxDIRA_XY(float dira_max) { m_config.motherMaxDIRAXY = dira_max; }

  void setMotherPT(float mother_pt) { m_config.motherMinPt = mother_pt; }

  void setMotherPV_DCA(float mother_PV_dca) { m_config.motherMaxPVDCA = mother_PV_dca; }

  void setMotherPV_DCA_XY(float mother_PV_dca) { m_config.motherMaxPVDCAXY = mother_PV_dca; }

  void setMotherPV_DCA_StdDev(float mother_PV_dca_stddev) { m_config.motherMaxPVDCAStdDev = mother_PV_dca_stddev; }

  void setMotherPV_DCA_StdDev_XY(float mother_PV_dca_stddev) { m_config.motherMaxPVDCAStdDevXY = mother_PV_dca_stddev; }

  void setMaximumMotherVertexVolume(float vertexvol) { m_config.motherMaxVertexVolume = vertexvol; }

  void constrainToPrimaryVertex(bool constrain_to_vertex = true) { m_config.constrainToPrimaryVertex = constrain_to_vertex; }

  void useMbdVertex(bool use = true)
  {
    m_config.useMbdVertex = use;
  }

  void useCentrality(bool use = true)
  {
    m_config.useCentrality = use;
  }

  void dontUseGlobalVertex(bool dont = true) { m_config.dontUseGlobalVertex = dont; }

  void useFakePrimaryVertex(bool use_fake = true)
  {
    m_config.useFakePrimaryVertex = use_fake;
  }

  void allowZeroMassTracks(bool allow = true) { m_config.allowZeroMassTracks = allow; }

  void extrapolateTracksToSV(bool extrapolate = true)
  {
    m_config.extrapolateTracksToSV = extrapolate;
  }

  void constrainIntermediateMasses(bool constrain_int_mass = true) { m_config.constrainIntermediateMasses = constrain_int_mass; }

  void setIntermediateMassRange(const std::vector<std::pair<float, float>> &intermediate_mass_range)
  {
    m_config.intermediateMassRange.insert(m_config.intermediateMassRange.end(), intermediate_mass_range.begin(), intermediate_mass_range.end());
  }

  void setIntermediateMinPT(const std::vector<float> &intermediate_min_pt)
  {
    m_config.intermediateMinPt = intermediate_min_pt;
  }

  void setIntermediateMinPV_DCA_XY(const std::vector<float> &intermediate_min_PV_DCA)
  {
    m_config.intermediateMinPVDCAXY.insert(m_config.intermediateMinPVDCAXY.end(), intermediate_min_PV_DCA.begin(), intermediate_min_PV_DCA.end());
  }

  void setIntermediatePV_DCARange_XY(const std::vector<std::pair<float, float> /*unused*/> &intermediate_PV_DCA_range)
  {
    for (const auto &entry : intermediate_PV_DCA_range)
    {
      m_config.intermediateMinPVDCAXY.push_back(entry.first);
      m_config.intermediateMaxPVDCAXY.push_back(entry.second);
    }
  }

  void setIntermediateMinPV_DCA(const std::vector<float> &intermediate_min_PV_DCA)
  {
    m_config.intermediateMinPVDCA.insert(m_config.intermediateMinPVDCA.end(), intermediate_min_PV_DCA.begin(), intermediate_min_PV_DCA.end());
  }

  void setIntermediatePV_DCARange(const std::vector<std::pair<float, float> /*unused*/> &intermediate_PV_DCA_range)
  {
    for (const auto &entry : intermediate_PV_DCA_range)
    {
      m_config.intermediateMinPVDCA.push_back(entry.first);
      m_config.intermediateMaxPVDCA.push_back(entry.second);
    }
  }

  void setIntermediateMinPV_DCA_StdDev_XY(const std::vector<float> &intermediate_min_PV_DCA_StdDev)
  {
    m_config.intermediateMinPVDCAStdDevXY.insert(m_config.intermediateMinPVDCAStdDevXY.end(), intermediate_min_PV_DCA_StdDev.begin(), intermediate_min_PV_DCA_StdDev.end());
  }

  void setIntermediatePV_DCA_StdDevRange_XY(const std::vector<std::pair<float, float> /*unused*/> &intermediate_PV_DCA_StdDev_range)
  {
    for (const auto &entry : intermediate_PV_DCA_StdDev_range)
    {
      m_config.intermediateMinPVDCAStdDevXY.push_back(entry.first);
      m_config.intermediateMaxPVDCAStdDevXY.push_back(entry.second);
    }
  }

  void setIntermediateMinPV_DCA_StdDev(const std::vector<float> &intermediate_min_PV_DCA_StdDev)
  {
    m_config.intermediateMinPVDCAStdDev.insert(m_config.intermediateMinPVDCAStdDev.end(), intermediate_min_PV_DCA_StdDev.begin(), intermediate_min_PV_DCA_StdDev.end());
  }

  void setIntermediatePV_DCA_StdDevRange(const std::vector<std::pair<float, float> /*unused*/> &intermediate_PV_DCA_StdDev_range)
  {
    for (const auto &entry : intermediate_PV_DCA_StdDev_range)
    {
      m_config.intermediateMinPVDCAStdDev.push_back(entry.first);
      m_config.intermediateMaxPVDCAStdDev.push_back(entry.second);
    }
  }

  void setIntermediateMinDIRA(const std::vector<float> &intermediate_min_DIRA)
  {
    m_config.intermediateMinDIRA.insert(m_config.intermediateMinDIRA.end(), intermediate_min_DIRA.begin(), intermediate_min_DIRA.end());
  }

  void setIntermediateMinFDchi2(const std::vector<float> &intermediate_min_FDchi2)
  {
    m_config.intermediateMinFDChi2.insert(m_config.intermediateMinFDChi2.end(), intermediate_min_FDchi2.begin(), intermediate_min_FDchi2.end());
  }

  /// Decay length [cm] of each intermediate state, 3D, from the decay vertex of the mother: one (min, max) per intermediate
  void setIntermediateDecayLengthRange(const std::vector<std::pair<float, float>> &intermediate_decayLength_range)
  {
    m_config.intermediateDecayLengthRange = intermediate_decayLength_range;
  }

  /// Same in the transverse plane
  void setIntermediateDecayLengthRange_XY(const std::vector<std::pair<float, float>> &intermediate_decayLength_range)
  {
    m_config.intermediateDecayLengthRangeXY = intermediate_decayLength_range;
  }

  /// Decay time [ps] of each intermediate state, 3D, from the decay vertex of the mother
  void setIntermediateDecayTimeRange(const std::vector<std::pair<float, float>> &intermediate_decayTime_range)
  {
    m_config.intermediateDecayTimeRange = intermediate_decayTime_range;
  }

  /// Same in the transverse plane (pseudo-proper decay time [ps])
  void setIntermediateDecayTimeRange_XY(const std::vector<std::pair<float, float>> &intermediate_decayTime_range)
  {
    m_config.intermediateDecayTimeRangeXY = intermediate_decayTime_range;
  }

  void setIntermediateMaxVertexVolume(const std::vector<float> &intermediate_max_vertexvol)
  {
    m_config.intermediateMaxVertexVolume.insert(m_config.intermediateMaxVertexVolume.end(), intermediate_max_vertexvol.begin(), intermediate_max_vertexvol.end());
  }

  void use2Dmatching(bool use_2D_matching_tools = true) { m_config.use2DMatching = use_2D_matching_tools; }

  void saveDST(bool save = true) { m_config.saveDST = save; }

  void saveTrackContainer(bool save = true) { m_config.writeTrackContainer = save; }

  void saveParticleContainer(bool save = true) { m_config.writeParticleContainer = save; }

  void setContainerName(const std::string &name) { m_config.containerName = name; }

  void saveOutput(bool save = true) { m_config.saveOutput = save; }

  void setOutputName(const std::string &name) { m_config.outputName = name; }

  void doTruthMatching(bool truth = true) { m_config.truthMatching = truth; }

  void getTriggerInfo(bool get = true) { m_config.triggerInfo = get; }

  void getDetectorInfo(bool detinfo = true) { m_config.detectorInfo = detinfo; }

  void getCaloInfo(bool caloinfo = true) { m_config.caloInfo = caloinfo; }

  void requireTrackEMCalMatch(bool require = true) { m_config.requireTrackEMCalMatch = require; }

  void getAllPVInfo(bool pvinfo = true) { m_config.allPVInfo = pvinfo; }

  // QA output is a sub-mode of the lightweight nTuple
  void saveLightweightNTuple(bool minimal = true)
  {
    m_config.lightweightNTuple = minimal;
    if (!minimal)
    {
      m_config.lightweightNTupleToQA = false;
    }
  }
  void saveLightweightNTupleToQA(bool minimal = true) { m_config.lightweightNTuple = m_config.lightweightNTupleToQA = minimal; }  // also turns on the lightweight nTuple

  /// Functions to get detailed subsystem information
  void GetDetailedTracking(bool set_variable = true)
  {
    m_config.detailedTracking = set_variable;
    if (set_variable)
    {
      m_config.detectorInfo = true;
    }
  }
  // Function that will add additional calorimetry info to the output TTree
  void GetDetailedCalorimetry(bool set_variable = true)
  {
    m_config.detailedCalorimetry = set_variable;
    if (set_variable)
    {
      m_config.caloInfo = true;
    }
  }
  // Function that will add info from the 5x5 cell of towers around the cluster center to the output TTree
  void GetCalo5x5CellInfo(bool set_variable = true)
  {
    m_config.calo5x5CellInfo = set_variable;
    if (set_variable)
    {
      m_config.caloInfo = true;
    }
  }
  void get_dEdx_info(bool set_variable = true) { m_config.dEdxInfo = set_variable; }

  // Track-to-EMCal matching
  void setVerbosityTrack2Calo(int verb) { m_config.track2CaloVerbosity = verb; }
  void set_emcal_radius_user(float set_variable) { m_config.emcalRadius = set_variable; }
  void set_emcal_e_low_cut(float set_variable) { m_config.emcalEnergyLowCut = set_variable; }
  void set_dphi_cut_low(float set_variable) { m_config.emcalDPhiCutLow = set_variable; }
  void set_dphi_cut_high(float set_variable) { m_config.emcalDPhiCutHigh = set_variable; }
  void set_dz_cut_low(float set_variable) { m_config.emcalDZCutLow = set_variable; }
  void set_dz_cut_high(float set_variable) { m_config.emcalDZCutHigh = set_variable; }

  void bunchCrossingZeroOnly(bool bcZeroOnly = true) { m_config.bunchCrossingZeroOnly = bcZeroOnly; }

  void forceMixedEvent(bool force = true) { m_config.forceMixedEvent = force; }

  void requireBunchCrossingMatch(bool require = true) { m_config.requireBunchCrossingMatch = require; }

  /// On by default. Builds each N-prong combination (3+ body decays) once, from its (N-1)-prong subset plus a track with a
  /// larger index, instead of once for every track that can be added last. About 2.7x faster, but the DCA test is made
  /// after the vertex refit, so it can reject combinations that the exhaustive search (use = false) would keep.
  void useOrderedProngCombinations(bool use = true) { m_config.useOrderedProngCombinations = use; }

  void requireTrackVertexBunchCrossingMatch(bool require = true) { m_config.requireTrackVertexBunchCrossingMatch = require; }

  void selectMotherByMassError(bool select = true) { m_config.selectMotherByMassError = select; }

  void usePID(bool use = true){ m_config.usePID = use; }

  void useLocalPIDFile(bool use = true){ m_config.useLocalPIDFile = use; }

  void setLocalPIDFilename(const std::string &filename){ m_config.localPIDFilename = filename; }
 
  void setPIDacceptFraction(float frac = 0.2){ m_config.dEdxBandWidth = frac; }

  /// Use alternate vertex and track fitters
  void setVertexMapNodeName(const std::string &vtx_map_node_name) { m_config.vertexMapNodeName = vtx_map_node_name; }

  /// Use alternate vertex and track fitters
  void setTrackMapNodeName(const std::string &trk_map_node_name) { m_config.trackMapNodeName = trk_map_node_name; }

  void magFieldFile(const std::string &fname) { m_config.magFieldFile = fname; }

  void getField();

  void incrementCandidateCounter(){ candidateCounter += 1; }
  void setCandidateCounter(int countNum) { candidateCounter = countNum; }
  int getCandidateCounter() const { return candidateCounter; }

 private:
  KFParticle_Config m_config;  // the only writable configuration; base classes hold a const reference to it
  int candidateCounter = 0;
  //Adding member variables for BCO matching
  uint64_t m_this_event_bco{0};
  uint64_t m_last_event_bco{0};
  uint64_t m_prev_event_bco{0};
  int64_t m_prev_runNumber{-1};
  int64_t m_prev_eventNumber{-1}; //till here
  TFile *m_outfile;
};

#endif  // KFPARTICLESPHENIX_KFPARTICLESPHENIX_H
