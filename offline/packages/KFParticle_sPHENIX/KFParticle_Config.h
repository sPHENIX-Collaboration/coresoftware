#ifndef KFPARTICLESPHENIX_KFPARTICLECONFIG_H
#define KFPARTICLESPHENIX_KFPARTICLECONFIG_H

#include <limits>
#include <string>
#include <utility>
#include <vector>

// Configuration for KFParticle, as plain read-only data.
// KFParticle_sPHENIX owns the only writable instance (a private member, set through its setters).
// Every other class holds a private const reference to it, so there the fields can be read but not written.
struct KFParticle_Config
{
  int verbosity = 0;  // copy of SubsysReco::Verbosity(), taken in KFParticle_sPHENIX::Init

  // Decay definition: filled from decayDescriptor by parseDecayDescriptor in Init (or by the deprecated setters).
  std::string decayDescriptor;
  std::string motherName;
  bool hasIntermediates = false;
  int numTracks = 0;
  int numIntermediateStates = 0;
  std::vector<int> numTracksFromIntermediate;  // one entry per intermediate; the setter appends
  std::vector<std::string> intermediateNames;  // one entry per intermediate; the setter appends
  std::vector<int> intermediateCharges;        // parallel to intermediateNames
  std::vector<std::string> daughterNames;      // one entry per final-state track; the setter appends
  std::vector<int> daughterCharges;            // parallel to daughterNames
  bool chargeConjugate = false;

  // DST output
  bool saveDST = false;
  bool writeTrackContainer = true;
  bool writeParticleContainer = true;
  std::string containerName;

  // PID (dE/dx)
  bool usePID = false;
  bool useLocalPIDFile = false;
  std::string localPIDFilename;
  float dEdxBandWidth = 0.2;  // accepted fraction around the expected dE/dx

  // nTuple output
  bool saveOutput = true;
  std::string outputName = "outputData.root";
  bool lightweightNTuple = false;      // reduced set of branches
  bool lightweightNTupleToQA = false;  // reduced nTuple registered with the QA histogram manager instead of a file
  bool truthMatching = false;
  bool triggerInfo = false;
  bool detectorInfo = false;
  bool detailedTracking = true;        // residuals, layers and sensor IDs, when detectorInfo is on
  bool caloInfo = false;
  bool detailedCalorimetry = false;    // EMCal projection and cluster position, chi2, nTowers, when caloInfo is on
  bool calo5x5CellInfo = false;        // towers in the 5x5 cell around the matched EMCal cluster, when caloInfo is on
  bool requireTrackEMCalMatch = false; // only fill candidates whose tracks all match an EMCal cluster
  bool allPVInfo = false;
  bool useCentrality = false;
  bool dEdxInfo = false;

  // Input nodes and primary vertex source
  std::string trackMapNodeName = "SvtxTrackMap";  // also the map the DST output copies tracks from
  std::string vertexMapNodeName = "SvtxVertexMap";
  bool useMbdVertex = false;         // MBD vertex instead of the silicon (SVTX) vertex
  bool dontUseGlobalVertex = false;  // read the MBD or SVTX vertex map directly instead of the GlobalVertexMap

  // Track selection (KFParticle_Tools::isGoodTrack and the track-quality check before it)
  float trackMinPt = 0;
  float trackMaxPt = 5e3;
  float trackMaxPtChi2 = std::numeric_limits<float>::max();
  float trackMinPVDCAXY = -100;       // minimum over all PVs
  float trackMinPVDCAStdDevXY = -1;   // XY significance
  float trackMinPVDCA = -1;
  float trackMinPVDCAStdDev = -1;     // 3D: chi2 units (9 = 3 sigma), see audit D
  float trackMaxChi2nDOF = 100;
  int minMVTXStates = 3;              // track states, although the setters say hits
  int minINTTStates = 1;
  int minTPCStates = 20;
  int minTPOTStates = 0;
  bool useTPCOnlyTracks = false;      // skip the MVTX/INTT/TPOT minimums and reject tracks with MVTX or INTT states

  // Vertex combination: applied to every track pair and N-prong vertex (findTwoProngs, findNProngs)
  float maxDaughterDCA = 0.5;
  float maxDaughterDCAXY = 1;
  float maxVertexChi2nDOF = 20;
  float minRadialSV = -1;

  // Mother selection (buildMother, constrainToVertex); the intermediates have their own vectors
  float motherMinMass = 0;  // with the default 0-0 window nothing passes: a decay needs a mass range
  float motherMaxMass = 0;
  float motherMinPt = 0;
  float motherMaxVertexVolume = std::numeric_limits<float>::max();
  float motherMinDecayTime = -std::numeric_limits<float>::max();  // [ps]
  float motherMaxDecayTime = std::numeric_limits<float>::max();
  float motherMinDecayTimeXY = -1000;
  float motherMaxDecayTimeXY = std::numeric_limits<float>::max();
  float motherMinDecayLength = -std::numeric_limits<float>::max();  // [cm]
  float motherMaxDecayLength = std::numeric_limits<float>::max();
  float motherMinDecayLengthXY = -1000;
  float motherMaxDecayLengthXY = std::numeric_limits<float>::max();
  float motherMinDecayTimeSignificance = -1;
  float motherMinDecayLengthSignificance = -1;
  float motherMinDecayLengthSignificanceXY = -1;
  float motherMinFDChi2 = 0;
  float motherMinDIRA = -1.01;  // just outside [-1, 1], so every value passes by default
  float motherMaxDIRA = 1.01;
  float motherMinDIRAXY = -1;
  float motherMaxDIRAXY = 1;
  float motherMaxPVDCA = std::numeric_limits<float>::max();
  float motherMaxPVDCAStdDev = std::numeric_limits<float>::max();  // 3D: chi2 units (9 = 3 sigma), see audit D
  float motherMaxPVDCAXY = std::numeric_limits<float>::max();
  float motherMaxPVDCAStdDevXY = std::numeric_limits<float>::max();

  // Intermediate selection: one entry per intermediate, in decay-descriptor order.
  // "appends" = the setter adds to what is there, so calling it twice concatenates; "replaces" = assigned.
  // Read with [intermediate index] and no size check unless noted (audit B).
  std::vector<std::pair<float, float>> intermediateMassRange;  // appends
  std::vector<float> intermediateMinPt;                        // replaces
  std::vector<float> intermediateMinDIRA;                      // appends
  std::vector<float> intermediateMinFDChi2;                    // appends
  std::vector<float> intermediateMaxVertexVolume;              // appends
  // (min, max) decay length [cm] / decay time [ps] from the mother decay vertex; replaces.
  // Size-checked: an intermediate with no entry is not cut on
  std::vector<std::pair<float, float>> intermediateDecayLengthRange;
  std::vector<std::pair<float, float>> intermediateDecayLengthRangeXY;
  std::vector<std::pair<float, float>> intermediateDecayTimeRange;
  std::vector<std::pair<float, float>> intermediateDecayTimeRangeXY;
  // PV DCA windows, read by evtReco; all append. The setIntermediateMin… setters fill only the minimum,
  // so the max vectors can be shorter than the min ones (audit B)
  std::vector<float> intermediateMinPVDCA;
  std::vector<float> intermediateMaxPVDCA;
  std::vector<float> intermediateMinPVDCAStdDev;  // 3D: chi2 units (9 = 3 sigma), see audit D
  std::vector<float> intermediateMaxPVDCAStdDev;
  std::vector<float> intermediateMinPVDCAXY;
  std::vector<float> intermediateMaxPVDCAXY;
  std::vector<float> intermediateMinPVDCAStdDevXY;
  std::vector<float> intermediateMaxPVDCAStdDevXY;
  bool constrainIntermediateMasses = false;  // build the mother with each intermediate at its PDG mass instead of its fitted mass

  // Vertexing and candidate selection
  bool constrainToPrimaryVertex = false;
  bool selectMotherByMassError = true;  // best candidate = smallest mass error; false = closest to its PV (only with the PV constraint, mother only)
  bool useFakePrimaryVertex = false;   // use one fake PV (createFakePV) instead of the vertex map
  bool extrapolateTracksToSV = true;   // set each daughter's production vertex to the candidate's decay vertex
  bool allowZeroMassTracks = false;    // otherwise a daughter with mass 0 fails the candidate
  bool use2DMatching = false;          // best-candidate PV association uses the XY deviation (selectBestCombination)
  bool bunchCrossingZeroOnly = false;  // only tracks with crossing 0
  bool useOrderedProngCombinations = true;  // see the setter for the speed/acceptance trade-off
  // Bunch-crossing requirements. With forceMixedEvent, KFParticle_sPHENIX::Init switches both match flags off
  bool forceMixedEvent = false;            // keep only combinations whose tracks come from different crossings
  bool requireBunchCrossingMatch = true;   // all tracks of a combination from the same crossing
  bool requireTrackVertexBunchCrossingMatch = false;  // tracks and PV from the same crossing (with the PV constraint)

  // Magnetic field
  std::string magFieldFile = "FIELDMAP_TRACKING";  // file path or CDB tag; getField() replaces a tag with the resolved path

  // Track-EMCal matching
  int track2CaloVerbosity = 0;
  float emcalRadius = 100.70;      // radius of the track state used for the match
  float emcalEnergyLowCut = 0.2;   // minimum cluster energy
  float emcalDPhiCutLow = -0.15;
  float emcalDPhiCutHigh = 0.15;
  float emcalDZCutLow = -10;
  float emcalDZCutHigh = 10;
};

#endif  // KFPARTICLESPHENIX_KFPARTICLECONFIG_H
