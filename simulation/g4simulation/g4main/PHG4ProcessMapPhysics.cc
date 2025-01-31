#include "PHG4ProcessMapPhysics.h"

#include "PHG4MCProcessDefs.h"
#include "PHG4ProcessMap.h"

// G4 process code headers
#include <Geant4/G4DecayProcessType.hh>
#include <Geant4/G4EmProcessSubType.hh>
#include <Geant4/G4FastSimulationProcessType.hh>
#include <Geant4/G4HadronicProcessType.hh>
#include <Geant4/G4OpProcessSubType.hh>
#include <Geant4/G4ProcessType.hh>
#include <Geant4/G4TransportationProcessType.hh>

PHG4ProcessMapPhysics::PHG4ProcessMapPhysics()
{
  FillMap();
}

void PHG4ProcessMapPhysics::FillMap()
{
  auto& pMap = PHG4ProcessMap::Instance();

  // clang-format off
  // G4EmProcessSubType: 1 - 26; 40; 49
  pMap.Add(fCoulombScattering, kPCoulombScattering);      // G4 value:  1
  pMap.Add(fIonisation, kPEnergyLoss);                  // G4 value:  2
  pMap.Add(fBremsstrahlung, kPBrem);                     // G4 value:  3
  pMap.Add(fPairProdByCharged, kPPair);                   // G4 value:  4
  pMap.Add(fAnnihilation, kPAnnihilation);                // G4 value:  5
  pMap.Add(fAnnihilationToMuMu, kPAnnihilation);          // G4 value:  6
    // Add code
  pMap.Add(fAnnihilationToHadrons, kPAnnihilation);       // G4 value:  7
    // Add code
  pMap.Add(fNuclearStopping, kPCoulombScattering);        // G4 value:  8
    // CHECK ??
  pMap.Add(fElectronGeneralProcess, kPNull);      // G4 value:  9

  pMap.Add(fMultipleScattering, kPMultipleScattering);    // G4 value: 10

  pMap.Add(fRayleigh, kPRayleigh);                        // G4 value: 11
  pMap.Add(fPhotoElectricEffect, kPPhotoelectric);        // G4 value: 12
  pMap.Add(fComptonScattering, kPCompton);                // G4 value: 13
  pMap.Add(fGammaConversion, kPPair);                     // G4 value: 14
  pMap.Add(fGammaConversionToMuMu, kPPair);               // G4 value: 15
    // Add code
  pMap.Add(fGammaGeneralProcess, kPNull);          // G4 value: 16
  pMap.Add(fPositronGeneralProcess, kPNull);       // G4 value: 17
  // pMap.Add(fAnnihilationToTauTau, kPAnnihilation);         // G4 value: 18
    // Add code

  pMap.Add(fCerenkov, kPCerenkov);                         // G4 value: 21
  pMap.Add(fScintillation, kPScintillation);       // G4 value: 22
  pMap.Add(fSynchrotronRadiation, kPSynchrotron);          // G4 value: 23
  pMap.Add(fTransitionRadiation, kPTransitionRadiation); // G4 value: 24

  pMap.Add(fSurfaceReflection, kPNull);            // G4 value: 25
     // low energy G4MicroElecSurface process
  // pMap.Add(fGammaReflection, kPNull);              // G4 value: 26 (added in Geant4 v11.0.0, sphenix version 10.7.4)

  // G4OpProcessSubType: 31 .36
  pMap.Add(fOpAbsorption, kPLightAbsorption);              // G4 value: 31
  pMap.Add(fOpBoundary, kPLightScattering);                // G4 value: 32
  pMap.Add(fOpRayleigh, kPRayleigh);                       // G4 value: 33
  pMap.Add(fOpWLS, kPLightWLShifting);             // G4 value: 34
  pMap.Add(fOpMieHG, kPLightScattering);           // G4 value: 35
    // Add code
  pMap.Add(fOpWLS2, kPLightWLShifting);            // G4 value: 36
    // Add code

  // G4EmProcessSubType: 40; 49
  // pMap.Add(fDarkBremsstrahlung, kPNull);           // G4 value: 40 (added in Geant4 v11.1.0, sphenix version 10.7.4)
  // pMap.Add(fMuonPairProdByCharged, kPNull);        // G4 value: 49 (added in Geant4 v11.1.0, sphenix version 10.7.4)

  // G4HadronicProcessType: 111 .161; 210; 310
  pMap.Add(fHadronElastic, kPHElastic);                    // G4 value: 111
  // pMap.Add(fNeutronGeneral, kPNull);               // G4 value: 116 (added in Geant4 v11.1.0, sphenix version 10.7.4)
  pMap.Add(fHadronInelastic, kPHInhelastic);               // G4 value: 121
  pMap.Add(fCapture, kPNCapture);                          // G4 value: 131
  pMap.Add(fMuAtomicCapture, kPHadronic);                  // G4 value: 132
    // Add code
    // was: kPNCapture, kHADR
  pMap.Add(fFission, kPNuclearFission);                    // G4 value: 141
  pMap.Add(fHadronAtRest, kPHadronic);                     // G4 value: 151
    // Add code .G4HadronStoppingProcess, "hadronCaptureAtRest"
    // was: kPNCapture, kHADR
  pMap.Add(fLeptonAtRest, kPHadronic);                     // G4 value: 152
    // No process defined with this code
  pMap.Add(fChargeExchange, kPHadronic);                   // G4 value: 161
  // pMap.Add(fNuOscillation, kPHadronic);                    // G4 value: 165 (added in Geant4 v11.2.0, sphenix version 10.7.4)
  // pMap.Add(fNuElectron, kPHadronic);                       // G4 value: 166 (added in Geant4 v11.2.0, sphenix version 10.7.4)
  // pMap.Add(fNuNucleus, kPHadronic);                        // G4 value: 167 (added in Geant4 v11.2.0, sphenix version 10.7.4)
  pMap.Add(fRadioactiveDecay, kPDecay);                    // G4 value: 210
  // pMap.Add(fEMDissociation, kPHadronic);                   // G4 value: 310 (no found in sPHENIX)
    // No process defined with this code

  // TG4HadronicProcessType: 171 .174 (VMC Hadronic process, VMC is not available in sPHENIX)
  // pMap.Add(fElectronNuclear, kPElectronNuclear);           // TG4 value: 171
  // pMap.Add(fPositronNuclear, kPPositronNuclear);           // TG4 value: 172
  // pMap.Add(fMuonNuclear, kPMuonNuclear);                   // TG4 value: 173
  // pMap.Add(fPhotoNuclear, kPPhotoNuclear);                 // TG4 value: 174

  // G4DecayProcessType: 201 .231
  pMap.Add(DECAY, kPDecay);                                // G4 value: 201
  pMap.Add(DECAY_WithSpin, kPDecay);                       // G4 value: 202
  pMap.Add(DECAY_PionMakeSpin, kPDecay);                   // G4 value: 203
  // DECAY_Radioactive ( G4 value: 210) .already added with G4HadronicProcessType
  pMap.Add(DECAY_Unknown, kPDecay);                        // G4 value: 211
  pMap.Add(DECAY_MuAtom, kPDecay);                         // G4 value: 221
  pMap.Add(DECAY_External, kPDecay);                       // G4 value: 231

  // G4TransportationProcessType: 91, 92
  pMap.Add(TRANSPORTATION, kPTransportation);         // G4 value: 91
  pMap.Add(COUPLED_TRANSPORTATION, kPTransportation); // G4 value: 92

  // G4FastSimulationProcessType: 301
  pMap.Add(FASTSIM_ManagerProcess, kPNull);        // G4 value: 301

  // G4TransportationProcessType: 401 .491
  // following processes belong to 'General' type
  pMap.Add(STEP_LIMITER, kStepMax);                // G4 value: 401
  pMap.Add(USER_SPECIAL_CUTS, kPStop);             // G4 value: 402
  pMap.Add(NEUTRON_KILLER, kPStop);                // G4 value: 403
    // was kPHadronic, kHADR
  // pMap.Add(PARALLEL_WORLD_PROCESS, kPNull);        // G4 value: 491 (added in Geant4 v11.0.0, sphenix version 10.7.4)

  // TG4BiasingProcessType: 501
  // pMap.Add(fBiasing, kPNull);                      // TG4 value: 501 (VMC Hadronic process, VMC is not available in sPHENIX)

  // TG4StackPopperProcessType: 601
  // pMap.Add(fStackPopper, kPNull);                  // TG4 value: 601 (VMC Hadronic process, VMC is not available in sPHENIX)
}

PHG4MCProcess PHG4ProcessMapPhysics::GetMCProcess(const G4VProcess* process) const
{
  return PHG4ProcessMap::Instance().GetMCProcess(process);
}
  std::string_view PHG4ProcessMapPhysics::GetMCProcessName(const G4VProcess* process) const
{
  return PHG4ProcessMap::Instance().GetMCProcessName(process);
}
