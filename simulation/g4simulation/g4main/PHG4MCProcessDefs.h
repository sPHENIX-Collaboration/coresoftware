#ifndef G4MAIN_PHG4MCPROCESSDEFS_H
#define G4MAIN_PHG4MCPROCESSDEFS_H

#include <array>
#include <cstdint>
#include <string_view>

const int kMaxMCProcess = 50;

/// physics process codes (taken from VMC )
enum PHG4MCProcess : uint8_t
{
  kPPrimary = 0,  ///< Primary interaction

  kPMultipleScattering = 1,   ///< multiple scattering
  kPCoulombScattering = 45,   ///< single Coulomb scattering
  kPEnergyLoss = 2,           ///< continuous energy loss
  kPMagneticFieldL = 3,       ///< bending in mag. field
  kPDecay = 4,                ///< particle decay
  kPPair = 5,                 ///< photon pair production or
                              ///< muon direct pair production
  kPCompton = 6,              ///< Compton scattering
  kPPhotoelectric = 7,        ///< photoelectric effect
  kPBrem = 8,                 ///< bremsstrahlung
  kPDeltaRay = 9,             ///< delta-ray production
  kPAnnihilation = 10,        ///< positron annihilation
  kPAnnihilationRest = 11,    ///< positron annihilation at rest
  kPAnnihilationFlight = 12,  ///< positron annihilation in flight

  kPHadronic = 13,           ///< hadronic interaction
  kPEvaporation = 14,        ///< nuclear evaporation
  kPNuclearFission = 15,     ///< nuclear fission
  kPNuclearAbsorption = 16,  ///< nuclear absorption
  kPPbarAnnihilation = 17,   ///< antiproton annihilation
  kPNbarAnnihilation = 18,   ///< antineutron annihilation
  kPNCapture = 19,           ///< neutron capture
  kPHElastic = 20,           ///< hadronic elastic scattering
  kPHIElastic = 21,          ///< hadronic elastic incoherent scattering
  kPHCElastic = 22,          ///< hadronic elastic coherent scattering
  kPHInhelastic = 23,        ///< hadronic inelastic scattering
  kPPhotonInhelastic = 24,   ///< photon inelastic scattering

  kPMuonNuclear = 25,      ///< muon nuclear interaction
  kPElectronNuclear = 26,  ///< electron nuclear interaction
  kPPositronNuclear = 27,  ///< positron nuclear interaction
  kPPhotoNuclear = 46,     ///< photo nuclear interaction

  kPTOFlimit = 28,      ///< exceeded time of flight cut
  kPPhotoFission = 29,  ///< nuclear photofission

  kPRayleigh = 30,  ///< Rayleigh scattering

  kPNull = 31,  ///< no mechanism is active, usually at the entrance
                ///< of a new volume
  kPStop = 32,  ///< particle has fallen below energy threshold
                ///< and tracking stops

  kPLightAbsorption = 33,  ///< Cerenkov photon absorption
  kPLightDetection = 34,   ///< Optical photon detection (absorption + photoelectron production)
  kPLightScattering = 35,  ///< Cerenkov photon reflection/refraction
  kPLightWLShifting = 48,  ///< Optical photon wavelength shifting
  kStepMax = 36,           ///< step limited by STEMAX

  kPCerenkov = 37,             ///< Cerenkov photon generation
  kPFeedBackPhoton = 38,       ///< Feed back photon in RICH -- ALICE specific
  kPLightReflection = 39,      ///< Cerenkov photon reflection
  kPLightRefraction = 40,      ///< Cerenkov photon refraction or
                               /// dichroic mirror transmission
  kPSynchrotron = 41,          ///< synchrotron radiation generation
  kPScintillation = 42,        ///< scintillation
  kPTransitionRadiation = 49,  ///< transition radiation

  kPTransportation = 43,  ///< Transportation
  kPUserDefined = 47,     ///< User defined process

  kPNoProcess = 44  ///< unknown process
};

static const std::array<std::string_view, kMaxMCProcess> PHG4MCProcessName =
    {
        {"Primary particle emission",
         "Multiple scattering",
         "Energy loss",
         "Bending in magnetic field",
         "Decay",
         "Lepton pair production",
         "Compton scattering",
         "Photoelectric effect",
         "Bremstrahlung",
         "Delta ray",
         "Positron annihilation",
         "Positron annihilation at rest",
         "Positron annihilation in flight",
         "Hadronic interaction",
         "Nuclear evaporation",
         "Nuclear fission",
         "Nuclear absorbtion",
         "Antiproton annihilation",
         "Antineutron annihilation",
         "Neutron capture",
         "Hadronic elastic",
         "Hadronic incoherent elastic",
         "Hadronic coherent elastic",
         "Hadronic inelastic",
         "Photon inelastic",
         "Muon nuclear interaction",
         "Electron nuclear interaction",
         "Positron nuclear interaction",
         "Time of flight limit",
         "Nuclear photofission",
         "Rayleigh effect",
         "No active process",
         "Energy threshold",
         "Light absorption",
         "Light detection",
         "Light scattering",
         "Maximum allowed step",
         "Cerenkov production",
         "Cerenkov feed back photon",
         "Cerenkov photon reflection",
         "Cerenkov photon refraction",
         "Synchrotron radiation",
         "Scintillation",
         "Transportation",
         "Unknown process",
         "Coulomb scattering",
         "Photo nuclear interaction",
         "User defined process",
         "Optical photon wavelength shifting",
         "Transition radiation"}};

#endif  //  G4MAIN_PHG4MCPROCESSDEFS_H
