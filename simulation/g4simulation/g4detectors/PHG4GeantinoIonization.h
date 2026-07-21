// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4GEANTINOIONIZATION_H
#define G4DETECTORS_PHG4GEANTINOIONIZATION_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>

class PHCompositeNode;
class PHG4HitContainer;
class PHG4TruthInfoContainer;

/**
 * Replaces the negative energy-deposition sentinel stored for charged
 * geantinos with a detector-dependent mean MIP ionization.
 *
 * This module must run after G4HIT_* nodes are loaded and before detector hit
 * reconstruction. It intentionally models only the mean energy deposition;
 * downstream detector modules retain their existing fluctuations.
 */
class PHG4GeantinoIonization : public SubsysReco
{
 public:
  explicit PHG4GeantinoIonization(const std::string& name = "PHG4GeantinoIonization");
  ~PHG4GeantinoIonization() override = default;

  int process_event(PHCompositeNode* topNode) override;

  void set_particle_name(const std::string& name) { m_particleName = name; }

  void set_mvtx_enabled(bool value) { m_detectorConfigs[0].enabled = value; }
  void set_intt_enabled(bool value) { m_detectorConfigs[1].enabled = value; }
  void set_tpc_enabled(bool value) { m_detectorConfigs[2].enabled = value; }
  void set_micromegas_enabled(bool value) { m_detectorConfigs[3].enabled = value; }
  void set_tpot_enabled(bool value) { set_micromegas_enabled(value); }

  void set_mvtx_mip_dedx(double value);
  void set_intt_mip_dedx(double value);

  /**
   * Configure the TPC mean MIP stopping power from the same gas fractions
   * passed to PHG4TpcElectronDrift. Fractions are expected to sum to one.
   */
  void set_tpc_gas_fractions(
      double neon,
      double argon,
      double cf4,
      double nitrogen,
      double isobutane);

  void set_tpot_gas_fractions(
      double neon,
      double argon,
      double cf4,
      double nitrogen,
      double isobutane);

  void set_tpot_gas_fractions(double argon, double isobutane)
  {
    set_tpot_gas_fractions(0, argon, 0, 0, isobutane);
  }

 private:
  enum class DetectorId
  {
    mvtx,
    intt,
    tpc,
    tpot
  };

  struct DetectorConfig
  {
    DetectorId detector;
    std::string name;
    std::string hitNodeName;
    bool enabled = true;
  };

  void process_detector(
      PHG4HitContainer* hits,
      const PHG4TruthInfoContainer* truthInfo,
      const DetectorConfig& config) const;

  static double mip_dedx(DetectorId detector);

  std::array<DetectorConfig, 4> m_detectorConfigs;
  std::string m_particleName = "chargedgeantino";
};

#endif  // G4DETECTORS_PHG4GEANTINOIONIZATION_H
