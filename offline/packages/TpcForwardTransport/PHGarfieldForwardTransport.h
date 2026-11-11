#ifndef PHGARFIELDFORWARDTRANSPORT_H
#define PHGARFIELDFORWARDTRANSPORT_H

#include <phgarfield/PHGarfield.h>

#include <array>
#include <iosfwd>
#include <string>

// Opt-in forward transport. The inherited PHGarfield reverse reconstruction
// remains unchanged; the forward integrator and pinned resource setup live here.
class PHGarfieldForwardTransport : public PHGarfield
{
 public:
  using PHGarfield::PHGarfield;
  PHGarfieldForwardTransport() = default;
  PHGarfieldForwardTransport(const PHGarfieldForwardTransport &) = delete;
  PHGarfieldForwardTransport &operator=(const PHGarfieldForwardTransport &) = delete;
  PHGarfieldForwardTransport(PHGarfieldForwardTransport &&) = delete;
  PHGarfieldForwardTransport &operator=(PHGarfieldForwardTransport &&) = delete;
  ~PHGarfieldForwardTransport() override;
  int InitRun(PHCompositeNode *topNode) override;

  static constexpr unsigned int MaxForwardDriftSteps = 200000;

  enum class ForwardDriftStatus
  {
    NotRun,
    ReachedReadout,
    NotInitialized,
    InvalidInput,
    OutsideDriftVolume,
    InvalidVelocity,
    WrongDirection,
    LeftDriftVolume,
    LeftMagneticFieldMap,
    MaxSteps
  };

  struct ForwardDriftResult
  {
    ForwardDriftStatus status{ForwardDriftStatus::NotRun};
    bool reached_readout{false};
    double x_cm{0.0};
    double y_cm{0.0};
    double z_cm{0.0};
    double drift_time_ns{0.0};
    unsigned int steps{0};
  };

  // Select the standalone, strictly pinned forward-transport initialization.
  // Set before InitRun. Ordinary PHGarfield initialization remains unchanged.
  void SetForwardTransportMode(bool value) { m_forwardTransportMode = value; }
  bool GetForwardTransportMode() const { return m_forwardTransportMode; }
  void SetTpcRotation(double theta_x, double theta_y, double theta_z);
  void SetCentralMembraneHalfWidth(double half_width_cm) { m_centralMembraneHalfWidth = half_width_cm; }
  void SetDriftVolumeRadii(double inner_radius_cm, double outer_radius_cm)
  {
    m_driftInnerRadius = inner_radius_cm;
    m_driftOuterRadius = outer_radius_cm;
  }
  void SetGasFile(const std::string &filename) { m_gasFileOverride = filename; }
  void SetLoadTpcFeeMap(bool value) { m_loadTpcFeeMap = value; }
  void SetMagneticFieldMap(const std::string &filename, double scale = 1.0)
  {
    m_magneticFieldMapOverride = filename;
    m_magneticFieldMapScale = scale;
  }
  void SetZeroMagneticField(bool value) { SetZeroField(value); }

  // Deterministic physical electron drift from an ionization point to the
  // readout plane. Input and output positions use local physical TPC-envelope
  // coordinates in cm; time and step size are in ns. readout_abs_z_cm is the
  // positive magnitude of the physical readout-plane z coordinate. Disabling
  // electric corrections supplies the nominal timing reference. Enable forward
  // transport mode before InitRun. Forward transport always stops at radial
  // field cages, including the correction-disabled reference trajectory.
  bool ForwardDrift(double x_cm,
                    double y_cm,
                    double z_cm,
                    double readout_abs_z_cm,
                    double step_ns,
                    ForwardDriftResult &result,
                    bool include_electric_field_corrections = true) const;

  bool IsReady() const
  {
    return (m_zerofield || m_field) &&
           m_gas &&
           m_GasFilesLoaded;
  }
  bool HasElectricFieldCorrections() const { return m_erCorrection && m_ezCorrection; }
  bool ElectricFieldMapCovers(double inner_radius_cm,
                              double outer_radius_cm,
                              double max_drift_distance_cm) const;
  // Check the selected space-charge maps in forward gas-surface coordinates.
  // Histogram edges define coverage; the edge-to-center clamp is intentional.
  bool ForwardElectricFieldMapsCover(double inner_radius_cm,
                                     double outer_radius_cm,
                                     double max_drift_distance_cm,
                                     std::ostream &error) const;
  bool HasKnownMagneticFieldDomain() const
  {
    return m_zerofield || m_hasMagneticFieldDomain;
  }
  bool MagneticFieldMapCoversTpcDriftVolume(
      double outer_radius_cm,
      double readout_abs_z_cm) const;
  static const char *ForwardDriftStatusName(ForwardDriftStatus status);

  double GetCentralMembraneHalfWidth() const { return m_centralMembraneHalfWidth; }
  double GetDriftInnerRadius() const { return m_driftInnerRadius; }
  double GetDriftOuterRadius() const { return m_driftOuterRadius; }
  double GetMagneticFieldMapScale() const { return m_magneticFieldMapScale; }
  bool GetZeroMagneticField() const { return m_zerofield; }
  const std::string &GetElectricFieldMap() const { return m_electricFieldMap; }
  const std::string &GetElectricFieldMap3D(std::size_t side) const { return m_electricFieldMap3D.at(side); }
  const std::string &GetFrameElectricFieldMap() const { return m_frameElectricFieldMap; }
  const std::string &GetFrameElectricFieldMap3D(std::size_t side) const { return m_frameElectricFieldMap3D.at(side); }
  double GetFrameChargeScale(std::size_t side) const { return side == 0 ? m_frameChargeScale_side0 : m_frameChargeScale_side1; }
  bool HasLoadedElectricFieldMap3D(std::size_t side) const { return HasElectricFieldCorrections3D(side); }
  bool HasLoadedFrameElectricFieldMap3D(std::size_t side) const { return HasFrameElectricFieldCorrections3D(side); }
  const std::string &GetResolvedGasFile() const { return m_resolvedGasFile; }
  const std::string &GetResolvedMagneticFieldMap() const { return m_resolvedMagneticFieldMap; }

  bool GetUseIFCVoltageDistortion() const { return m_useIFCVoltageDistortion; }
  bool GetUseOFCVoltageDistortion() const { return m_useOFCVoltageDistortion; }

 private:
  friend struct PHGarfieldForwardFieldTest;
  int InitializeForwardTransport(PHCompositeNode *topNode);
  void GetForwardElectricFieldVcm(double x_cm, double y_cm, double z_cm,
                                  double &ex_vcm, double &ey_vcm, double &ez_vcm,
                                  double correction_scale = 1.0) const;
  bool LoadMagneticFieldDomain(const std::string &filename);
  bool MagneticFieldDomainContains(double x_cm, double y_cm, double z_cm) const;

  bool m_ownsField{false};
  bool m_forwardTransportMode{false};
  std::string m_gasFileOverride;
  std::string m_resolvedGasFile;
  bool m_loadTpcFeeMap{true};
  double m_centralMembraneHalfWidth{0.0};  // cm; field-map |z| starts at the CM gas surface
  double m_driftInnerRadius{18.0};         // cm
  double m_driftOuterRadius{82.0};         // cm
  std::string m_magneticFieldMapOverride;
  std::string m_resolvedMagneticFieldMap;
  double m_magneticFieldMapScale{1.0};
  bool m_hasMagneticFieldDomain{false};
  std::array<double, 3> m_magneticFieldMinimumCm{{0.0, 0.0, 0.0}};
  std::array<double, 3> m_magneticFieldMaximumCm{{0.0, 0.0, 0.0}};
};

#endif
