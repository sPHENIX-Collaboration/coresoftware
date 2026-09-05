#ifndef PHGARFIELD__H
#define PHGARFIELD__H

#include <fun4all/SubsysReco.h>

#include <array>
#include <cstddef>
#include <numbers>
#include <string>
#include <vector>

#include <TRotation.h>
#include <TVector3.h>

class CDBTTree;
class PHField;
class TPolyLine3D;
class TH2;
class TH3;

namespace Garfield
{
  class ComponentUser;
  class MediumMagboltz;
}  // namespace Garfield

class PHGarfield : public SubsysReco
{
 public:
  // Backward-compatible constructor used by existing sPHENIX/Devon TPC code.
  PHGarfield(const std::string &name = "PHGarfield",
             const std::string &electricFieldMap = "",
             double spaceChargeScale_side0 = 1.0,
             double spaceChargeScale_side1 = 1.0);

  // Extended constructor for direct side-separated 3D field-map configuration.
  PHGarfield(const std::string &name,
             const std::string &electricFieldMap,
             double spaceChargeScale_side0,
             double spaceChargeScale_side1,
             const std::string &electricFieldMap3D_side0,
             const std::string &electricFieldMap3D_side1);
  ~PHGarfield() override;

  int InitRun(PHCompositeNode *) override;
  int process_event(PHCompositeNode *topNode) override;

  enum class ReverseDriftStatus
  {
    Running,
    CentralMembrane,
    RadialBoundary,
    ZBoundary,
    Stuck,
    InvalidVelocity,
    MaxSteps
  };

  ReverseDriftStatus StopHere(const double x, const double y, const double z, const double zPrevious);

  void PrintMaps() const;
  void PrintGarfield(double x, double y, double z) const;
  void PrintGasSummary() const;

  // Rigid-body placement of the TPC and magnetic-field map.
  // Units are cm and radians.  Garfield drift coordinates remain in the
  // local TPC frame; these transforms are used only to sample B(x).
  void MoveMagnet(double x_cm, double y_cm, double z_cm);
  void RotateMagnet(double theta_x, double theta_y, double theta_z);
  void MoveTpc(double x_cm, double y_cm, double z_cm);
  void RotateTpc(double theta_x, double theta_y, double theta_z);
  void SetCMVoltageDefault(double voltage) { m_CMVoltageDefault = voltage; }

  //  These are left in public namespace for easy plotting macros...
  //  The user is encouraged to add more routine to fit their analysis goals...
  // Existing macros should call this one.  Input and returned polyline are in
  // local TPC/Garfield coordinates.
  TPolyLine3D *ReverseDrift(double x_cm, double y_cm, double z_cm, double step_ns = 50.0, ReverseDriftStatus *status = nullptr);

  // Debug/visualization helper.  Input and returned polyline are in global
  // detector coordinates.  Internally the drift is still computed in local TPC
  // coordinates to keep the Garfield gas tables valid.
  TPolyLine3D *ReverseDriftGlobalCoords(double x_cm, double y_cm, double z_cm, double step_ns = 50.0, ReverseDriftStatus *status = nullptr);

  double GetRadius(size_t index) const { return radii.at(index); }

  // Axisymmetric ROOT map must contain QA/hErDefault and QA/hEzDefault.
  // The histograms are expected in cm on the axes and V/m in the bins.
  void SetElectricFieldMap(const std::string &filename) { m_electricFieldMap = filename; }

  // Side-separated 3D ROOT maps must contain Field3D/hEx, Field3D/hEy,
  // and Field3D/hEz. Axes are (r [cm], phi [rad], |z| [cm]); bin contents
  // are V/m. hEz is expressed along the +|z| solver coordinate.
  void SetElectricFieldMap3D(const std::string &side0_filename, const std::string &side1_filename)
  {
    m_electricFieldMap3D[0] = side0_filename;
    m_electricFieldMap3D[1] = side1_filename;
  }
  void SetElectricFieldMap3DSide0(const std::string &filename) { m_electricFieldMap3D[0] = filename; }
  void SetElectricFieldMap3DSide1(const std::string &filename) { m_electricFieldMap3D[1] = filename; }

  // Optional frame-charge correction maps are added on top of the existing
  // space-charge correction. The 2D format is QA/hErDefault + QA/hEzDefault;
  // the side-separated 3D format may be either Field3D/hEx + hEy + hEz
  // or the Rossegger cylindrical format hEr + hEphi + hEz at file root.
  void SetFrameElectricFieldMap(const std::string &filename) { m_frameElectricFieldMap = filename; }
  void SetFrameElectricFieldMap3D(const std::string &side0_filename, const std::string &side1_filename)
  {
    m_frameElectricFieldMap3D[0] = side0_filename;
    m_frameElectricFieldMap3D[1] = side1_filename;
  }
  void SetFrameElectricFieldMap3DSide0(const std::string &filename) { m_frameElectricFieldMap3D[0] = filename; }
  void SetFrameElectricFieldMap3DSide1(const std::string &filename) { m_frameElectricFieldMap3D[1] = filename; }
  void SetFrameChargeScale(double value)
  {
    m_frameChargeScale_side0 = value;
    m_frameChargeScale_side1 = value;
  }
  void SetFrameChargeScaleSide0(double value) { m_frameChargeScale_side0 = value; }
  void SetFrameChargeScaleSide1(double value) { m_frameChargeScale_side1 = value; }

  void SetSpaceChargeScale(double value)
  {
    m_spaceChargeScale_side0 = value;
    m_spaceChargeScale_side1 = value;
  }

  void SetSpaceChargeScaleSide0(double value) { m_spaceChargeScale_side0 = value; }
  void SetSpaceChargeScaleSide1(double value) { m_spaceChargeScale_side1 = value; }

  double GetSpaceChargeScaleSide0() const { return m_spaceChargeScale_side0; }
  double GetSpaceChargeScaleSide1() const { return m_spaceChargeScale_side1; }

  void SetZeroField(bool zerofield) { m_zerofield = zerofield; }

  // Static field-cage boundary-voltage distortions.
  // side0 = south (z < 0), side1 = north (z > 0).
  // Offsets are endpoint perturbations relative to the nominal resistor-chain
  // boundary voltage, in volts. IFC and OFC are independent, giving four
  // tunable parameters: IFC South/North and OFC South/North.
  void SetUseIFCVoltageDistortion(bool value) { m_useIFCVoltageDistortion = value; }
  void SetIFCVoltageOffset(double side0_south_v, double side1_north_v)
  {
    m_ifcVoltageOffset_side0 = side0_south_v;
    m_ifcVoltageOffset_side1 = side1_north_v;
  }
  void SetIFCVoltageOffsetSide0(double value_v) { m_ifcVoltageOffset_side0 = value_v; }
  void SetIFCVoltageOffsetSide1(double value_v) { m_ifcVoltageOffset_side1 = value_v; }
  double GetIFCVoltageOffsetSide0() const { return m_ifcVoltageOffset_side0; }
  double GetIFCVoltageOffsetSide1() const { return m_ifcVoltageOffset_side1; }

  void SetUseOFCVoltageDistortion(bool value) { m_useOFCVoltageDistortion = value; }
  void SetOFCVoltageOffset(double side0_south_v, double side1_north_v)
  {
    m_ofcVoltageOffset_side0 = side0_south_v;
    m_ofcVoltageOffset_side1 = side1_north_v;
  }
  void SetOFCVoltageOffsetSide0(double value_v) { m_ofcVoltageOffset_side0 = value_v; }
  void SetOFCVoltageOffsetSide1(double value_v) { m_ofcVoltageOffset_side1 = value_v; }
  double GetOFCVoltageOffsetSide0() const { return m_ofcVoltageOffset_side0; }
  double GetOFCVoltageOffsetSide1() const { return m_ofcVoltageOffset_side1; }

  // Convenience setter for all four fit parameters.
  void SetFieldCageVoltageOffsets(double ifc_side0_south_v, double ifc_side1_north_v,
                                  double ofc_side0_south_v, double ofc_side1_north_v)
  {
    SetIFCVoltageOffset(ifc_side0_south_v, ifc_side1_north_v);
    SetOFCVoltageOffset(ofc_side0_south_v, ofc_side1_north_v);
  }

  // Geometry of the common axisymmetric IFC/OFC boundary solver, in cm.
  void SetIFCGeometry(double inner_radius_cm, double outer_radius_cm, double half_length_cm)
  {
    m_ifcInnerRadius_cm = inner_radius_cm;
    m_ifcOuterRadius_cm = outer_radius_cm;
    m_ifcHalfLength_cm = half_length_cm;
  }
  void SetIFCVoltageModes(unsigned int modes) { m_ifcVoltageModes = modes; }
  void SetIFCGridSize(unsigned int nr, unsigned int nz)
  {
    m_ifcGridNR = nr;
    m_ifcGridNZ = nz;
  }

 private:
  void GetMagneticFieldTesla(double x_cm, double y_cm, double z_cm, double &bx_t, double &by_t, double &bz_t) const;      // Feeds magnetic field to Garfield
  void GetElectricFieldVcm(double x_cm, double y_cm, double z_cm, double &ex_vcm, double &ey_vcm, double &ez_vcm) const;  // Feeds electric field to Garfield
  void InitializeGas(const std::string &name);                                                                            // Accepts a file or a directory
  bool LoadElectricFieldCorrections(const std::string &filename);
  bool LoadElectricFieldCorrections3D(const std::string &filename, std::size_t side);
  bool HasElectricFieldCorrections3D(std::size_t side) const;
  void ClearElectricFieldCorrections3D(std::size_t side);

  bool LoadFrameElectricFieldCorrections(const std::string &filename);
  bool LoadFrameElectricFieldCorrections3D(const std::string &filename, std::size_t side);
  bool HasFrameElectricFieldCorrections3D(std::size_t side) const;
  void ClearFrameElectricFieldCorrections3D(std::size_t side);
  void AddFrameElectricFieldCorrections(double r_cm, double phi_rad, double z_cm,
                                        double &ex_vcm, double &ey_vcm, double &ez_vcm) const;
  void AddFieldCageVoltageDistortion(double x_cm, double y_cm, double z_cm,
                                     double &ex_vcm, double &ey_vcm, double &ez_vcm) const;
  bool BuildFieldCageVoltageFieldGrids();
  double InterpolateFieldCageGrid(const std::vector<double> &grid, double r_cm, double abs_z_cm) const;

  double InterpolateCorrectionVcm(const TH2 *hist, double r_cm, double abs_z_cm) const;
  double InterpolateCorrectionVcm(const TH3 *hist, double r_cm, double phi_rad, double abs_z_cm) const;
  TVector3 TpcPointToGlobalPoint(double x_cm, double y_cm, double z_cm) const;
  TVector3 GlobalPointToTpcPoint(double x_cm, double y_cm, double z_cm) const;
  TVector3 TpcPointToMagnetFieldMapPoint(double x_cm, double y_cm, double z_cm) const;
  TVector3 MagnetFieldMapVectorToTpcVector(double bx, double by, double bz) const;
  void FillRadii();
  static double bounder(double phi, double phi_min);

  CDBTTree *m_cdbTPCMAPttree{nullptr};            // Locations of the pads from CDB...
  //PHField3DCartesian *m_field{nullptr};           // The standard sPHENIX field holding container.
  PHField *m_field{nullptr};
  Garfield::ComponentUser *m_component{nullptr};  // This handles the interface of the electric and magnetic fields as handed to Garfield
  Garfield::MediumMagboltz *m_gas{nullptr};       // This is the pre-tabulated gas properties required by Garfield...
  std::string m_defaultGasfile{"/cvmfs/sphenix.sdcc.bnl.gov/files/gasfiles"};
  bool m_GasFilesLoaded{false};

  // Transform convention:
  //   global = rotation * local + translation
  // tpcrot/tpcpos place the local TPC frame in the global detector frame.
  // magrot/magpos place the magnetic-field map frame in the same global frame.
  TVector3 m_magpos{0.0, 0.0, 0.0};
  TRotation m_magrot;
  TVector3 m_tpcpos{0.0, 0.0, 0.0};
  TRotation m_tpcrot;

  // Space-charge correction maps. Histograms are cloned from the input ROOT
  // files and owned here. If a valid 3D map is loaded for a side, it takes
  // precedence over the optional axisymmetric map on that side.
  std::string m_electricFieldMap;
  std::array<std::string, 2> m_electricFieldMap3D{};
  double m_spaceChargeScale_side0{1.0};  // south, z < 0
  double m_spaceChargeScale_side1{1.0};  // north, z > 0
  double m_CMVoltageDefault{380.0};      // V/cm, nominal TPC field
  bool m_zerofield{false};
  TH2 *m_erCorrection{nullptr};          // radial correction, input bins in V/m
  TH2 *m_ezCorrection{nullptr};          // local longitudinal correction, input bins in V/m
  // Component order is Ex, Ey, Ez. Ez is along +|z| in the map.
  std::array<std::array<TH3 *, 3>, 2> m_field3DCorrection{};

  // Optional frame-charge field correction, added independently on top of the
  // ordinary space-charge field. A side-specific 3D frame map takes precedence
  // over the optional axisymmetric frame map on that side.
  std::string m_frameElectricFieldMap;
  std::array<std::string, 2> m_frameElectricFieldMap3D{};
  double m_frameChargeScale_side0{1.0};
  double m_frameChargeScale_side1{1.0};
  TH2 *m_frameErCorrection{nullptr};
  TH2 *m_frameEzCorrection{nullptr};
  std::array<std::array<TH3 *, 3>, 2> m_frameField3DCorrection{};
  // false: Ex,Ey,Ez Cartesian format; true: Er,Ephi,Ez cylindrical Rossegger format.
  std::array<bool, 2> m_frameField3DIsCylindrical{{false, false}};

  // Static IFC/OFC boundary-condition distortions (Laplace solution, not space charge).
  // The unit grids correspond to +1 V at the readout-plane endpoint of the
  // selected cage, ramping linearly from 0 V at the CM to +1 V at |z|=L.
  // IFC and OFC contributions are then superposed linearly.
  bool m_useIFCVoltageDistortion{false};
  bool m_useOFCVoltageDistortion{false};
  double m_ifcVoltageOffset_side0{0.0};  // south, z < 0, V relative to nominal
  double m_ifcVoltageOffset_side1{0.0};  // north, z > 0, V relative to nominal
  double m_ofcVoltageOffset_side0{0.0};  // south, z < 0, V relative to nominal
  double m_ofcVoltageOffset_side1{0.0};  // north, z > 0, V relative to nominal
  double m_ifcInnerRadius_cm{20.0};
  double m_ifcOuterRadius_cm{78.0};
  double m_ifcHalfLength_cm{102.325};
  unsigned int m_ifcVoltageModes{60};
  unsigned int m_ifcGridNR{241};
  unsigned int m_ifcGridNZ{321};
  double m_ifcGridDr_cm{0.0};
  double m_ifcGridDz_cm{0.0};
  bool m_fieldCageGridReady{false};
  std::vector<double> m_ifcUnitErGrid;  // V/cm for +1 V IFC endpoint perturbation
  std::vector<double> m_ifcUnitEsGrid;  // V/cm along +|z| for +1 V IFC endpoint perturbation
  std::vector<double> m_ofcUnitErGrid;  // V/cm for +1 V OFC endpoint perturbation
  std::vector<double> m_ofcUnitEsGrid;  // V/cm along +|z| for +1 V OFC endpoint perturbation

  //  These are utilities for a spot check of the overall routine:
  // std::string calibdir;
  // std::string m_DiodeContainerName;
  double PHI_MIN{-std::numbers::pi};
  std::array<double, 48> radii{};  // Radius on each layer just for test purposes...need to be cm!
};

#endif