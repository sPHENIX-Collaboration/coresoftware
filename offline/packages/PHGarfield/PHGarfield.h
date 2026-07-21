#ifndef PHGARFIELD__H
#define PHGARFIELD__H

#include <fun4all/SubsysReco.h>

#include <array>
#include <numbers>
#include <string>

#include <TRotation.h>
#include <TVector3.h>

class CDBTTree;
class PHField3DCartesian;
class TPolyLine3D;
class TH2;

namespace Garfield
{
  class ComponentUser;
  class MediumMagboltz;
}  // namespace Garfield

class PHGarfield : public SubsysReco
{
 public:
  PHGarfield(const std::string &name = "PHGarfield",
             const std::string &electricFieldMap = "",
             double spaceChargeScale_side0 = 1.0,
             double spaceChargeScale_side1 = 1.0);
  ~PHGarfield() override;

  int InitRun(PHCompositeNode *) override;
  int process_event(PHCompositeNode * topNode) override;

  bool StopHere(const double x, const double y, const double z, const double zPrevious);

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
  TPolyLine3D *ReverseDrift(double x_cm, double y_cm, double z_cm, double step_ns = 50.0);

  // Debug/visualization helper.  Input and returned polyline are in global
  // detector coordinates.  Internally the drift is still computed in local TPC
  // coordinates to keep the Garfield gas tables valid.
  TPolyLine3D *ReverseDriftGlobalCoords(double x_cm, double y_cm, double z_cm, double step_ns = 50.0);

  double GetRadius(size_t index) const {return radii.at(index);}

  // ROOT map must contain QA/hErDefault and QA/hEzDefault.
  // The histograms are expected in cm on the axes and V/m in the bins.
  void SetElectricFieldMap(const std::string &filename) { m_electricFieldMap = filename; }void SetSpaceChargeScale(double value)
  {
    m_spaceChargeScale_side0 = value;
    m_spaceChargeScale_side1 = value;
  }

  void SetSpaceChargeScaleSide0(double value) { m_spaceChargeScale_side0 = value; }
  void SetSpaceChargeScaleSide1(double value) { m_spaceChargeScale_side1 = value; }

  double GetSpaceChargeScaleSide0() const { return m_spaceChargeScale_side0; }
  double GetSpaceChargeScaleSide1() const { return m_spaceChargeScale_side1; }

 private:
  void GetMagneticFieldTesla(double x_cm, double y_cm, double z_cm, double &bx_t, double &by_t, double &bz_t) const;      // Feeds magnetic field to Garfield
  void GetElectricFieldVcm(double x_cm, double y_cm, double z_cm, double &ex_vcm, double &ey_vcm, double &ez_vcm) const;  // Feeds electric field to Garfield
  void InitializeGas(const std::string &name);  // Accepts a file or a directory
  bool LoadElectricFieldCorrections(const std::string &filename);
  double InterpolateCorrectionVcm(const TH2 *hist, double r_cm, double abs_z_cm) const;
  TVector3 TpcPointToGlobalPoint(double x_cm, double y_cm, double z_cm) const;
  TVector3 GlobalPointToTpcPoint(double x_cm, double y_cm, double z_cm) const;
  TVector3 TpcPointToMagnetFieldMapPoint(double x_cm, double y_cm, double z_cm) const;
  TVector3 MagnetFieldMapVectorToTpcVector(double bx, double by, double bz) const;
  void FillRadii();
  static double bounder(double phi, double phi_min);

  CDBTTree *m_cdbTPCMAPttree{nullptr};            // Locations of the pads from CDB...
  PHField3DCartesian *m_field{nullptr};           // The standard sPHENIX field holding container.
  Garfield::ComponentUser *m_component{nullptr};  // This handles the interface of the electric and magnetic fields as handed to Garfield
  Garfield::MediumMagboltz *m_gas{nullptr};       // This is the pre-tabulated gas properties required by Garfield...
  std::string m_defaultGasfile;
  bool m_GasFilesLoaded{false};

  // Transform convention:
  //   global = rotation * local + translation
  // tpcrot/tpcpos place the local TPC frame in the global detector frame.
  // magrot/magpos place the magnetic-field map frame in the same global frame.
  TVector3 m_magpos{0.0, 0.0, 0.0};
  TRotation m_magrot;
  TVector3 m_tpcpos{0.0, 0.0, 0.0};
  TRotation m_tpcrot;

  // Axisymmetric space-charge correction maps.
  // Histograms are cloned from the input ROOT file and owned here.
  std::string m_electricFieldMap;
  double m_spaceChargeScale_side0{1.0};  // south, z < 0
  double m_spaceChargeScale_side1{1.0};  // north, z > 0
  double m_CMVoltageDefault{432.8};  // V/cm, nominal TPC field
  TH2 *m_erCorrection{nullptr};  // radial correction, input bins in V/m
  TH2 *m_ezCorrection{nullptr};  // local longitudinal correction, input bins in V/m

  //  These are utilities for a spot check of the overall routine:
  // std::string calibdir;
  // std::string m_DiodeContainerName;
  double PHI_MIN{-std::numbers::pi};
  std::array<double, 48> radii{};  // Radius on each layer just for test purposes...need to be cm!
};

#endif
