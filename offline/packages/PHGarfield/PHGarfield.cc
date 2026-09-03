#include "PHGarfield.h"
#include <cdbobjects/CDBTTree.h>
#include <phool/phool.h>

#include <phfield/PHField3DCartesian.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <TAxis.h>
#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TPolyLine3D.h>
#include <TRotation.h>
#include <TVector3.h>

#include <CLHEP/Units/SystemOfUnits.h>

#include <Garfield/ComponentUser.hh>
#include <Garfield/MediumMagboltz.hh>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>  // for basic_ostream, operat...
#include <regex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

PHGarfield::PHGarfield(const std::string& name,
                       const std::string& electricFieldMap,
                       double spaceChargeScale_side0,
                       double spaceChargeScale_side1)
  : PHGarfield(name, electricFieldMap, spaceChargeScale_side0, spaceChargeScale_side1, "", "")
{
}

PHGarfield::PHGarfield(const std::string& name,
                       const std::string& electricFieldMap,
                       double spaceChargeScale_side0,
                       double spaceChargeScale_side1,
                       const std::string& electricFieldMap3D_side0,
                       const std::string& electricFieldMap3D_side1)
  : SubsysReco(name)
  , m_electricFieldMap(electricFieldMap)
  , m_electricFieldMap3D{{electricFieldMap3D_side0, electricFieldMap3D_side1}}
  , m_spaceChargeScale_side0(spaceChargeScale_side0)
  , m_spaceChargeScale_side1(spaceChargeScale_side1)
{
}

PHGarfield::~PHGarfield()
{
  //  Housekeeping.
  delete m_field;
  delete m_cdbTPCMAPttree;
  delete m_component;
  delete m_gas;
  delete m_erCorrection;
  delete m_ezCorrection;
  ClearElectricFieldCorrections3D(0);
  ClearElectricFieldCorrections3D(1);
  delete m_frameErCorrection;
  delete m_frameEzCorrection;
  ClearFrameElectricFieldCorrections3D(0);
  ClearFrameElectricFieldCorrections3D(1);
}

int PHGarfield::InitRun(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "PHGarfield::InitRun(PHCompositeNode *topNode) Initializing" << std::endl;
  }
  CDBInterface* m_cdb = CDBInterface::instance();

  //  Here we use the CDBInterface to set up the magnetic field map:
  std::string url = m_cdb->getUrl("FIELDMAP_TRACKING");
  m_field = new PHField3DCartesian(url, 1.0);

  //  Here we use the CDBInterface to set up the channel making of the TPC:
  std::string text = m_cdb->getUrl("TPC_FEE_CHANNEL_MAP");
  m_cdbTPCMAPttree = new CDBTTree(text);
  m_cdbTPCMAPttree->LoadCalibrations();

  // Load the optional axisymmetric space-charge correction map.
  // Failure is non-fatal: Garfield then uses only the nominal 400 V/cm field.
  if (!m_electricFieldMap.empty())
  {
    if (!LoadElectricFieldCorrections(m_electricFieldMap))
    {
      std::cout << PHWHERE << " Failed to load electric-field correction map: "
                << m_electricFieldMap << std::endl;
    }
  }

  // Load optional side-separated 3D space-charge correction maps. A valid
  // 3D map takes precedence over the axisymmetric map on that side.
  for (std::size_t side = 0; side < m_electricFieldMap3D.size(); ++side)
  {
    if (!m_electricFieldMap3D[side].empty() &&
        !LoadElectricFieldCorrections3D(m_electricFieldMap3D[side], side))
    {
      std::cout << PHWHERE << " Failed to load side " << side
                << " 3D electric-field correction map: "
                << m_electricFieldMap3D[side] << std::endl;
    }
  }

  // Load optional frame-charge maps independently. They are additive and do
  // not alter the existing space-charge-map selection or k_eff flow.
  if (!m_frameElectricFieldMap.empty() &&
      !LoadFrameElectricFieldCorrections(m_frameElectricFieldMap))
  {
    std::cout << PHWHERE << " Failed to load frame electric-field correction map: "
              << m_frameElectricFieldMap << std::endl;
  }

  for (std::size_t side = 0; side < m_frameElectricFieldMap3D.size(); ++side)
  {
    if (!m_frameElectricFieldMap3D[side].empty() &&
        !LoadFrameElectricFieldCorrections3D(m_frameElectricFieldMap3D[side], side))
    {
      std::cout << PHWHERE << " Failed to load side " << side
                << " 3D frame electric-field correction map: "
                << m_frameElectricFieldMap3D[side] << std::endl;
    }
  }

  // Precompute the two axisymmetric unit boundary fields once. The IFC grid is
  // +1 V on the inner cage endpoint with the OFC held at nominal; the OFC grid
  // is +1 V on the outer cage endpoint with the IFC held at nominal. Because
  // Laplace's equation is linear, South/North IFC/OFC offsets are four
  // independent coefficients applied at field-query time.
  const bool needIFCGrid = m_useIFCVoltageDistortion &&
                           (m_ifcVoltageOffset_side0 != 0.0 || m_ifcVoltageOffset_side1 != 0.0);
  const bool needOFCGrid = m_useOFCVoltageDistortion &&
                           (m_ofcVoltageOffset_side0 != 0.0 || m_ofcVoltageOffset_side1 != 0.0);
  if ((needIFCGrid || needOFCGrid) && !BuildFieldCageVoltageFieldGrids())
  {
    std::cout << PHWHERE << " Failed to build IFC/OFC boundary-voltage field grids; "
              << "field-cage voltage corrections will be disabled." << std::endl;
    m_useIFCVoltageDistortion = false;
    m_useOFCVoltageDistortion = false;
  }

  //  Make the Garfield Component and register the methods that will interface to our fields...
  m_component = new Garfield::ComponentUser();
  m_component->SetMagneticField([this](double x, double y, double z, double& bx, double& by, double& bz)
                                { GetMagneticFieldTesla(x, y, z, bx, by, bz); });
  m_component->SetElectricField([this](double x, double y, double z, double& ex, double& ey, double& ez)
                                { GetElectricFieldVcm(x, y, z, ex, ey, ez); });

  // Here we fetch the gas from the CDB
  std::string gasfile = m_cdb->getUrl("PHGARFIELD_GAS");

  if (gasfile.empty() || !fs::exists(gasfile))
  {
    std::cout << PHWHERE << " Missing CDB gasfile: " << gasfile << std::endl;
    std::cout << PHWHERE << " Using default gasfile: " << m_defaultGasfile << std::endl;
    gasfile = m_defaultGasfile;
  }
  InitializeGas(gasfile);

  //  Diagnostic during code development...
  FillRadii();
  if (Verbosity() > 1)
  {
    PrintMaps();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHGarfield::FillRadii()
{
  //  Unload the pad map to get the radii in a handy location:
  for (unsigned int side = 0; side < 2; side++)
  {
    for (unsigned int sector = 0; sector < 12; sector++)
    {
      for (unsigned int fee = 0; fee < 26; fee++)
      {
        for (unsigned int channel = 0; channel < 256; channel++)
        {
          unsigned int key = (256 * (fee)) + channel;
          int layer = m_cdbTPCMAPttree->GetIntValue(key, "layer");
          double r = m_cdbTPCMAPttree->GetDoubleValue(key, "R") / CLHEP::cm;
          if (layer > 6)
          {
            radii[layer - 7] = r;
          }
        }
      }
    }
  }
}

void PHGarfield::PrintGarfield(double x, double y, double z) const
{
  double ex;
  double ey;
  double ez;
  double bx;
  double by;
  double bz;
  double vx;
  double vy;
  double vz;
  GetElectricFieldVcm(x, y, z, ex, ey, ez);
  GetMagneticFieldTesla(x, y, z, bx, by, bz);
  m_gas->ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  std::cout << " x:" << x
            << " y:" << y
            << " z:" << z
            << " ex:" << ex
            << " ey:" << ey
            << " ez:" << ez
            << " bx:" << bx
            << " by:" << by
            << " bz:" << bz
            << " vx:" << vx
            << " vy:" << vy
            << " vz:" << vz
            << std::endl;
}

void PHGarfield::PrintGasSummary() const
{
  if (!m_GasFilesLoaded)
  {
    std::cout << PHWHERE << "No Gas File(s) have been successfully loaded." << std::endl;
    return;
  }

  std::vector<double> nE;
  std::vector<double> nB;
  std::vector<double> nA;
  m_gas->GetFieldGrid(nE, nB, nA);

  std::cout << "Gas File Grid Dimensions: " << std::endl;
  std::cout << nE.size() << " E-fields ranging from " << nE.front() << " to " << nE.back() << std::endl;
  std::cout << nB.size() << " B-fields ranging from " << nB.front() << " to " << nB.back() << std::endl;
  std::cout << nA.size() << " Angles   ranging from " << nA.front() << " to " << nA.back() << std::endl;
}

void PHGarfield::PrintMaps() const
{
  //  Print out a few test points of the Garfield information
  PrintGarfield(0.0, 0.0, 0.1);
  PrintGarfield(0.0, 0.0, 100.0);
  PrintGarfield(0.0, 40.0, 100.1);
  PrintGarfield(0.0, 78.0, 010.1);

  //  Print out the pad coordinate map:
  int MAX = 10;
  int prints = 0;
  for (unsigned int side = 0; side < 2; side++)
  {
    for (unsigned int sector = 0; sector < 12; sector++)
    {
      for (unsigned int fee = 0; fee < 26; fee++)
      {
        for (unsigned int channel = 0; channel < 256; channel++)
        {
          unsigned int key = (256 * (fee)) + channel;
          int layer = m_cdbTPCMAPttree->GetIntValue(key, "layer");
          double phi = ((side == 1 ? 1 : -1) * (m_cdbTPCMAPttree->GetDoubleValue(key, "phi") - std::numbers::pi / 2.)) + ((sector % 12) * std::numbers::pi / 6);
          double r = m_cdbTPCMAPttree->GetDoubleValue(key, "R") / CLHEP::cm;

          phi = bounder(phi, PHI_MIN);

          if (layer > 6)
          {
            if (prints < MAX)
            {
              prints++;
              std::cout << " side: " << side;
              std::cout << " sector: " << sector;
              std::cout << " fee: " << fee;
              std::cout << " channel: " << channel;
              std::cout << " layer: " << layer;
              std::cout << " phi: " << phi;
              std::cout << " r: " << r;
              std::cout << std::endl;
            }
          }
        }
      }
    }
  }
}

void PHGarfield::MoveMagnet(double x_cm, double y_cm, double z_cm)
{
  m_magpos.SetXYZ(x_cm, y_cm, z_cm);
  if (Verbosity() > 0)
  {
    std::cout << "PHGarfield: magnetic-field map translation = ("
              << x_cm << ", " << y_cm << ", " << z_cm << ") cm" << std::endl;
  }
}

void PHGarfield::RotateMagnet(double theta_x, double theta_y, double theta_z)
{
  m_magrot.RotateX(theta_x);
  m_magrot.RotateY(theta_y);
  m_magrot.RotateZ(theta_z);
  if (Verbosity() > 0)
  {
    std::cout << "PHGarfield: magnetic-field map rotation increment = ("
              << theta_x << ", " << theta_y << ", " << theta_z << ") rad" << std::endl;
  }
}

void PHGarfield::MoveTpc(double x_cm, double y_cm, double z_cm)
{
  m_tpcpos.SetXYZ(x_cm, y_cm, z_cm);
  if (Verbosity() > 0)
  {
    std::cout << "PHGarfield: TPC translation = ("
              << x_cm << ", " << y_cm << ", " << z_cm << ") cm" << std::endl;
  }
}

void PHGarfield::RotateTpc(double theta_x, double theta_y, double theta_z)
{
  m_tpcrot.RotateX(theta_x);
  m_tpcrot.RotateY(theta_y);
  m_tpcrot.RotateZ(theta_z);
  if (Verbosity() > 0)
  {
    std::cout << "PHGarfield: TPC rotation increment = ("
              << theta_x << ", " << theta_y << ", " << theta_z << ") rad" << std::endl;
  }
}

TVector3 PHGarfield::TpcPointToGlobalPoint(double x_cm, double y_cm, double z_cm) const
{
  // Transform a point from the local TPC/Garfield frame to the global detector frame.
  return m_tpcrot * TVector3(x_cm, y_cm, z_cm) + m_tpcpos;
}

TVector3 PHGarfield::GlobalPointToTpcPoint(double x_cm, double y_cm, double z_cm) const
{
  // Transform a point from the global detector frame to the local TPC/Garfield frame.
  return m_tpcrot.Inverse() * (TVector3(x_cm, y_cm, z_cm) - m_tpcpos);
}

TVector3 PHGarfield::TpcPointToMagnetFieldMapPoint(double x_cm, double y_cm, double z_cm) const
{
  // Garfield coordinates are kept in the local TPC frame.  To sample the
  // sPHENIX magnetic-field map, first place that point into the global frame,
  // then express the global point in the magnetic-field-map frame.
  const TVector3 p_global = TpcPointToGlobalPoint(x_cm, y_cm, z_cm);
  return m_magrot.Inverse() * (p_global - m_magpos);
}

TVector3 PHGarfield::MagnetFieldMapVectorToTpcVector(double bx, double by, double bz) const
{
  // GetFieldValue returns the B vector in the magnetic-field-map axes.  Rotate
  // it to global axes and then into the local TPC axes before handing it to
  // Garfield.  No translation is applied to vectors.
  const TVector3 b_map(bx, by, bz);
  const TVector3 b_global = m_magrot * b_map;
  return m_tpcrot.Inverse() * b_global;
}

void PHGarfield::GetMagneticFieldTesla(double x_cm, double y_cm, double z_cm, double& bx_t, double& by_t, double& bz_t) const
{
  // NOTE: Garfield uses cm, V/cm, and Tesla.
  //       PHField3DCartesian/CLHEP coordinates use CLHEP length units and
  //       magnetic fields in CLHEP field units.
  //
  // Important: x_cm,y_cm,z_cm are local TPC coordinates.  We keep Garfield and
  // the gas tables in this local frame.  Only the point used to query the
  // magnetic-field map is transformed.

  const TVector3 p_map_cm = TpcPointToMagnetFieldMapPoint(x_cm, y_cm, z_cm);

  double point[4] =
      {
          p_map_cm.X() * CLHEP::cm,
          p_map_cm.Y() * CLHEP::cm,
          p_map_cm.Z() * CLHEP::cm,
          0.0};

  double bfield_map[3] = {0.0, 0.0, 0.0};

  // Get the magnetic field via the PHField3DCartesian object constructed using
  // the CDB url reference.
  if (!m_zerofield)
  {
    m_field->GetFieldValue(point, bfield_map);
  }

  const TVector3 b_tpc = MagnetFieldMapVectorToTpcVector(
      bfield_map[0], bfield_map[1], bfield_map[2]);

  bx_t = b_tpc.X() / CLHEP::tesla;
  by_t = b_tpc.Y() / CLHEP::tesla;
  bz_t = b_tpc.Z() / CLHEP::tesla;
}

void PHGarfield::GetElectricFieldVcm(double x_cm, double y_cm, double z_cm, double& ex_vcm, double& ey_vcm, double& ez_vcm) const
{
  // NOTE: Garfield uses cm, V/cm, and Tesla.
  // The correction maps use cm (and rad for phi) on their axes and V/m in
  // their bins. The maps are produced for one TPC half using s = |z|,
  // measured from the central membrane toward the readout plane.
  //
  // SetZeroField only disables B. The electric field, including space-charge
  // corrections, is intentionally identical in zero-field and field-on runs.

  const double r_cm = std::hypot(x_cm, y_cm);
  const double phi_rad = std::atan2(y_cm, x_cm);
  const double abs_z_cm = std::abs(z_cm);
  const std::size_t side = z_cm < 0.0 ? 0 : 1;

  // Nominal drift field. Electrons drift away from the central membrane,
  // therefore E_z points toward the central membrane on both sides.
  ex_vcm = 0.0;
  ey_vcm = 0.0;
  ez_vcm = z_cm > 0.0 ? -m_CMVoltageDefault : m_CMVoltageDefault;

  const double spaceChargeScale = side == 0 ? m_spaceChargeScale_side0
                                            : m_spaceChargeScale_side1;
  if (spaceChargeScale == 0.0)
  {
    AddFrameElectricFieldCorrections(r_cm, phi_rad, z_cm, ex_vcm, ey_vcm, ez_vcm);
    AddFieldCageVoltageDistortion(x_cm, y_cm, z_cm, ex_vcm, ey_vcm, ez_vcm);
    return;
  }

  // Prefer the side-specific 3D map when all three Cartesian components are
  // available. Ex and Ey are already local-TPC Cartesian components. hEz is
  // stored along the +|z| solver coordinate and must be reflected on side 0.
  if (HasElectricFieldCorrections3D(side))
  {
    ex_vcm += spaceChargeScale *
              InterpolateCorrectionVcm(m_field3DCorrection[side][0], r_cm, phi_rad, abs_z_cm);
    ey_vcm += spaceChargeScale *
              InterpolateCorrectionVcm(m_field3DCorrection[side][1], r_cm, phi_rad, abs_z_cm);

    const double delta_ez_local_vcm = spaceChargeScale *
                                      InterpolateCorrectionVcm(m_field3DCorrection[side][2], r_cm, phi_rad, abs_z_cm);
    ez_vcm += z_cm >= 0.0 ? delta_ez_local_vcm : -delta_ez_local_vcm;
    AddFrameElectricFieldCorrections(r_cm, phi_rad, z_cm, ex_vcm, ey_vcm, ez_vcm);
    AddFieldCageVoltageDistortion(x_cm, y_cm, z_cm, ex_vcm, ey_vcm, ez_vcm);
    return;
  }

  // Fall back to the legacy axisymmetric map if no 3D map is loaded for this
  // side. This preserves the existing PHGarfield behavior and API.
  if (!m_erCorrection || !m_ezCorrection)
  {
    AddFrameElectricFieldCorrections(r_cm, phi_rad, z_cm, ex_vcm, ey_vcm, ez_vcm);
    AddFieldCageVoltageDistortion(x_cm, y_cm, z_cm, ex_vcm, ey_vcm, ez_vcm);
    return;
  }

  const double delta_er_vcm = spaceChargeScale *
                              InterpolateCorrectionVcm(m_erCorrection, r_cm, abs_z_cm);
  const double delta_ez_local_vcm = spaceChargeScale *
                                    InterpolateCorrectionVcm(m_ezCorrection, r_cm, abs_z_cm);

  // Convert the cylindrical radial correction to Cartesian components.
  if (r_cm > 0.0)
  {
    ex_vcm += delta_er_vcm * x_cm / r_cm;
    ey_vcm += delta_er_vcm * y_cm / r_cm;
  }

  // hEzDefault is expressed along the local coordinate s = |z|.
  // Convert it to the local TPC Cartesian z direction.
  ez_vcm += z_cm >= 0.0 ? delta_ez_local_vcm : -delta_ez_local_vcm;

  AddFrameElectricFieldCorrections(r_cm, phi_rad, z_cm, ex_vcm, ey_vcm, ez_vcm);
  AddFieldCageVoltageDistortion(x_cm, y_cm, z_cm, ex_vcm, ey_vcm, ez_vcm);
}

// Build a fast lookup table for the axisymmetric Laplace solution generated by
// a +1 V IFC endpoint perturbation.  The field is linear in the endpoint
// voltage, so the same table is reused for independently tunable South/North
// offsets during drift propagation.
bool PHGarfield::BuildFieldCageVoltageFieldGrids()
{
  m_fieldCageGridReady = false;
  m_ifcUnitErGrid.clear();
  m_ifcUnitEsGrid.clear();
  m_ofcUnitErGrid.clear();
  m_ofcUnitEsGrid.clear();

  const double a = m_ifcInnerRadius_cm;
  const double b = m_ifcOuterRadius_cm;
  const double L = m_ifcHalfLength_cm;
  if (!(a > 0.0 && b > a && L > 0.0) ||
      m_ifcVoltageModes == 0 || m_ifcGridNR < 2 || m_ifcGridNZ < 2)
  {
    std::cout << PHWHERE << " Invalid IFC/OFC boundary solver configuration." << std::endl;
    return false;
  }

  m_ifcGridDr_cm = (b - a) / static_cast<double>(m_ifcGridNR - 1);
  m_ifcGridDz_cm = L / static_cast<double>(m_ifcGridNZ - 1);

  const std::size_t ngrid = static_cast<std::size_t>(m_ifcGridNR) * m_ifcGridNZ;
  m_ifcUnitErGrid.assign(ngrid, 0.0);
  m_ifcUnitEsGrid.assign(ngrid, 0.0);
  m_ofcUnitErGrid.assign(ngrid, 0.0);
  m_ofcUnitEsGrid.assign(ngrid, 0.0);

  // Radial eigenfunctions for the two independent Dirichlet boundaries:
  //   F_IFC(a)=1, F_IFC(b)=0
  //   F_OFC(a)=0, F_OFC(b)=1
  // and their radial derivatives.
  const std::size_t radialSize = static_cast<std::size_t>(m_ifcVoltageModes) * m_ifcGridNR;
  std::vector<double> radialIFC(radialSize, 0.0);
  std::vector<double> dradialIFC(radialSize, 0.0);
  std::vector<double> radialOFC(radialSize, 0.0);
  std::vector<double> dradialOFC(radialSize, 0.0);

  for (unsigned int n = 1; n <= m_ifcVoltageModes; ++n)
  {
    const double nd = static_cast<double>(n);
    const double k = nd * std::numbers::pi / L;
    const double i0a = std::cyl_bessel_i(0, k * a);
    const double i0b = std::cyl_bessel_i(0, k * b);
    const double k0a = std::cyl_bessel_k(0, k * a);
    const double k0b = std::cyl_bessel_k(0, k * b);

    const double denomIFC = k0b * i0a - i0b * k0a;
    const double denomOFC = k0a * i0b - i0a * k0b;

    if (!std::isfinite(denomIFC) || !std::isfinite(denomOFC) ||
        denomIFC == 0.0 || denomOFC == 0.0)
    {
      std::cout << PHWHERE << " Invalid field-cage radial mode n=" << n << std::endl;
      return false;
    }

    for (unsigned int ir = 0; ir < m_ifcGridNR; ++ir)
    {
      const double r = a + static_cast<double>(ir) * m_ifcGridDr_cm;
      const double i0r = std::cyl_bessel_i(0, k * r);
      const double k0r = std::cyl_bessel_k(0, k * r);
      const double i1r = std::cyl_bessel_i(1, k * r);
      const double k1r = std::cyl_bessel_k(1, k * r);
      const std::size_t idx = static_cast<std::size_t>(n - 1) * m_ifcGridNR + ir;

      radialIFC[idx] = (k0b * i0r - i0b * k0r) / denomIFC;
      dradialIFC[idx] = k * (k0b * i1r + i0b * k1r) / denomIFC;

      radialOFC[idx] = (k0a * i0r - i0a * k0r) / denomOFC;
      dradialOFC[idx] = k * (k0a * i1r + i0a * k1r) / denomOFC;

      if (!std::isfinite(radialIFC[idx]) || !std::isfinite(dradialIFC[idx]) ||
          !std::isfinite(radialOFC[idx]) || !std::isfinite(dradialOFC[idx]))
      {
        std::cout << PHWHERE << " Non-finite field-cage radial coefficient for mode n=" << n << std::endl;
        return false;
      }
    }
  }

  // A +1 V endpoint perturbation is taken to ramp linearly along the cage:
  // dV(s)=s/L, with s=|z|. Its sine coefficient is
  // b_n = 2*(-1)^(n+1)/(n*pi). The same longitudinal coefficients apply to
  // IFC and OFC; only the radial eigenfunction changes.
  for (unsigned int iz = 0; iz < m_ifcGridNZ; ++iz)
  {
    const double s = static_cast<double>(iz) * m_ifcGridDz_cm;
    for (unsigned int ir = 0; ir < m_ifcGridNR; ++ir)
    {
      double erIFC = 0.0;
      double esIFC = 0.0;
      double erOFC = 0.0;
      double esOFC = 0.0;

      for (unsigned int n = 1; n <= m_ifcVoltageModes; ++n)
      {
        const double nd = static_cast<double>(n);
        const double k = nd * std::numbers::pi / L;
        const double alternating = (n % 2 == 1) ? 1.0 : -1.0;
        const double bn = 2.0 * alternating / (nd * std::numbers::pi);
        const double phase = k * s;
        const std::size_t ridx = static_cast<std::size_t>(n - 1) * m_ifcGridNR + ir;
        const double sinPhase = std::sin(phase);
        const double cosPhase = std::cos(phase);

        erIFC -= bn * dradialIFC[ridx] * sinPhase;
        esIFC -= bn * radialIFC[ridx] * k * cosPhase;
        erOFC -= bn * dradialOFC[ridx] * sinPhase;
        esOFC -= bn * radialOFC[ridx] * k * cosPhase;
      }

      const std::size_t gidx = static_cast<std::size_t>(iz) * m_ifcGridNR + ir;
      m_ifcUnitErGrid[gidx] = erIFC;
      m_ifcUnitEsGrid[gidx] = esIFC;
      m_ofcUnitErGrid[gidx] = erOFC;
      m_ofcUnitEsGrid[gidx] = esOFC;
    }
  }

  m_fieldCageGridReady = true;
  if (Verbosity() > 0)
  {
    std::cout << "PHGarfield: built IFC/OFC boundary field grids: "
              << m_ifcGridNR << " x " << m_ifcGridNZ
              << ", modes=" << m_ifcVoltageModes
              << ", r=[" << a << ", " << b << "] cm"
              << ", |z|=[0, " << L << "] cm"
              << "\n  IFC dV endpoint South/North = ("
              << m_ifcVoltageOffset_side0 << ", "
              << m_ifcVoltageOffset_side1 << ") V"
              << "\n  OFC dV endpoint South/North = ("
              << m_ofcVoltageOffset_side0 << ", "
              << m_ofcVoltageOffset_side1 << ") V" << std::endl;
  }
  return true;
}

double PHGarfield::InterpolateFieldCageGrid(const std::vector<double>& grid,
                                            const double r_cm,
                                            const double abs_z_cm) const
{
  if (!m_fieldCageGridReady || grid.empty() || m_ifcGridNR < 2 || m_ifcGridNZ < 2 ||
      m_ifcGridDr_cm <= 0.0 || m_ifcGridDz_cm <= 0.0)
  {
    return 0.0;
  }

  const double r = std::clamp(r_cm, m_ifcInnerRadius_cm, m_ifcOuterRadius_cm);
  const double s = std::clamp(std::abs(abs_z_cm), 0.0, m_ifcHalfLength_cm);

  const double fr = (r - m_ifcInnerRadius_cm) / m_ifcGridDr_cm;
  const double fz = s / m_ifcGridDz_cm;

  unsigned int ir0 = static_cast<unsigned int>(std::floor(fr));
  unsigned int iz0 = static_cast<unsigned int>(std::floor(fz));
  if (ir0 >= m_ifcGridNR - 1)
  {
    ir0 = m_ifcGridNR - 2;
  }
  if (iz0 >= m_ifcGridNZ - 1)
  {
    iz0 = m_ifcGridNZ - 2;
  }

  const unsigned int ir1 = ir0 + 1;
  const unsigned int iz1 = iz0 + 1;
  const double tr = std::clamp(fr - static_cast<double>(ir0), 0.0, 1.0);
  const double tz = std::clamp(fz - static_cast<double>(iz0), 0.0, 1.0);

  const auto at = [&](const unsigned int ir, const unsigned int iz)
  {
    return grid[static_cast<std::size_t>(iz) * m_ifcGridNR + ir];
  };

  const double v00 = at(ir0, iz0);
  const double v10 = at(ir1, iz0);
  const double v01 = at(ir0, iz1);
  const double v11 = at(ir1, iz1);
  const double v0 = v00 + tr * (v10 - v00);
  const double v1 = v01 + tr * (v11 - v01);
  return v0 + tz * (v1 - v0);
}

// Add the static fields from changed IFC and/or OFC boundary conditions.
// These are independent of the Rossegger volume-space-charge correction and
// are superposed linearly. For each side there are therefore two independent
// parameters: dV_IFC and dV_OFC.
void PHGarfield::AddFieldCageVoltageDistortion(const double x_cm,
                                               const double y_cm,
                                               const double z_cm,
                                               double& ex_vcm,
                                               double& ey_vcm,
                                               double& ez_vcm) const
{
  if (!m_fieldCageGridReady ||
      (!m_useIFCVoltageDistortion && !m_useOFCVoltageDistortion))
  {
    return;
  }

  const std::size_t side = z_cm < 0.0 ? 0 : 1;
  double deltaVIFC = 0.0;
  if (m_useIFCVoltageDistortion)
  {
    deltaVIFC = (side == 0) ? m_ifcVoltageOffset_side0 : m_ifcVoltageOffset_side1;
  }

  double deltaVOFC = 0.0;
  if (m_useOFCVoltageDistortion)
  {
    deltaVOFC = (side == 0) ? m_ofcVoltageOffset_side0 : m_ofcVoltageOffset_side1;
  }

  if (deltaVIFC == 0.0 && deltaVOFC == 0.0)
  {
    return;
  }

  const double r = std::hypot(x_cm, y_cm);
  const double s = std::abs(z_cm);
  if (r < m_ifcInnerRadius_cm || r > m_ifcOuterRadius_cm ||
      s > m_ifcHalfLength_cm)
  {
    return;
  }

  const double delta_er_vcm =
      deltaVIFC * InterpolateFieldCageGrid(m_ifcUnitErGrid, r, s) +
      deltaVOFC * InterpolateFieldCageGrid(m_ofcUnitErGrid, r, s);
  const double delta_es_vcm =
      deltaVIFC * InterpolateFieldCageGrid(m_ifcUnitEsGrid, r, s) +
      deltaVOFC * InterpolateFieldCageGrid(m_ofcUnitEsGrid, r, s);

  if (r > 0.0)
  {
    ex_vcm += delta_er_vcm * x_cm / r;
    ey_vcm += delta_er_vcm * y_cm / r;
  }

  // s=|z|; reflect the local +s field on the south side.
  ez_vcm += z_cm >= 0.0 ? delta_es_vcm : -delta_es_vcm;
}

bool PHGarfield::LoadElectricFieldCorrections(const std::string& filename)
{
  std::unique_ptr<TFile> input(TFile::Open(filename.c_str(), "READ"));
  if (!input || input->IsZombie())
  {
    std::cout << PHWHERE << " Could not open electric-field map: "
              << filename << std::endl;
    return false;
  }

  auto* er = dynamic_cast<TH2*>(input->Get("QA/hErDefault"));
  auto* ez = dynamic_cast<TH2*>(input->Get("QA/hEzDefault"));

  // Also allow maps written at the ROOT-file top level.
  if (!er)
  {
    er = dynamic_cast<TH2*>(input->Get("hErDefault"));
  }
  if (!ez)
  {
    ez = dynamic_cast<TH2*>(input->Get("hEzDefault"));
  }

  if (!er || !ez)
  {
    std::cout << PHWHERE
              << " Missing QA/hErDefault or QA/hEzDefault in "
              << filename << std::endl;
    return false;
  }

  delete m_erCorrection;
  delete m_ezCorrection;
  m_erCorrection = dynamic_cast<TH2*>(er->Clone("PHGarfield_ErCorrection"));
  m_ezCorrection = dynamic_cast<TH2*>(ez->Clone("PHGarfield_EzCorrection"));

  if (!m_erCorrection || !m_ezCorrection)
  {
    delete m_erCorrection;
    delete m_ezCorrection;
    m_erCorrection = nullptr;
    m_ezCorrection = nullptr;
    return false;
  }

  m_erCorrection->SetDirectory(nullptr);
  m_ezCorrection->SetDirectory(nullptr);

  std::cout << "Loaded axisymmetric electric-field corrections from "
            << filename << std::endl;
  std::cout << "  scale k_eff side0/south/z<0 = "
            << m_spaceChargeScale_side0 << std::endl;
  std::cout << "  scale k_eff side1/north/z>0 = "
            << m_spaceChargeScale_side1 << std::endl;
  std::cout << "  r range [cm] = ["
            << m_erCorrection->GetXaxis()->GetXmin() << ", "
            << m_erCorrection->GetXaxis()->GetXmax() << "]" << std::endl;
  std::cout << "  |z| range [cm] = ["
            << m_erCorrection->GetYaxis()->GetXmin() << ", "
            << m_erCorrection->GetYaxis()->GetXmax() << "]" << std::endl;

  return true;
}

double PHGarfield::InterpolateCorrectionVcm(
    const TH2* hist,
    const double r_cm,
    const double abs_z_cm) const
{
  if (!hist)
  {
    return 0.0;
  }

  const TAxis* xaxis = hist->GetXaxis();
  const TAxis* yaxis = hist->GetYaxis();

  if (!xaxis || !yaxis)
  {
    return 0.0;
  }

  // Interpolate safely between bin centers, not at histogram edges.
  const double r_min = xaxis->GetBinCenter(1);
  const double r_max = xaxis->GetBinCenter(xaxis->GetNbins());

  const double z_min = yaxis->GetBinCenter(1);
  const double z_max = yaxis->GetBinCenter(yaxis->GetNbins());

  constexpr double epsilon = 1.0e-6;

  const double r_eval =
      std::clamp(r_cm, r_min + epsilon, r_max - epsilon);

  const double z_eval =
      std::clamp(std::abs(abs_z_cm),
                 z_min + epsilon,
                 z_max - epsilon);

  // Input histogram is in V/m. Garfield expects V/cm.
  return 0.01 * hist->Interpolate(r_eval, z_eval);
}

bool PHGarfield::LoadElectricFieldCorrections3D(const std::string& filename, const std::size_t side)
{
  if (side >= m_field3DCorrection.size())
  {
    std::cout << PHWHERE << " Invalid TPC side " << side
              << " for 3D electric-field map " << filename << std::endl;
    return false;
  }

  std::unique_ptr<TFile> input(TFile::Open(filename.c_str(), "READ"));
  if (!input || input->IsZombie())
  {
    std::cout << PHWHERE << " Could not open 3D electric-field map: "
              << filename << std::endl;
    return false;
  }

  std::array<TH3*, 3> source{{dynamic_cast<TH3*>(input->Get("Field3D/hEx")),
                              dynamic_cast<TH3*>(input->Get("Field3D/hEy")),
                              dynamic_cast<TH3*>(input->Get("Field3D/hEz"))}};

  if (!source[0] || !source[1] || !source[2])
  {
    std::cout << PHWHERE
              << " Missing Field3D/hEx, Field3D/hEy, or Field3D/hEz in "
              << filename << std::endl;
    return false;
  }

  const auto sameAxis = [](const TAxis* lhs, const TAxis* rhs)
  {
    if (!lhs || !rhs || lhs->GetNbins() != rhs->GetNbins())
    {
      return false;
    }
    constexpr double tolerance = 1.0e-9;
    return std::abs(lhs->GetXmin() - rhs->GetXmin()) < tolerance &&
           std::abs(lhs->GetXmax() - rhs->GetXmax()) < tolerance;
  };

  for (std::size_t component = 1; component < source.size(); ++component)
  {
    if (!sameAxis(source[0]->GetXaxis(), source[component]->GetXaxis()) ||
        !sameAxis(source[0]->GetYaxis(), source[component]->GetYaxis()) ||
        !sameAxis(source[0]->GetZaxis(), source[component]->GetZaxis()))
    {
      std::cout << PHWHERE << " Inconsistent hEx/hEy/hEz axes in "
                << filename << std::endl;
      return false;
    }
  }

  const TAxis* phiAxis = source[0]->GetYaxis();
  constexpr double phiTolerance = 1.0e-3;
  if (std::abs((phiAxis->GetXmax() - phiAxis->GetXmin()) - 2.0 * std::numbers::pi) > phiTolerance)
  {
    std::cout << PHWHERE << " 3D field-map phi axis must span 2*pi radians in "
              << filename << std::endl;
    return false;
  }

  std::array<TH3*, 3> cloned{};
  constexpr std::array<const char*, 3> componentNames{{"Ex", "Ey", "Ez"}};
  for (std::size_t component = 0; component < cloned.size(); ++component)
  {
    std::ostringstream cloneName;
    cloneName << "PHGarfield_" << componentNames[component]
              << "Correction_side" << side;
    cloned[component] = dynamic_cast<TH3*>(source[component]->Clone(cloneName.str().c_str()));
    if (!cloned[component])
    {
      for (TH3* histogram : cloned)
      {
        delete histogram;
      }
      return false;
    }
    cloned[component]->SetDirectory(nullptr);
  }

  ClearElectricFieldCorrections3D(side);
  m_field3DCorrection[side] = cloned;

  std::cout << "Loaded 3D electric-field corrections for side " << side
            << (side == 0 ? " (south/z<0)" : " (north/z>0)")
            << " from " << filename << std::endl;
  std::cout << "  scale k_eff = "
            << (side == 0 ? m_spaceChargeScale_side0 : m_spaceChargeScale_side1)
            << std::endl;
  std::cout << "  r range [cm] = ["
            << source[0]->GetXaxis()->GetXmin() << ", "
            << source[0]->GetXaxis()->GetXmax() << "]" << std::endl;
  std::cout << "  phi range [rad] = ["
            << source[0]->GetYaxis()->GetXmin() << ", "
            << source[0]->GetYaxis()->GetXmax() << "]" << std::endl;
  std::cout << "  |z| range [cm] = ["
            << source[0]->GetZaxis()->GetXmin() << ", "
            << source[0]->GetZaxis()->GetXmax() << "]" << std::endl;

  return true;
}

bool PHGarfield::HasElectricFieldCorrections3D(const std::size_t side) const
{
  return side < m_field3DCorrection.size() &&
         m_field3DCorrection[side][0] &&
         m_field3DCorrection[side][1] &&
         m_field3DCorrection[side][2];
}

void PHGarfield::ClearElectricFieldCorrections3D(const std::size_t side)
{
  if (side >= m_field3DCorrection.size())
  {
    return;
  }

  for (TH3*& histogram : m_field3DCorrection[side])
  {
    delete histogram;
    histogram = nullptr;
  }
}

bool PHGarfield::LoadFrameElectricFieldCorrections(const std::string& filename)
{
  std::unique_ptr<TFile> input(TFile::Open(filename.c_str(), "READ"));
  if (!input || input->IsZombie())
  {
    std::cout << PHWHERE << " Could not open frame electric-field map: "
              << filename << std::endl;
    return false;
  }

  auto* er = dynamic_cast<TH2*>(input->Get("QA/hErDefault"));
  auto* ez = dynamic_cast<TH2*>(input->Get("QA/hEzDefault"));
  if (!er)
  {
    er = dynamic_cast<TH2*>(input->Get("hErDefault"));
  }
  if (!ez)
  {
    ez = dynamic_cast<TH2*>(input->Get("hEzDefault"));
  }

  if (!er || !ez)
  {
    std::cout << PHWHERE
              << " Missing QA/hErDefault or QA/hEzDefault in frame map "
              << filename << std::endl;
    return false;
  }

  TH2* erClone = dynamic_cast<TH2*>(er->Clone("PHGarfield_FrameErCorrection"));
  TH2* ezClone = dynamic_cast<TH2*>(ez->Clone("PHGarfield_FrameEzCorrection"));
  if (!erClone || !ezClone)
  {
    delete erClone;
    delete ezClone;
    return false;
  }
  erClone->SetDirectory(nullptr);
  ezClone->SetDirectory(nullptr);

  delete m_frameErCorrection;
  delete m_frameEzCorrection;
  m_frameErCorrection = erClone;
  m_frameEzCorrection = ezClone;

  std::cout << "Loaded axisymmetric frame electric-field corrections from "
            << filename << std::endl;
  std::cout << "  frame k_eff side0/side1 = "
            << m_frameChargeScale_side0 << "/" << m_frameChargeScale_side1
            << std::endl;
  return true;
}

bool PHGarfield::LoadFrameElectricFieldCorrections3D(const std::string& filename, const std::size_t side)
{
  if (side >= m_frameField3DCorrection.size())
  {
    std::cout << PHWHERE << " Invalid TPC side " << side
              << " for 3D frame electric-field map " << filename << std::endl;
    return false;
  }

  std::unique_ptr<TFile> input(TFile::Open(filename.c_str(), "READ"));
  if (!input || input->IsZombie())
  {
    std::cout << PHWHERE << " Could not open 3D frame electric-field map: "
              << filename << std::endl;
    return false;
  }

  // Accept both PHGarfield Cartesian maps and the native Rossegger 3D output.
  // Rossegger writes hEr, hEphi, hEz at the ROOT-file top level.
  std::array<TH3*, 3> source{{dynamic_cast<TH3*>(input->Get("Field3D/hEx")),
                              dynamic_cast<TH3*>(input->Get("Field3D/hEy")),
                              dynamic_cast<TH3*>(input->Get("Field3D/hEz"))}};
  bool cylindrical = false;

  if (!source[0] || !source[1] || !source[2])
  {
    source = {{dynamic_cast<TH3*>(input->Get("hEr")),
               dynamic_cast<TH3*>(input->Get("hEphi")),
               dynamic_cast<TH3*>(input->Get("hEz"))}};
    cylindrical = true;
  }

  if (!source[0] || !source[1] || !source[2])
  {
    std::cout << PHWHERE
              << " Missing either Field3D/hEx,hEy,hEz or hEr,hEphi,hEz in frame map "
              << filename << std::endl;
    return false;
  }

  const auto sameAxis = [](const TAxis* lhs, const TAxis* rhs)
  {
    if (!lhs || !rhs || lhs->GetNbins() != rhs->GetNbins())
    {
      return false;
    }
    constexpr double tolerance = 1.0e-9;
    return std::abs(lhs->GetXmin() - rhs->GetXmin()) < tolerance &&
           std::abs(lhs->GetXmax() - rhs->GetXmax()) < tolerance;
  };

  for (std::size_t component = 1; component < source.size(); ++component)
  {
    if (!sameAxis(source[0]->GetXaxis(), source[component]->GetXaxis()) ||
        !sameAxis(source[0]->GetYaxis(), source[component]->GetYaxis()) ||
        !sameAxis(source[0]->GetZaxis(), source[component]->GetZaxis()))
    {
      std::cout << PHWHERE << " Inconsistent 3D frame-map axes in "
                << filename << std::endl;
      return false;
    }
  }

  const TAxis* phiAxis = source[0]->GetYaxis();
  constexpr double phiTolerance = 1.0e-3;
  if (std::abs((phiAxis->GetXmax() - phiAxis->GetXmin()) - 2.0 * std::numbers::pi) > phiTolerance)
  {
    std::cout << PHWHERE << " 3D frame-map phi axis must span 2*pi radians in "
              << filename << std::endl;
    return false;
  }

  std::array<TH3*, 3> cloned{};
  const std::array<const char*, 3> componentNames =
      cylindrical ? std::array<const char*, 3>{{"Er", "Ephi", "Ez"}}
                  : std::array<const char*, 3>{{"Ex", "Ey", "Ez"}};
  for (std::size_t component = 0; component < cloned.size(); ++component)
  {
    std::ostringstream cloneName;
    cloneName << "PHGarfield_Frame" << componentNames[component]
              << "Correction_side" << side;
    cloned[component] = dynamic_cast<TH3*>(source[component]->Clone(cloneName.str().c_str()));
    if (!cloned[component])
    {
      for (TH3* histogram : cloned)
      {
        delete histogram;
      }
      return false;
    }
    cloned[component]->SetDirectory(nullptr);
  }

  ClearFrameElectricFieldCorrections3D(side);
  m_frameField3DCorrection[side] = cloned;
  m_frameField3DIsCylindrical[side] = cylindrical;

  std::cout << "Loaded 3D frame electric-field corrections for side " << side
            << (side == 0 ? " (south/z<0)" : " (north/z>0)")
            << " from " << filename << std::endl;
  std::cout << "  format = " << (cylindrical ? "Rossegger Er/Ephi/Ez" : "Cartesian Ex/Ey/Ez")
            << ", frame k_eff = "
            << (side == 0 ? m_frameChargeScale_side0 : m_frameChargeScale_side1)
            << std::endl;
  return true;
}

bool PHGarfield::HasFrameElectricFieldCorrections3D(const std::size_t side) const
{
  return side < m_frameField3DCorrection.size() &&
         m_frameField3DCorrection[side][0] &&
         m_frameField3DCorrection[side][1] &&
         m_frameField3DCorrection[side][2];
}

void PHGarfield::ClearFrameElectricFieldCorrections3D(const std::size_t side)
{
  if (side >= m_frameField3DCorrection.size())
  {
    return;
  }
  for (TH3*& histogram : m_frameField3DCorrection[side])
  {
    delete histogram;
    histogram = nullptr;
  }
  m_frameField3DIsCylindrical[side] = false;
}

void PHGarfield::AddFrameElectricFieldCorrections(const double r_cm,
                                                  const double phi_rad,
                                                  const double z_cm,
                                                  double& ex_vcm,
                                                  double& ey_vcm,
                                                  double& ez_vcm) const
{
  const std::size_t side = z_cm < 0.0 ? 0 : 1;
  const double frameScale = side == 0 ? m_frameChargeScale_side0
                                      : m_frameChargeScale_side1;
  if (frameScale == 0.0)
  {
    return;
  }

  const double abs_z_cm = std::abs(z_cm);

  // As for the ordinary space-charge field, prefer a valid side-specific 3D
  // frame map over the optional axisymmetric frame map.
  if (HasFrameElectricFieldCorrections3D(side))
  {
    if (m_frameField3DIsCylindrical[side])
    {
      const double deltaErVcm = frameScale *
                                InterpolateCorrectionVcm(m_frameField3DCorrection[side][0], r_cm, phi_rad, abs_z_cm);
      const double deltaEphiVcm = frameScale *
                                  InterpolateCorrectionVcm(m_frameField3DCorrection[side][1], r_cm, phi_rad, abs_z_cm);
      ex_vcm += deltaErVcm * std::cos(phi_rad) - deltaEphiVcm * std::sin(phi_rad);
      ey_vcm += deltaErVcm * std::sin(phi_rad) + deltaEphiVcm * std::cos(phi_rad);
    }
    else
    {
      ex_vcm += frameScale *
                InterpolateCorrectionVcm(m_frameField3DCorrection[side][0], r_cm, phi_rad, abs_z_cm);
      ey_vcm += frameScale *
                InterpolateCorrectionVcm(m_frameField3DCorrection[side][1], r_cm, phi_rad, abs_z_cm);
    }

    const double deltaEzLocalVcm = frameScale *
                                   InterpolateCorrectionVcm(m_frameField3DCorrection[side][2], r_cm, phi_rad, abs_z_cm);
    ez_vcm += z_cm >= 0.0 ? deltaEzLocalVcm : -deltaEzLocalVcm;
    return;
  }

  if (!m_frameErCorrection || !m_frameEzCorrection)
  {
    return;
  }

  const double deltaErVcm = frameScale *
                            InterpolateCorrectionVcm(m_frameErCorrection, r_cm, abs_z_cm);
  const double deltaEzLocalVcm = frameScale *
                                 InterpolateCorrectionVcm(m_frameEzCorrection, r_cm, abs_z_cm);

  if (r_cm > 0.0)
  {
    ex_vcm += deltaErVcm * std::cos(phi_rad);
    ey_vcm += deltaErVcm * std::sin(phi_rad);
  }
  ez_vcm += z_cm >= 0.0 ? deltaEzLocalVcm : -deltaEzLocalVcm;
}

double PHGarfield::InterpolateCorrectionVcm(
    const TH3* hist,
    const double r_cm,
    const double phi_rad,
    const double abs_z_cm) const
{
  if (!hist)
  {
    return 0.0;
  }

  const TAxis* rAxis = hist->GetXaxis();
  const TAxis* phiAxis = hist->GetYaxis();
  const TAxis* zAxis = hist->GetZaxis();
  if (!rAxis || !phiAxis || !zAxis)
  {
    return 0.0;
  }

  struct AxisWeights
  {
    int low{1};
    int high{1};
    double highWeight{0.0};
  };

  const auto clampedWeights = [](const TAxis* axis, const double value)
  {
    AxisWeights result;
    const int nBins = axis->GetNbins();
    if (nBins <= 1)
    {
      return result;
    }

    const double firstCenter = axis->GetBinCenter(1);
    const double lastCenter = axis->GetBinCenter(nBins);
    if (value <= firstCenter)
    {
      result.low = 1;
      result.high = 1;
      return result;
    }
    if (value >= lastCenter)
    {
      result.low = nBins;
      result.high = nBins;
      return result;
    }

    int bin = axis->FindFixBin(value);
    bin = std::clamp(bin, 1, nBins);
    const double center = axis->GetBinCenter(bin);
    if (value >= center)
    {
      result.low = bin;
      result.high = std::min(bin + 1, nBins);
    }
    else
    {
      result.low = std::max(bin - 1, 1);
      result.high = bin;
    }

    const double lowCenter = axis->GetBinCenter(result.low);
    const double highCenter = axis->GetBinCenter(result.high);
    if (result.low != result.high)
    {
      result.highWeight = (value - lowCenter) / (highCenter - lowCenter);
    }
    return result;
  };

  const auto periodicWeights = [](const TAxis* axis, const double value)
  {
    AxisWeights result;
    const int nBins = axis->GetNbins();
    if (nBins <= 1)
    {
      return result;
    }

    const double minimum = axis->GetXmin();
    const double maximum = axis->GetXmax();
    const double period = maximum - minimum;
    if (!(period > 0.0))
    {
      return result;
    }

    double wrapped = std::fmod(value - minimum, period);
    if (wrapped < 0.0)
    {
      wrapped += period;
    }
    wrapped += minimum;

    int bin = axis->FindFixBin(wrapped);
    bin = std::clamp(bin, 1, nBins);
    const double center = axis->GetBinCenter(bin);

    double lowCenter = 0.0;
    double highCenter = 0.0;
    if (wrapped >= center)
    {
      result.low = bin;
      result.high = bin == nBins ? 1 : bin + 1;
      lowCenter = center;
      highCenter = result.high == 1 ? axis->GetBinCenter(1) + period
                                    : axis->GetBinCenter(result.high);
      if (result.high == 1 && wrapped < lowCenter)
      {
        wrapped += period;
      }
    }
    else
    {
      result.high = bin;
      result.low = bin == 1 ? nBins : bin - 1;
      highCenter = center;
      lowCenter = result.low == nBins ? axis->GetBinCenter(nBins) - period
                                      : axis->GetBinCenter(result.low);
    }

    result.highWeight = (wrapped - lowCenter) / (highCenter - lowCenter);
    result.highWeight = std::clamp(result.highWeight, 0.0, 1.0);
    return result;
  };

  const AxisWeights r = clampedWeights(rAxis, r_cm);
  const AxisWeights phi = periodicWeights(phiAxis, phi_rad);
  const AxisWeights z = clampedWeights(zAxis, std::abs(abs_z_cm));

  const std::array<int, 2> rBins{{r.low, r.high}};
  const std::array<int, 2> phiBins{{phi.low, phi.high}};
  const std::array<int, 2> zBins{{z.low, z.high}};
  const std::array<double, 2> rWeights{{1.0 - r.highWeight, r.highWeight}};
  const std::array<double, 2> phiWeights{{1.0 - phi.highWeight, phi.highWeight}};
  const std::array<double, 2> zWeights{{1.0 - z.highWeight, z.highWeight}};

  double value_vpm = 0.0;
  for (std::size_t ir = 0; ir < rBins.size(); ++ir)
  {
    for (std::size_t iphi = 0; iphi < phiBins.size(); ++iphi)
    {
      for (std::size_t iz = 0; iz < zBins.size(); ++iz)
      {
        value_vpm += rWeights[ir] * phiWeights[iphi] * zWeights[iz] *
                     hist->GetBinContent(rBins[ir], phiBins[iphi], zBins[iz]);
      }
    }
  }

  // Input histogram is in V/m. Garfield expects V/cm.
  return 0.01 * value_vpm;
}

void PHGarfield::InitializeGas(const std::string& name)
{
  //  Create and fill the gas object so that we can trace particles through the gas...
  m_gas = new Garfield::MediumMagboltz();

  if (!std::filesystem::exists(name))
  {
    std::cout << "Missing gas file or gas directory: " << name << std::endl;
    return;
  }

  if (fs::is_regular_file(name))
  {
    std::cout << "Loading Garfield gas from file: " << name << std::endl;
    if (!m_gas->LoadGasFile(name))
    {
      std::cout << "Failed to load " << name << std::endl;
      return;
    }
    m_GasFilesLoaded = true;
  }
  else if (fs::is_directory(name))
  {
    std::cout << "Loading Garfield gas from directory: " << name << std::endl;
    std::regex filePattern(R"(^MERGED_E([0-9]{3})\.gas$)");
    std::smatch matchResults;

    // Iterate through all items in the directory
    // NOTE:  Map assures that files are properly ordered when merged...
    std::map<unsigned int, std::string> FilesToMerge;
    for (const auto& entry : fs::directory_iterator(name))
    {
      // Only process regular files
      if (entry.is_regular_file())
      {
        std::string filepath = entry.path().string();
        std::string filename = entry.path().filename().string();

        // Check if the filename matches our target pattern
        if (std::regex_match(filename, matchResults, filePattern))
        {
          // std::cout << "matchResults: " << matchResults[1].str() << std::endl;
          FilesToMerge[std::stoul(matchResults[1].str())] = filepath;
        }
      }
    }
    bool firstE = true;
    for (const auto& [key, filepath] : FilesToMerge)
    {
      if (firstE)
      {
        m_gas->LoadGasFile(filepath);
        firstE = false;
        m_GasFilesLoaded = true;
      }
      else
      {
        m_gas->MergeGasFile(filepath, true);
        m_GasFilesLoaded = true;
      }
    }
  }

  PrintGasSummary();
}

int PHGarfield::process_event(PHCompositeNode* topNode)
{
  // Initial implementation doesn't do anything event-by-event.
  // Nonetheless, a future user might want do do something here...
  (void) topNode;
  return Fun4AllReturnCodes::EVENT_OK;
}

double PHGarfield::bounder(double phi, double phi_min)
{
  double phi_max = phi_min + 2.0 * std::numbers::pi;
  while (phi < phi_min)
  {
    phi = phi + 2.0 * std::numbers::pi;
  }
  while (phi >= phi_max)
  {
    phi = phi - 2.0 * std::numbers::pi;
  }

  return phi;
}

TPolyLine3D* PHGarfield::ReverseDrift(double x, double y, double z, double step_ns, ReverseDriftStatus* driftStatus)
{
  // Public/default drift function: input and output are local TPC coordinates.
  // Do not transform these coordinates before asking Garfield for E/gas values;
  // otherwise the pre-generated gas tables and nominal TPC E field can be sampled
  // outside their intended frame.
  std::vector<double> xlist;
  std::vector<double> ylist;
  std::vector<double> zlist;

  xlist.push_back(x);
  ylist.push_back(y);
  zlist.push_back(z);

  double ex = 0.0;
  double ey = 0.0;
  double ez = 0.0;
  double bx = 0.0;
  double by = 0.0;
  double bz = 0.0;
  double vx = 0.0;
  double vy = 0.0;
  double vz = 0.0;

  // Use NaN so the first StopHere call does not trigger the z==zPrevious guard.
  double zPrevious = std::numeric_limits<double>::quiet_NaN();

  constexpr int maxSteps = 20000;
  int step = 0;

  ReverseDriftStatus status = StopHere(x, y, z, zPrevious);
  while (status == ReverseDriftStatus::Running && step < maxSteps)
  {
    zPrevious = z;

    GetMagneticFieldTesla(x, y, z, bx, by, bz);
    GetElectricFieldVcm(x, y, z, ex, ey, ez);
    m_gas->ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);

    if (!std::isfinite(vx) || !std::isfinite(vy) || !std::isfinite(vz) ||
        (vx == 0.0 && vy == 0.0 && vz == 0.0))
    {
      std::cout << PHWHERE
                << " Garfield returned invalid/zero electron velocity; stopping drift. "
                << " p_tpc=(" << x << ", " << y << ", " << z << ") cm"
                << " E=(" << ex << ", " << ey << ", " << ez << ") V/cm"
                << " B=(" << bx << ", " << by << ", " << bz << ") T"
                << std::endl;
      status = ReverseDriftStatus::InvalidVelocity;
      break;
    }

    x = x - vx * step_ns;
    y = y - vy * step_ns;
    z = z - vz * step_ns;

    xlist.push_back(x);
    ylist.push_back(y);
    zlist.push_back(z);
    ++step;

    status = StopHere(x, y, z, zPrevious);
  }

  if (step >= maxSteps && status == ReverseDriftStatus::Running)
  {
    std::cout << PHWHERE << " ReverseDrift reached maxSteps=" << maxSteps
              << "; stopping drift at p_tpc=(" << x << ", " << y << ", " << z
              << ") cm" << std::endl;

    status = ReverseDriftStatus::MaxSteps;
  }

  if (driftStatus)
  {
    *driftStatus = status;
  }

  // Keep all stored points, including the initial point.  This avoids returning
  // an empty TPolyLine3D if the starting point is already outside the active volume.
  // Drop the final point. It is the first point that triggered StopHere
  // on the next loop iteration, for example after crossing the central membrane.
  // This restores the old PHGarfield behavior.
  const size_t nGoodPoints = (xlist.size() > 1) ? xlist.size() - 1 : xlist.size();

  TPolyLine3D* poly = new TPolyLine3D(static_cast<int>(nGoodPoints));
  for (size_t i = 0; i < nGoodPoints; ++i)
  {
    poly->SetPoint(static_cast<int>(i), xlist[i], ylist[i], zlist[i]);
  }

  return poly;
}

TPolyLine3D* PHGarfield::ReverseDriftGlobalCoords(double x_cm, double y_cm, double z_cm, double step_ns, ReverseDriftStatus* status)
{
  // Caller gives a GLOBAL point and receives a GLOBAL-coordinate polyline.
  // Actual Garfield drift is still done in local TPC coordinates.

  const TVector3 p_tpc = GlobalPointToTpcPoint(x_cm, y_cm, z_cm);

  TPolyLine3D* poly = ReverseDrift(p_tpc.X(), p_tpc.Y(), p_tpc.Z(), step_ns, status);
  if (!poly)
  {
    return nullptr;
  }

  const int nPoints = poly->GetN();
  Float_t* points = poly->GetP();

  for (int i = 0; i < nPoints; ++i)
  {
    const double px = points[3 * i + 0];
    const double py = points[3 * i + 1];
    const double pz = points[3 * i + 2];

    const TVector3 p_global = TpcPointToGlobalPoint(px, py, pz);
    poly->SetPoint(i, p_global.X(), p_global.Y(), p_global.Z());
  }

  return poly;
}

PHGarfield::ReverseDriftStatus PHGarfield::StopHere(const double x, const double y, const double z,
                                                    const double zPrevious)
{
  const double r = std::hypot(x, y);

  if (r < 18.0)
  {
    return ReverseDriftStatus::RadialBoundary;
  }
  if (r > 82.0)
  {
    return ReverseDriftStatus::RadialBoundary;
  }
  if (z > 120.0)
  {
    return ReverseDriftStatus::ZBoundary;
  }
  if (z < -120.0)
  {
    return ReverseDriftStatus::ZBoundary;
  }

  // The first step passes zPrevious=NaN.  Do not apply crossing/stuck checks then.
  if (std::isfinite(zPrevious))
  {
    // z crossed the central membrane.
    if (z * zPrevious < 0.0)
    {
      return ReverseDriftStatus::CentralMembrane;
    }

    // Protect against infinite loops if Garfield returns no longitudinal motion.
    if (z == zPrevious)
    {
      return ReverseDriftStatus::Stuck;
    }
  }

  return ReverseDriftStatus::Running;
}
