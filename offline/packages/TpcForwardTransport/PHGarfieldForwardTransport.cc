#include "PHGarfieldForwardTransport.h"
#include "PHGarfieldFieldMapValidation.h"
#include "PHGarfieldForwardSetup.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phfield/PHField.h>
#include <phfield/PHField3DCartesian.h>
#include <phfield/PHFieldUtility.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TAxis.h>
#include <TH2.h>
#include <TH3.h>
#include <TRotation.h>
#include <TVector3.h>

#include <Garfield/ComponentUser.hh>
#include <Garfield/MediumMagboltz.hh>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

class PHGarfieldForwardNode final : public PHDataNode<PHGarfield>
{
 public:
  explicit PHGarfieldForwardNode(PHGarfield* garfield)
    : PHDataNode<PHGarfield>(garfield, "PHGarfield", "Data")
  {
  }

  PHGarfieldForwardNode(const PHGarfieldForwardNode&) = delete;
  PHGarfieldForwardNode& operator=(const PHGarfieldForwardNode&) = delete;
  PHGarfieldForwardNode(PHGarfieldForwardNode&&) = delete;
  PHGarfieldForwardNode& operator=(PHGarfieldForwardNode&&) = delete;

  ~PHGarfieldForwardNode() override
  {
    setData(nullptr);  // Fun4AllServer deletes PHGarfield.
  }
};

PHGarfieldForwardTransport::~PHGarfieldForwardTransport()
{
  // The base owns its gas/maps, but normally borrows the magnetic field.
  if (m_ownsField)
  {
    delete m_field;
    m_field = nullptr;
  }
}

int PHGarfieldForwardTransport::InitRun(PHCompositeNode* topNode)
{
  if (!m_forwardTransportMode)
  {
    return PHGarfield::InitRun(topNode);
  }
  return InitializeForwardTransport(topNode);
}

int PHGarfieldForwardTransport::InitializeForwardTransport(PHCompositeNode* topNode)
{
  PHGarfieldForwardSetup::Request request;
  request.has_node_tree = topNode != nullptr;
  request.space_charge_scale_override = m_spaceChargeScaleOverride;
  request.zero_magnetic_field = m_zerofield;
  request.magnetic_field_override = m_magneticFieldMapOverride;
  request.load_tpc_fee_map = m_loadTpcFeeMap;
  request.electric_field_map_override = m_electricFieldMapOverride;
  request.use_2d_electric_field_map = m_use2DElectricFieldMap;
  request.electric_field_maps_3d = m_electricFieldMap3D;
  request.electric_field_map_3d_override = m_electricFieldMap3DOverride;
  request.frame_electric_field_map_override = m_frameElectricFieldMapOverride;
  request.frame_electric_field_maps_3d = m_frameElectricFieldMap3D;
  request.frame_electric_field_map_3d_override = m_frameElectricFieldMap3DOverride;
  const auto inputs = PHGarfieldForwardSetup::SelectInputs(request);
  if (!inputs.valid)
  {
    std::cerr << PHWHERE << " Standalone forward transport requires explicit side coefficients"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (!topNode && m_usePTCorrection && !m_CMVoltageOverride)
  {
    std::cerr << PHWHERE << " Standalone forward transport requires a pinned CM field"
              << " (SetCMVoltageDefault) or SetUsePTCorrection(false)" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Rebind all run-dependent payloads. The magnetic field borrowed from a node
  // tree is not owned here; the standalone fallback is.
  if (m_ownsField)
  {
    delete m_field;
  }
  m_field = nullptr;
  m_ownsField = false;
  m_hasMagneticFieldDomain = false;
  delete m_cdbTPCMAPttree;
  m_cdbTPCMAPttree = nullptr;
  delete m_component;
  m_component = nullptr;
  delete m_gas;
  m_gas = nullptr;
  m_GasFilesLoaded = false;
  delete m_erCorrection;
  m_erCorrection = nullptr;
  delete m_ezCorrection;
  m_ezCorrection = nullptr;
  ClearElectricFieldCorrections3D(0);
  ClearElectricFieldCorrections3D(1);
  delete m_frameErCorrection;
  m_frameErCorrection = nullptr;
  delete m_frameEzCorrection;
  m_frameEzCorrection = nullptr;
  ClearFrameElectricFieldCorrections3D(0);
  ClearFrameElectricFieldCorrections3D(1);

  m_electricFieldMap3D = inputs.electric_field_maps_3d;
  m_frameElectricFieldMap3D = inputs.frame_electric_field_maps_3d;

  if (Verbosity() > 1)
  {
    std::cout << "PHGarfieldForwardTransport::InitRun(PHCompositeNode *topNode) Initializing" << std::endl;
  }
  CDBInterface* m_cdb = CDBInterface::instance();
  const auto lookup = [m_cdb](const std::string& payload)
  { return m_cdb->getUrl(payload); };

  // In a Fun4All job, use the existing field node so that the map type,
  // rescaling, and uniform/custom field configuration exactly match the
  // simulation. Standalone PHGarfield tests retain the historical CDB fallback.
  switch (inputs.magnetic_field_source)
  {
  case PHGarfieldForwardSetup::MagneticFieldSource::Zero:
    m_resolvedMagneticFieldMap = PHGarfieldForwardSetup::ResolveMagneticField(inputs, lookup);
    break;
  case PHGarfieldForwardSetup::MagneticFieldSource::ExplicitFile:
  case PHGarfieldForwardSetup::MagneticFieldSource::CdbFile:
    m_resolvedMagneticFieldMap = PHGarfieldForwardSetup::ResolveMagneticField(inputs, lookup);
    if (!LoadMagneticFieldDomain(m_resolvedMagneticFieldMap))
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }
    m_field = new PHField3DCartesian(
        m_resolvedMagneticFieldMap,
        m_magneticFieldMapScale);
    m_ownsField = true;
    break;
  case PHGarfieldForwardSetup::MagneticFieldSource::NodeTree:
    m_field = PHFieldUtility::GetFieldMapNode(nullptr, topNode, Verbosity());
    m_ownsField = false;
    m_resolvedMagneticFieldMap = PHGarfieldForwardSetup::ResolveMagneticField(inputs, lookup);
    break;
  }

  // The pad map is used only by the legacy plotting helpers. Standalone
  // forward-map generation can omit this unrelated CDB dependency.
  if (inputs.load_tpc_fee_map)
  {
    const std::string text = PHGarfieldForwardSetup::ResolveFeeMap(inputs, lookup);
    m_cdbTPCMAPttree = new CDBTTree(text);
    m_cdbTPCMAPttree->LoadCalibrations();
  }

  // The standalone forward-map generator pins its field model and geometry
  // without an event node tree. Gas is resolved independently below; avoid
  // adding unrelated reconstruction CDB inputs to this configuration.
  if (inputs.use_run_dependent_inputs && !LoadCDBInputs(topNode))
  {
    std::cout << PHWHERE << " Failed to load PHGarfield CDB inputs"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (inputs.use_run_dependent_inputs)
  {
    ConfigureRunDependentFieldCage();
  }
  ConfigurePTCorrectedField(topNode);
  ConfigureRunDriftCorrections();

  // Load the optional axisymmetric space-charge correction map. A nonempty
  // filename is an explicit request; never silently replace a requested map
  // with the nominal electric field when the payload cannot be loaded.
  if (!m_electricFieldMap.empty())
  {
    if (!LoadElectricFieldCorrections(m_electricFieldMap))
    {
      std::cout << PHWHERE << " Failed to load electric-field correction map: "
                << m_electricFieldMap << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
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
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }

  // Load optional frame-charge maps independently. They are additive and do
  // not alter the existing space-charge-map selection or k_eff flow.
  if (!m_frameElectricFieldMap.empty() &&
      !LoadFrameElectricFieldCorrections(m_frameElectricFieldMap))
  {
    std::cout << PHWHERE << " Failed to load frame electric-field correction map: "
              << m_frameElectricFieldMap << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  for (std::size_t side = 0; side < m_frameElectricFieldMap3D.size(); ++side)
  {
    if (!m_frameElectricFieldMap3D[side].empty() &&
        !LoadFrameElectricFieldCorrections3D(m_frameElectricFieldMap3D[side], side))
    {
      std::cout << PHWHERE << " Failed to load side " << side
                << " 3D frame electric-field correction map: "
                << m_frameElectricFieldMap3D[side] << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
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
    std::cerr << PHWHERE << " Failed to build requested IFC/OFC boundary-voltage field grids"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  //  Make the Garfield Component and register the methods that will interface to our fields...
  m_component = new Garfield::ComponentUser();
  m_component->SetMagneticField([this](double x, double y, double z, double& bx, double& by, double& bz)
                                { GetMagneticFieldTesla(x, y, z, bx, by, bz); });
  m_component->SetElectricField([this](double x, double y, double z, double& ex, double& ey, double& ez)
                                { GetForwardElectricFieldVcm(x, y, z, ex, ey, ez); });

  // Resolve gas after the run-dependent inputs and electric-field resources,
  // preserving the original initialization and CDB lookup order.
  const auto gas = PHGarfieldForwardSetup::ResolveGas(
      {m_gasFileOverride, m_defaultGasfile}, lookup,
      [](const std::string& filename)
      { return fs::exists(filename); });
  if (!gas.valid)
  {
    std::cerr << PHWHERE << " Explicit Garfield gas path does not exist: "
              << gas.filename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (gas.source == PHGarfieldForwardSetup::GasSource::DefaultFile)
  {
    std::cout << PHWHERE << " Missing CDB gasfile: " << gas.unavailable_cdb_filename << std::endl;
    std::cout << PHWHERE << " Using default gasfile: " << m_defaultGasfile << std::endl;
  }
  m_resolvedGasFile = gas.filename;
  InitializeGas(gas.filename, true);
  if (!IsReady())
  {
    std::cerr << PHWHERE << " Failed to initialize the magnetic field or Garfield gas"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  //  Diagnostic during code development...
  if (m_cdbTPCMAPttree)
  {
    FillRadii();
    if (Verbosity() > 1)
    {
      PrintMaps();
    }
  }

  // Standalone forward transport has no node tree to publish into.
  if (!topNode)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  auto* existing = findNode::getClass<PHGarfield>(topNode, "PHGarfield");
  if (existing && existing != this)
  {
    // A shared instance is already registered; this one stays private
    // (e.g. Tpc_PolyClusterizer with setUsePHGarfieldDefaults(false)).
    std::cout << Name() << "::InitRun - another PHGarfield is already in the node tree; "
              << "this instance stays private" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (!existing && m_publishToNodeTree)
  {
    PHNodeIterator iter(topNode);
    auto* parNode = dynamic_cast<PHCompositeNode*>(
        iter.findFirst("PHCompositeNode", "PAR"));

    if (!parNode)
    {
      parNode = new PHCompositeNode("PAR");
      topNode->addNode(parNode);
    }

    parNode->addNode(new PHGarfieldForwardNode(this));
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHGarfieldForwardTransport::SetTpcRotation(double theta_x, double theta_y, double theta_z)
{
  m_tpcrot = TRotation();
  RotateTpc(theta_x, theta_y, theta_z);
}

void PHGarfieldForwardTransport::GetForwardElectricFieldVcm(double x_cm, double y_cm, double z_cm, double& ex_vcm, double& ey_vcm, double& ez_vcm, const double correction_scale) const
{
  // NOTE: Garfield uses cm, V/cm, and Tesla.
  // The notebook maps use cm on their axes and V/m in their bins.
  // Forward space-charge maps use distance from the CM gas surface,
  // max(0, |z|-CM_half_width), for both the axisymmetric and 3D models.
  // Frame/cage corrections and upstream reverse reconstruction retain their
  // existing physical-coordinate conventions.

  const double r_cm = std::hypot(x_cm, y_cm);
  const double drift_z_cm =
      std::max(0.0, std::abs(z_cm) - m_centralMembraneHalfWidth);
  const double phi_rad = std::atan2(y_cm, x_cm);
  const std::size_t side = z_cm < 0.0 ? 0 : 1;

  // Nominal drift field. Electrons drift away from the central membrane,
  // therefore E_z points toward the central membrane on both sides.
  ex_vcm = 0.0;
  ey_vcm = 0.0;
  const double nominal_vcm = m_CMVoltageDefault * m_cmFieldFactor[side];
  ez_vcm = z_cm > 0.0 ? -nominal_vcm : nominal_vcm;

  // A forward scale-zero reference must disable all optional corrections,
  // including the upstream frame and field-cage terms.
  if (correction_scale == 0.0)
  {
    return;
  }

  const double spaceChargeScale = correction_scale *
                                  (side == 0 ? m_spaceChargeScale_side0 : m_spaceChargeScale_side1);
  if (spaceChargeScale == 0.0)
  {
    AddFrameElectricFieldCorrections(r_cm, phi_rad, z_cm, ex_vcm, ey_vcm, ez_vcm);
    AddFieldCageVoltageDistortion(x_cm, y_cm, z_cm, ex_vcm, ey_vcm, ez_vcm);
    return;
  }

  // Prefer the side-specific 3D map when all three Cartesian components are
  // available. Ex and Ey are already local-TPC Cartesian components. hEz is
  // stored along the longitudinal solver coordinate and reflected on side 0.
  if (HasElectricFieldCorrections3D(side))
  {
    ex_vcm += spaceChargeScale *
              InterpolateCorrectionVcm(m_field3DCorrection[side][0], r_cm, phi_rad, drift_z_cm);
    ey_vcm += spaceChargeScale *
              InterpolateCorrectionVcm(m_field3DCorrection[side][1], r_cm, phi_rad, drift_z_cm);

    const double delta_ez_local_vcm = spaceChargeScale *
                                      InterpolateCorrectionVcm(m_field3DCorrection[side][2], r_cm, phi_rad, drift_z_cm);
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
                              InterpolateCorrectionVcm(m_erCorrection, r_cm, drift_z_cm);
  const double delta_ez_local_vcm = spaceChargeScale *
                                    InterpolateCorrectionVcm(m_ezCorrection, r_cm, drift_z_cm);

  // Convert the cylindrical radial correction to Cartesian components.
  if (r_cm > 0.0)
  {
    ex_vcm += delta_er_vcm * x_cm / r_cm;
    ey_vcm += delta_er_vcm * y_cm / r_cm;
  }

  // Reflect the solver's longitudinal component onto the selected TPC side.
  // Convert it to the local TPC Cartesian z direction.
  ez_vcm += z_cm >= 0.0 ? delta_ez_local_vcm : -delta_ez_local_vcm;

  AddFrameElectricFieldCorrections(r_cm, phi_rad, z_cm, ex_vcm, ey_vcm, ez_vcm);
  AddFieldCageVoltageDistortion(x_cm, y_cm, z_cm, ex_vcm, ey_vcm, ez_vcm);
}

bool PHGarfieldForwardTransport::ElectricFieldMapCovers(
    const double inner_radius_cm,
    const double outer_radius_cm,
    const double max_drift_distance_cm) const
{
  if (!m_erCorrection || !m_ezCorrection ||
      !std::isfinite(inner_radius_cm) ||
      !std::isfinite(outer_radius_cm) ||
      !std::isfinite(max_drift_distance_cm))
  {
    return false;
  }

  constexpr double tolerance_cm = 1.0e-6;
  const auto covers = [=](const TH2* hist)
  {
    const TAxis* r_axis = hist->GetXaxis();
    const TAxis* z_axis = hist->GetYaxis();
    return r_axis && z_axis &&
           r_axis->GetXmin() <= inner_radius_cm + tolerance_cm &&
           r_axis->GetXmax() >= outer_radius_cm - tolerance_cm &&
           z_axis->GetXmin() <= tolerance_cm &&
           z_axis->GetXmax() >=
               max_drift_distance_cm - tolerance_cm;
  };

  return covers(m_erCorrection) && covers(m_ezCorrection);
}

bool PHGarfieldForwardTransport::ForwardElectricFieldMapsCover(
    const double inner_radius_cm,
    const double outer_radius_cm,
    const double max_drift_distance_cm,
    std::ostream& error) const
{
  if (!std::isfinite(inner_radius_cm) || !std::isfinite(outer_radius_cm) ||
      !std::isfinite(max_drift_distance_cm) || inner_radius_cm <= 0. ||
      outer_radius_cm <= inner_radius_cm || max_drift_distance_cm <= 0.)
  {
    error << "Invalid forward electric-field coverage geometry\n";
    return false;
  }
  constexpr double tolerance_cm = 1.e-6;
  const auto covers = [=](const TAxis* radius, const TAxis* longitudinal)
  {
    return radius && longitudinal &&
           radius->GetXmin() <= inner_radius_cm + tolerance_cm &&
           radius->GetXmax() >= outer_radius_cm - tolerance_cm &&
           longitudinal->GetXmin() <= tolerance_cm &&
           longitudinal->GetXmax() >= max_drift_distance_cm - tolerance_cm;
  };
  for (std::size_t side = 0; side < m_field3DCorrection.size(); ++side)
  {
    if (HasElectricFieldCorrections3D(side))
    {
      for (const auto* component : m_field3DCorrection[side])
      {
        if (!covers(component->GetXaxis(), component->GetZaxis()))
        {
          error << "Forward 3D space-charge map for side " << side;
          error << " does not cover r=[" << inner_radius_cm << ", " << outer_radius_cm
                << "] cm and gas-surface z=[0, " << max_drift_distance_cm << "] cm\n";
          return false;
        }
      }
    }
    else if (HasElectricFieldCorrections() &&
             !ElectricFieldMapCovers(inner_radius_cm, outer_radius_cm, max_drift_distance_cm))
    {
      error << "Forward 2D space-charge map does not cover r=["
            << inner_radius_cm << ", " << outer_radius_cm << "] cm and gas-surface z=[0, "
            << max_drift_distance_cm << "] cm\n";
      return false;
    }
  }
  return true;
}

bool PHGarfieldForwardTransport::LoadMagneticFieldDomain(const std::string& filename)
{
  m_hasMagneticFieldDomain = false;
  PHGarfieldFieldMapValidation::CartesianDomain domain;
  if (!PHGarfieldFieldMapValidation::InspectCartesianMap(filename, domain, std::cerr))
  {
    return false;
  }
  m_magneticFieldMinimumCm = domain.minimum_cm;
  m_magneticFieldMaximumCm = domain.maximum_cm;
  m_hasMagneticFieldDomain = true;
  return true;
}

bool PHGarfieldForwardTransport::MagneticFieldDomainContains(
    const double x_cm,
    const double y_cm,
    const double z_cm) const
{
  if (m_zerofield || !m_hasMagneticFieldDomain)
  {
    return true;
  }

  const TVector3 point =
      TpcPointToMagnetFieldMapPoint(x_cm, y_cm, z_cm);
  const std::array<double, 3> coordinate = {
      {point.X(), point.Y(), point.Z()}};
  for (unsigned int axis = 0; axis < coordinate.size(); ++axis)
  {
    if (coordinate[axis] < m_magneticFieldMinimumCm[axis] ||
        coordinate[axis] > m_magneticFieldMaximumCm[axis])
    {
      return false;
    }
  }
  return true;
}

bool PHGarfieldForwardTransport::MagneticFieldMapCoversTpcDriftVolume(
    const double outer_radius_cm,
    const double readout_abs_z_cm) const
{
  if (!std::isfinite(outer_radius_cm) ||
      !std::isfinite(readout_abs_z_cm) ||
      outer_radius_cm <= 0.0 ||
      readout_abs_z_cm <= 0.0 ||
      !HasKnownMagneticFieldDomain())
  {
    return false;
  }
  if (m_zerofield)
  {
    return true;
  }

  // The local box encloses the full cylindrical gas volume. Requiring every
  // corner to lie inside the Cartesian field grid is conservative and, since
  // both boxes and the rigid transform are convex, guarantees full coverage.
  for (const double x_cm : {-outer_radius_cm, outer_radius_cm})
  {
    for (const double y_cm : {-outer_radius_cm, outer_radius_cm})
    {
      for (const double z_cm :
           {-readout_abs_z_cm, readout_abs_z_cm})
      {
        if (!MagneticFieldDomainContains(x_cm, y_cm, z_cm))
        {
          return false;
        }
      }
    }
  }
  return true;
}

bool PHGarfieldForwardTransport::ForwardDrift(double x_cm,
                                              double y_cm,
                                              double z_cm,
                                              const double readout_abs_z_cm,
                                              const double step_ns,
                                              ForwardDriftResult& result,
                                              const bool include_electric_field_corrections) const
{
  result = ForwardDriftResult{};
  result.x_cm = x_cm;
  result.y_cm = y_cm;
  result.z_cm = z_cm;

  if (!m_forwardTransportMode || !IsReady())
  {
    result.status = ForwardDriftStatus::NotInitialized;
    return false;
  }

  if (!std::isfinite(x_cm) ||
      !std::isfinite(y_cm) ||
      !std::isfinite(z_cm) ||
      !std::isfinite(readout_abs_z_cm) ||
      !std::isfinite(step_ns) ||
      readout_abs_z_cm <= m_centralMembraneHalfWidth ||
      step_ns <= 0.0)
  {
    result.status = ForwardDriftStatus::InvalidInput;
    return false;
  }

  // Keep the side fixed for the entire drift. The active gas does not include
  // the central membrane, so z == 0 is not a valid forward-drift start.
  constexpr double z_tolerance_cm = 1.0e-9;
  if (z_cm == 0.0 ||
      std::abs(z_cm) <
          m_centralMembraneHalfWidth - z_tolerance_cm)
  {
    result.status = ForwardDriftStatus::OutsideDriftVolume;
    return false;
  }

  const double side_sign = z_cm > 0.0 ? 1.0 : -1.0;
  const double target_z_cm = side_sign * readout_abs_z_cm;
  constexpr double radial_tolerance_cm = 1.0e-9;

  // Reject points beyond the readout plane and points outside the drift cage.
  const double start_radius_cm = std::hypot(x_cm, y_cm);
  if (side_sign * (z_cm - target_z_cm) > z_tolerance_cm ||
      start_radius_cm < m_driftInnerRadius - radial_tolerance_cm ||
      start_radius_cm > m_driftOuterRadius + radial_tolerance_cm)
  {
    result.status = ForwardDriftStatus::OutsideDriftVolume;
    return false;
  }
  if (!MagneticFieldDomainContains(x_cm, y_cm, z_cm))
  {
    result.status = ForwardDriftStatus::LeftMagneticFieldMap;
    return false;
  }

  if (std::abs(target_z_cm - z_cm) <= z_tolerance_cm)
  {
    result.z_cm = target_z_cm;
    result.status = ForwardDriftStatus::ReachedReadout;
    result.reached_readout = true;
    return true;
  }

  double drift_time_ns = 0.0;

  for (unsigned int step = 0; step < MaxForwardDriftSteps; ++step)
  {
    double ex = 0.0;
    double ey = 0.0;
    double ez = 0.0;
    double bx = 0.0;
    double by = 0.0;
    double bz = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;

    GetMagneticFieldTesla(x_cm, y_cm, z_cm, bx, by, bz);
    GetForwardElectricFieldVcm(
        x_cm,
        y_cm,
        z_cm,
        ex,
        ey,
        ez,
        include_electric_field_corrections ? 1.0 : 0.0);

    const bool velocity_ok =
        m_gas->ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
    if (!velocity_ok ||
        !std::isfinite(vx) ||
        !std::isfinite(vy) ||
        !std::isfinite(vz) ||
        (vx == 0.0 && vy == 0.0 && vz == 0.0))
    {
      result.status = ForwardDriftStatus::InvalidVelocity;
      return false;
    }

    const double remaining_z_cm = target_z_cm - z_cm;
    const double time_to_readout_ns = remaining_z_cm / vz;
    if (!std::isfinite(time_to_readout_ns) || time_to_readout_ns <= 0.0)
    {
      result.status = ForwardDriftStatus::WrongDirection;
      return false;
    }

    const bool reaches_readout = time_to_readout_ns <= step_ns;
    const double this_step_ns =
        reaches_readout ? time_to_readout_ns : step_ns;

    const double previous_x_cm = x_cm;
    const double previous_y_cm = y_cm;
    x_cm += vx * this_step_ns;
    y_cm += vy * this_step_ns;
    z_cm += vz * this_step_ns;
    drift_time_ns += this_step_ns;

    if (reaches_readout)
    {
      // Avoid leaving a small floating-point offset from the physical plane.
      z_cm = target_z_cm;
    }

    result.x_cm = x_cm;
    result.y_cm = y_cm;
    result.z_cm = z_cm;
    result.drift_time_ns = drift_time_ns;
    result.steps = step + 1;

    const double radius_cm = std::hypot(x_cm, y_cm);
    if (!std::isfinite(radius_cm))
    {
      result.status = ForwardDriftStatus::InvalidVelocity;
      return false;
    }
    if (!MagneticFieldDomainContains(x_cm, y_cm, z_cm))
    {
      result.status = ForwardDriftStatus::LeftMagneticFieldMap;
      return false;
    }
    // Forward transport always stops at a radial field cage.
    {
      // The Euler step is a straight segment. Test the whole segment rather
      // than only its endpoint: a chord can cross the inner cage and return
      // to an allowed endpoint.
      const double delta_x_cm = x_cm - previous_x_cm;
      const double delta_y_cm = y_cm - previous_y_cm;
      const double segment_length_squared =
          delta_x_cm * delta_x_cm + delta_y_cm * delta_y_cm;
      double closest_fraction = 0.0;
      if (segment_length_squared > 0.0)
      {
        closest_fraction = std::clamp(
            -(previous_x_cm * delta_x_cm +
              previous_y_cm * delta_y_cm) /
                segment_length_squared,
            0.0,
            1.0);
      }
      const double closest_x_cm =
          previous_x_cm + closest_fraction * delta_x_cm;
      const double closest_y_cm =
          previous_y_cm + closest_fraction * delta_y_cm;
      const double minimum_segment_radius_cm =
          std::hypot(closest_x_cm, closest_y_cm);
      const double previous_radius_cm =
          std::hypot(previous_x_cm, previous_y_cm);
      const double maximum_segment_radius_cm =
          std::max(previous_radius_cm, radius_cm);

      if (!std::isfinite(minimum_segment_radius_cm) ||
          !std::isfinite(maximum_segment_radius_cm))
      {
        result.status = ForwardDriftStatus::InvalidVelocity;
        return false;
      }
      if (minimum_segment_radius_cm <
              m_driftInnerRadius - radial_tolerance_cm ||
          maximum_segment_radius_cm >
              m_driftOuterRadius + radial_tolerance_cm)
      {
        result.status = ForwardDriftStatus::LeftDriftVolume;
        return false;
      }
    }

    if (reaches_readout)
    {
      result.status = ForwardDriftStatus::ReachedReadout;
      result.reached_readout = true;
      return true;
    }
  }

  if (Verbosity() > 0)
  {
    std::cerr << PHWHERE << " ForwardDrift reached max_steps="
              << MaxForwardDriftSteps
              << " at p_tpc=(" << x_cm << ", " << y_cm << ", " << z_cm
              << ") cm" << std::endl;
  }

  result.status = ForwardDriftStatus::MaxSteps;
  return false;
}

const char* PHGarfieldForwardTransport::ForwardDriftStatusName(
    const ForwardDriftStatus status)
{
  switch (status)
  {
  case ForwardDriftStatus::NotRun:
    return "not run";
  case ForwardDriftStatus::ReachedReadout:
    return "reached readout";
  case ForwardDriftStatus::NotInitialized:
    return "not initialized";
  case ForwardDriftStatus::InvalidInput:
    return "invalid input";
  case ForwardDriftStatus::OutsideDriftVolume:
    return "outside drift volume";
  case ForwardDriftStatus::InvalidVelocity:
    return "invalid Garfield velocity";
  case ForwardDriftStatus::WrongDirection:
    return "wrong longitudinal direction";
  case ForwardDriftStatus::LeftDriftVolume:
    return "left drift volume";
  case ForwardDriftStatus::LeftMagneticFieldMap:
    return "left magnetic-field map domain";
  case ForwardDriftStatus::MaxSteps:
    return "maximum steps reached";
  }
  return "unknown";
}
