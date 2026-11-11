#include "../PHGarfieldForwardTransport.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <iostream>
#include <limits>
#include <sstream>

int main()
{
  PHGarfieldForwardTransport garfield("forward_transport_config");
  if (garfield.GetForwardTransportMode())
  {
    std::cerr << "Forward transport must require explicit enablement\n";
    return 1;
  }

  garfield.SetZeroField(true);
  if (!garfield.GetZeroMagneticField())
  {
    std::cerr << "Upstream zero-field setter must control forward B as well\n";
    return 1;
  }
  garfield.SetZeroMagneticField(false);
  if (garfield.GetZeroMagneticField())
  {
    std::cerr << "Forward zero-field alias must use the same setting\n";
    return 1;
  }

  garfield.SetDriftVolumeRadii(40.0, 50.0);
  garfield.SetCentralMembraneHalfWidth(0.28);
  const double first_step = std::numeric_limits<double>::quiet_NaN();
  if (garfield.StopHere(30.0, 0.0, 0.1, first_step) !=
      PHGarfield::ReverseDriftStatus::Running)
  {
    std::cerr << "Forward geometry must not change upstream reverse boundaries\n";
    return 1;
  }
  if (garfield.StopHere(17.0, 0.0, 10.0, first_step) !=
      PHGarfield::ReverseDriftStatus::RadialBoundary)
  {
    std::cerr << "Upstream reverse inner-cage boundary changed\n";
    return 1;
  }

  PHGarfieldForwardTransport::ForwardDriftResult result;
  garfield.SetForwardTransportMode(true);
  if (!garfield.GetForwardTransportMode() ||
      garfield.ForwardDrift(45.0, 0.0, 50.0, 102.605, 12.5, result) ||
      result.status != PHGarfieldForwardTransport::ForwardDriftStatus::NotInitialized)
  {
    std::cerr << "Enabling forward mode must not bypass initialization\n";
    return 1;
  }

  if (garfield.InitRun(nullptr) != Fun4AllReturnCodes::ABORTRUN)
  {
    std::cerr << "Standalone forward initialization must reject unspecified side coefficients\n";
    return 1;
  }

  // No external gas/CDB payload is needed: allowed field policies reach the
  // deliberately missing explicit gas file, while an unpinned P/T policy must
  // stop before looking for conditions on a null node tree.
  for (int policy : {0, 1, 2})
  {
    PHGarfieldForwardTransport standalone("standalone_field_policy", "", 0., 0.);
    standalone.SetForwardTransportMode(true);
    standalone.SetZeroMagneticField(true);
    standalone.SetLoadTpcFeeMap(false);
    standalone.SetUseCDBDriftCorrections(false);
    standalone.SetPTReference(390., 1000., 293.15);
    standalone.SetGasFile("missing_forward_field_policy_test.gas");
    if (policy == 1)
    {
      standalone.SetCMVoltageDefault(375.);
    }
    else if (policy == 2)
    {
      standalone.SetUsePTCorrection(false);
    }
    std::ostringstream errors;
    auto* previous = std::cerr.rdbuf(errors.rdbuf());
    const int status = standalone.InitRun(nullptr);
    std::cerr.rdbuf(previous);
    const char* expected = policy == 0
                               ? "Standalone forward transport requires a pinned CM field"
                               : "Explicit Garfield gas path does not exist";
    if (status != Fun4AllReturnCodes::ABORTRUN || errors.str().find(expected) == std::string::npos ||
        (policy == 1 && standalone.GetCMVoltageDefault() != 375.) ||
        (policy == 2 && standalone.GetCMVoltageDefault() != 390.))
    {
      std::cerr << "Standalone P/T policy failed: " << policy << ' ' << errors.str();
      return 1;
    }
  }

  return 0;
}
