#include "../PHGarfield.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <iostream>
#include <limits>

int main()
{
  PHGarfield garfield("forward_transport_config");
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

  PHGarfield::ForwardDriftResult result;
  garfield.SetForwardTransportMode(true);
  if (!garfield.GetForwardTransportMode() ||
      garfield.ForwardDrift(45.0, 0.0, 50.0, 102.605, 12.5, result) ||
      result.status != PHGarfield::ForwardDriftStatus::NotInitialized)
  {
    std::cerr << "Enabling forward mode must not bypass initialization\n";
    return 1;
  }

  if (garfield.InitRun(nullptr) != Fun4AllReturnCodes::ABORTRUN)
  {
    std::cerr << "Standalone forward initialization must reject unspecified side coefficients\n";
    return 1;
  }

  return 0;
}
