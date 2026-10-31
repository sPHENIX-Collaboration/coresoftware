#include "../PHGarfield.h"

#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <numbers>
#include <sstream>

// Exercise the real loader and both field paths without external gas or CDB.
struct PHGarfieldForwardFieldTest
{
  static bool Coverage(const char* filename, bool three_dimensional)
  {
    PHGarfield field("field_coverage", filename, 1., 1.);
    if (three_dimensional
            ? (!field.LoadElectricFieldCorrections3D(filename, 0) || !field.LoadElectricFieldCorrections3D(filename, 1))
            : !field.LoadElectricFieldCorrections(filename))
    {
      return false;
    }
    std::ostringstream error;
    const auto accepts = [&](double inner, double outer, double distance)
    {
      error.str("");
      return field.ForwardElectricFieldMapsCover(inner, outer, distance, error);
    };
    // The map ends at gas-surface z=30, matching readout 30.28 and CM hw=.28.
    if (!accepts(40., 60., 30.28 - .28) || !error.str().empty() ||
        accepts(40., 60., 30.28) || error.str().empty() ||
        accepts(39.9, 60., 30.) || accepts(40., 60.1, 30.))
    {
      return false;
    }
    if (three_dimensional)
    {
      // Check every component on both sides, including a CM-center-origin
      // payload whose low edge is .28 rather than zero in gas coordinates.
      for (std::size_t side = 0; side < 2; ++side)
      {
        for (auto* component : field.m_field3DCorrection[side])
        {
          component->GetZaxis()->SetLimits(.28, 30.28);
          if (accepts(40., 60., 30.) || error.str().find("side " + std::to_string(side)) == std::string::npos)
          {
            return false;
          }
          component->GetZaxis()->SetLimits(0., 29.9);
          if (accepts(40., 60., 30.))
          {
            return false;
          }
          component->GetZaxis()->SetLimits(0., 30.);
        }
      }
    }
    else
    {
      field.m_erCorrection->GetYaxis()->SetLimits(.28, 30.28);
      if (accepts(40., 60., 30.))
      {
        return false;
      }
    }
    return true;
  }

  static bool Run(const char* filename, bool three_dimensional)
  {
    PHGarfield field("field_coordinates", filename, 1., 1.);
    field.SetCMVoltageDefault(375.);
    field.SetUseIFCVoltageDistortion(false);
    field.SetUseOFCVoltageDistortion(false);
    if (three_dimensional
            ? (!field.LoadElectricFieldCorrections3D(filename, 0) || !field.LoadElectricFieldCorrections3D(filename, 1))
            : !field.LoadElectricFieldCorrections(filename))
    {
      return false;
    }
    for (double half_width : {0., .28, 1.})
    {
      field.SetCentralMembraneHalfWidth(half_width);
      for (double height : {2., 10., 20.})
        for (double sign : {-1., 1.})
        {
          double fx = 0., fy = 0., fz = 0., rx = 0., ry = 0., rz = 0.;
          field.GetForwardElectricFieldVcm(45., 0., sign * height, fx, fy, fz);
          field.GetElectricFieldVcm(45., 0., sign * height, rx, ry, rz);
          const double forward_z = height - half_width;
          const auto close = [](double a, double b)
          { return std::abs(a - b) < 1.e-9; };
          if (!close(fx, 3. + 2. * forward_z) || !close(fy, three_dimensional ? 7. + 4. * forward_z : 0.) ||
              !close(fz, sign * (-375. + 5. + 3. * forward_z)) ||
              !close(rx, 3. + 2. * height) || !close(ry, three_dimensional ? 7. + 4. * height : 0.) ||
              !close(rz, sign * (-375. + 5. + 3. * height)))
          {
            std::cerr << "Forward/reverse space-charge coordinate convention changed\n";
            return false;
          }
        }
    }
    return true;
  }
};

int main()
{
  char filename[] = "forward-field-XXXXXX";
  const int descriptor = mkstemp(filename);
  if (descriptor < 0)
  {
    return 1;
  }
  close(descriptor);
  {
    TFile output(filename, "RECREATE");
    TH2D er("hErDefault", "Known linear field in V/m", 4, 40., 60., 60, 0., 30.);
    TH2D ez("hEzDefault", "Known linear field in V/m", 4, 40., 60., 60, 0., 30.);
    for (int r = 1; r <= er.GetNbinsX(); ++r)
      for (int z = 1; z <= er.GetNbinsY(); ++z)
      {
        const double coordinate = er.GetYaxis()->GetBinCenter(z);
        er.SetBinContent(r, z, 100. * (3. + 2. * coordinate));
        ez.SetBinContent(r, z, 100. * (5. + 3. * coordinate));
      }
    er.Write();
    ez.Write();
    output.mkdir("Field3D")->cd();
    TH3D ex3("hEx", "Known 3D field in V/m", 4, 40., 60., 4, 0., 2. * std::numbers::pi, 60, 0., 30.);
    TH3D ey3("hEy", "Known 3D field in V/m", 4, 40., 60., 4, 0., 2. * std::numbers::pi, 60, 0., 30.);
    TH3D ez3("hEz", "Known 3D field in V/m", 4, 40., 60., 4, 0., 2. * std::numbers::pi, 60, 0., 30.);
    for (int r = 1; r <= ex3.GetNbinsX(); ++r)
      for (int phi = 1; phi <= ex3.GetNbinsY(); ++phi)
        for (int z = 1; z <= ex3.GetNbinsZ(); ++z)
        {
          const double coordinate = ex3.GetZaxis()->GetBinCenter(z);
          ex3.SetBinContent(r, phi, z, 100. * (3. + 2. * coordinate));
          ey3.SetBinContent(r, phi, z, 100. * (7. + 4. * coordinate));
          ez3.SetBinContent(r, phi, z, 100. * (5. + 3. * coordinate));
        }
    ex3.Write();
    ey3.Write();
    ez3.Write();
  }
  const bool valid = PHGarfieldForwardFieldTest::Run(filename, false) && PHGarfieldForwardFieldTest::Run(filename, true) &&
                     PHGarfieldForwardFieldTest::Coverage(filename, false) && PHGarfieldForwardFieldTest::Coverage(filename, true);
  std::remove(filename);
  return valid ? 0 : 1;
}
