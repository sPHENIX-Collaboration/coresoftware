#include <cstdlib>
#include <iostream>
#include <string>

#include <Garfield/MediumMagboltz.hh>

//------------------------------------------------------------
//  This standalone executable makes gas calculations for
//  whatever mixture you specify and range of electric,
//  magnetic, and angle between values that you select.
//                           TKH   5/27/2026
//
//
//------------------------------------------------------------

int main(int argc, char* argv[])
{
  if (argc != 11)
  {
    std::cerr
        << "Usage:\n"
        << argv[0]
        << " Emin Emax nE Bmin Bmax nB Amin Amax nA output_file_name\n\n"
        << "Units:\n"
        << "  E: V/cm\n"
        << "  B: Tesla\n"
        << "  angle: radians\n";
    return 1;
  }

  const double Emin = std::atof(argv[1]);
  const double Emax = std::atof(argv[2]);
  const int nE = std::atoi(argv[3]);

  const double Bmin = std::atof(argv[4]);
  const double Bmax = std::atof(argv[5]);
  const int nB = std::atoi(argv[6]);

  const double Amin = std::atof(argv[7]);
  const double Amax = std::atof(argv[8]);
  const int nA = std::atoi(argv[9]);

  const std::string output_file(argv[10]);

  std::cout << "E grid: "
            << Emin << " -> " << Emax
            << " with " << nE << " points\n";

  std::cout << "B grid: "
            << Bmin << " -> " << Bmax
            << " with " << nB << " points\n";

  std::cout << "Angle grid: "
            << Amin << " -> " << Amax
            << " with " << nA << " points\n";

  std::cout << "Output File: " << output_file << std::endl;

  // ------------------------------------------------------------
  // Gas: Ar/CF4/isobutane = 75/20/5.
  // ------------------------------------------------------------
  Garfield::MediumMagboltz gas;
  gas.SetComposition("ar", 75., "cf4", 20., "isobutane", 5.);
  gas.SetTemperature(301.65);  // K     from Grafana
  gas.SetPressure(762.);       // Torr  from Grafana

  // Try to load an existing gas table.
  bool LogGrid = false;
  gas.SetFieldGrid(Emin, Emax, nE, LogGrid,
                   Bmin, Bmax, nB,
                   Amin, Amax, nA);

  gas.GenerateGasTable(10);
  gas.WriteGasFile(output_file);
}
