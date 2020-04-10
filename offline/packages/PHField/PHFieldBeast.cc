#include "PHFieldBeast.h"

#include <BeastMagneticField.h>

#include <Geant4/G4SystemOfUnits.hh>

#include <TSystem.h>

#include <iostream>

using namespace std;

PHFieldBeast::PHFieldBeast(const string &filename, const int verb, const float magfield_rescale)
  : PHField(verb)
  , m_MagFieldScale(magfield_rescale)
{
  m_BeastMagneticField = new BeastMagneticField(filename.c_str());
  if (m_BeastMagneticField->ValidMapImported())
  {
    // Turn linear interpolation on;
    m_BeastMagneticField->UseInterpolation();
  }
  else
  {
    cout << "error reading " << filename << endl;
    gSystem->Exit(1);
  }
}

void PHFieldBeast::GetFieldValue(const double point[4], double *Bfield) const
{
  double x = point[0] / cm;
  double y = point[1] / cm;
  double z = point[2] / cm;

  if (!m_BeastMagneticField->GetFieldValue(x, y, z, Bfield[0], Bfield[1], Bfield[2]))
  {
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;
  }
  else
  {
    Bfield[0] *= (tesla * m_MagFieldScale);
    Bfield[1] *= (tesla * m_MagFieldScale);
    Bfield[2] *= (tesla * m_MagFieldScale);
  }
  if (Verbosity() > 2)
  {
    cout << "PHFieldBeast::GetFieldValue: bx: " << Bfield[0] / tesla
         << ", by: " << Bfield[1] / tesla << ", bz: " << Bfield[2] / tesla << endl;
  }
  return;
}
