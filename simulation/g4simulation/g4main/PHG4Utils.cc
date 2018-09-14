#include "PHG4Utils.h"

#include <Geant4/G4VisAttributes.hh>

#include <cmath>

using namespace std;

double PHG4Utils::_eta_coverage = 1.;

double
PHG4Utils::GetLengthForRapidityCoverage(const double radius, const double eta)
{
  double length;
  double theta = 2.0 * std::atan(std::exp(-eta));
  length = radius / std::tan(theta);
  return length;
}

double
PHG4Utils::GetLengthForRapidityCoverage(const double radius)
{
  return GetLengthForRapidityCoverage(radius, _eta_coverage);
}

void PHG4Utils::SetPseudoRapidityCoverage(const double eta)
{
  _eta_coverage = eta;
}

double
PHG4Utils::get_theta(const double eta)
{
  double theta = 2 * atan(exp(-eta));
  return theta;
}

double
PHG4Utils::get_eta(const double theta)
{
  double eta = -log(tan(theta / 2.));
  return eta;
}

pair<double, double>
PHG4Utils::get_etaphi(const double x, const double y, const double z)
{
  double eta;
  double phi;
  double radius;
  double theta;
  radius = sqrt(x * x + y * y);
  phi = atan2(y, x);
  theta = atan2(radius, z);
  eta = -log(tan(theta / 2.));
  return make_pair(eta, phi);
}

double
PHG4Utils::get_eta(const double radius, const double z)
{
  double eta;
  double theta;
  theta = atan2(radius, fabs(z));
  eta = -log(tan(theta / 2.));
  if (z < 0)
  {
    eta = -eta;
  }
  return eta;
}

void PHG4Utils::SetColour(G4VisAttributes* att, const string& material)
{
  if (!att)
  {
    cout << "G4VisAttributes pointer is NULL" << endl;
    return;
  }
  if (material == "AL_BABAR_MAG")
  {
    att->SetColour(G4Colour::Blue());
  }
  else if (material == "BlackHole")
  {
    att->SetColour(G4Colour::Black());
  }
  else if (material == "C4F10")
  {
    att->SetColour(0., 0., 0.5, 0.25);
  }
  else if (material == "CF4")
  {
    att->SetColour(G4Colour::Magenta());
  }
  else if (material == "G4_AIR")
  {
    att->SetColour(G4Colour::Black());
  }
  else if (material == "G4_Al")
  {
    att->SetColour(G4Colour::Blue());
  }
  else if (material == "G4_Au")
  {
    att->SetColour(G4Colour::Yellow());
  }
  else if (material == "G4_CARBON_DIOXIDE")
  {
    att->SetColour(G4Colour::Green());
  }
  else if (material == "G4_CELLULOSE_CELLOPHANE")
  {
    att->SetColour(0.25, 0.25, 0.);
  }
  else if (material == "G4_Cu")
  {
    att->SetColour(1., 0.51, 0.278);
  }
  else if (material == "G4_Fe")
  {
    att->SetColour(0.29, 0.44, 0.54);
  }
  else if (material == "G4_KAPTON")
  {
    att->SetColour(G4Colour::Yellow());
  }
  else if (material == "G4_MYLAR")
  {
    att->SetColour(0.5, 0.5, 0.5, 0.25);
  }
  else if (material == "G4_METHANE")
  {
    att->SetColour(0., 1., 1., 0.25);
  }
  else if (material == "G4_Si")
  {
    att->SetColour(G4Colour::Yellow());
  }
  else if (material == "G4_TEFLON")
  {
    att->SetColour(G4Colour::White());
  }
  else if (material == "G4_W")
  {
    att->SetColour(0.36, 0.36, 0.36);
  }
  else if (material == "Quartz")
  {
    att->SetColour(G4Colour::Green());
  }
  else if (material == "Scintillator" || material == "G4_POLYSTYRENE")
  {
    att->SetColour(0., 1., 1.);
  }
  else if (material == "W_Epoxy")
  {
    att->SetColour(0.5, 0.5, 0.5);
  }
  else if (material == "G10")
  {
    att->SetColour(1., 1., 0., 0.5);
  }
  else
  {
    //cout << "default color red for material " << material << endl;
    att->SetColour(G4Colour::Cyan());
  }
  return;
}
