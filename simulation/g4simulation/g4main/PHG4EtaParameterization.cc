#include "PHG4EtaParameterization.h"

#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>  // for G4int
#include <Geant4/G4VPhysicalVolume.hh>

#include <algorithm>  // for copy
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>

using namespace std;

PHG4EtaParameterization::PHG4EtaParameterization(
    unsigned int neta,  // Binning in eta
    double minEta,      // "
    double maxEta,      // "
    double startPhi,
    double deltaPhi,
    double radiusIn,   // Radius of inner face of cylinder
    double radiusOut,  // Radius of outer face of cylinder
    double centerZ     // Z of center of set
    )
  : _neta(neta)
  , _minEta(minEta)
  , _maxEta(maxEta)
  , _startPhi(startPhi)
  , _deltaPhi(deltaPhi)
  , _radiusIn(radiusIn)
  , _radiusOut(radiusOut)
  , _centerZ(centerZ)
{
  if (_maxEta < _minEta)
  {
    cout << " invalid eta, max<min"
         << " etamin: " << _minEta
         << " etamax: " << _maxEta
         << endl;
    exit(1);
  }
  //    G4Exception("PHG4EtaParameterization::PHG4EtaParameterization", "invalid eta, max<min",G4ExceptionSeverity::FatalException);

  if ((_radiusIn < 0.0) || (_radiusOut < 0.0) || (_radiusOut < _radiusIn))
  {
    cout << " invalid radius parameters:"
         << " radiusIn: " << radiusIn
         << " radiusOut: " << radiusOut
         << endl;
    exit(1);
  }
  //    G4Exception("PHG4EtaParameterization::PHG4EtaParameterization: invalid radius parameters");

  double totalEta = _maxEta - _minEta;
  double dEta = totalEta / _neta;
  //double minZ = 1e6;
  //double maxZ = -1e6;
  for (unsigned int i = 0; i < neta; i++)
  {
    // Compute the edges of this eta bin
    double etaMin = _minEta + dEta * i;
    double etaMax = etaMin + dEta;
    // Compute the corresponding Z positions of the edges
    double zmin = _centerZ + _radiusIn * std::sinh(etaMin);
    double zmax = _centerZ + _radiusIn * std::sinh(etaMax);
    // Z positions is halfway between the edges
    double zpos = (zmin + zmax) / 2.0;
    double zhalf = (zmax - zmin) / 2.0;
    _zpos.push_back(zpos);
    _zhalf.push_back(zhalf);
    std::cout << zmin << " " << zmax << " " << zpos << " +/- " << zhalf << std::endl;
  }

  // Build lookup vectors for the copyNo->(ieta, iphi) translation
  //
  for (unsigned int i = 0; i < _neta; i++)
  {
    _ieta.push_back(i);
  }

  std::cout << "*********** Constructing PHG4EtaParameterization ***************" << std::endl;
  std::cout << std::endl;

  std::cout << "Radii = " << _radiusIn << ", " << _radiusOut << std::endl;
  std::cout << "Phi,dPhi = " << _startPhi << ", " << _deltaPhi << std::endl;
  std::cout << "Min/Max Z = " << _zpos.front() - _zhalf.front() << " / " << _zpos.back() + _zhalf.back() << std::endl;

  std::cout << std::endl;
  std::cout << "********* End Constructing PHG4EtaParameterization *************" << std::endl;
}

PHG4EtaParameterization::~PHG4EtaParameterization()
{
  std::cout << "PHG4EtaParameterization::~PHG4EtaParameterization: Alas, poor Yorick! I knew him, Horatio"
            << std::endl;
}

void PHG4EtaParameterization::Print(std::ostream& os) const
{
  os << "PhiEtaParameterization: NETA = " << _neta << std::endl;
  os << "Zpos: ";
  std::copy(_zpos.begin(), _zpos.end(), std::ostream_iterator<double>(os, " "));
  os << std::endl;
}

void PHG4EtaParameterization::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  int iring = copyNo;
  G4ThreeVector origin(0, 0, _zpos.at(iring));
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

void PHG4EtaParameterization::ComputeDimensions(G4Tubs& ring, const G4int copyNo,
                                                const G4VPhysicalVolume*) const
{
  //int ieta = GetIEta(copyNo);
  ring.SetZHalfLength(_zhalf.at(copyNo));
  ring.SetInnerRadius(_radiusIn);
  ring.SetOuterRadius(_radiusOut);
  ring.SetStartPhiAngle(_startPhi);
  ring.SetDeltaPhiAngle(_deltaPhi);
}
