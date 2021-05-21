#include "PHG4EtaPhiParameterization.h"

#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>  // for G4int
#include <Geant4/G4VPhysicalVolume.hh>

#include <algorithm>  // for copy
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>

PHG4EtaPhiParameterization::PHG4EtaPhiParameterization(
    unsigned int neta,  // Binning in eta
    double minEta,      // "
    double maxEta,      // "
    unsigned int nphi,  // number of phi bins
    double startPhi,
    double deltaPhi,
    double radiusIn,   // Radius of inner face of cylinder
    double radiusOut,  // Radius of outer face of cylinder
    double centerZ     // Z of center of set
    )
  : _neta(neta)
  , _minEta(minEta)
  , _maxEta(maxEta)
  , _nphi(nphi)
  , _startPhi(startPhi)
  , _deltaPhi(deltaPhi)
  , _radiusIn(radiusIn)
  , _radiusOut(radiusOut)
  , _centerZ(centerZ)
{
  if (_maxEta < _minEta)
  {
    std::cout << " invalid eta, max<min"
              << " etamin: " << _minEta
              << " etamax: " << _maxEta
              << std::endl;
    exit(1);
  }

  if ((_radiusIn < 0.0) || (_radiusOut < 0.0) || (_radiusOut < _radiusIn))
  {
    std::cout << " invalid radius parameters:"
              << " radiusIn: " << radiusIn
              << " radiusOut: " << radiusOut
              << std::endl;
    exit(1);
  }

  double totalEta = _maxEta - _minEta;
  double dEta = totalEta / _neta;
  for (unsigned int i = 0; i < neta; i++)
  {
    // Compute the edges of this eta bin
    double etaMin = _minEta + dEta * i;
    double etaMax = etaMin + dEta;
    //double eta = _minEta + dEta * (i + 0.5);
    // Compute the corresponding Z positions of the edges
    double zmin = _centerZ + _radiusIn * std::sinh(etaMin);
    double zmax = _centerZ + _radiusIn * std::sinh(etaMax);
    // Z positions is halfway between the edges
    double zpos = (zmin + zmax) / 2.0;
    double zhalf = (zmax - zmin) / 2.0;
    _zpos.push_back(zpos);
    _zhalf.push_back(zhalf);
  }

  for (unsigned int i = 0; i < _nphi; i++)
  {
    _phi0.push_back(_startPhi + i * _deltaPhi);
    _phi1.push_back(_startPhi + (i + 1) * _deltaPhi);
  }

  // Build lookup vectors for the copyNo->(ieta, iphi) translation
  //
  for (unsigned int i = 0; i < _neta * _nphi; i++)
  {
    div_t q = div((int)i, (int)_nphi);
    int ieta = q.quot;
    int iphi = q.rem;
    _ieta.push_back(ieta);
    _iphi.push_back(iphi);
  }
}

PHG4EtaPhiParameterization::~PHG4EtaPhiParameterization()
{
  std::cout << "PHG4EtaPhiParameterization::~PHG4EtaPhiParameterization: Alas, poor Yorick! I knew him, Horatio"
            << std::endl;
}

void PHG4EtaPhiParameterization::Print(std::ostream& os) const
{
  os << "PhiEtaPhiParameterization: NETA x NPHI = " << _neta << " x " << _nphi << std::endl;
  os << "Zpos: ";
  std::copy(_zpos.begin(), _zpos.end(), std::ostream_iterator<double>(os, " "));
  os << std::endl;
  os << "Phi0: ";
  std::copy(_phi0.begin(), _phi0.end(), std::ostream_iterator<double>(os, " "));
  os << std::endl;
  os << "Phi1: ";
  std::copy(_phi0.begin(), _phi0.end(), std::ostream_iterator<double>(os, " "));
  os << std::endl;
}

void PHG4EtaPhiParameterization::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  int iring = copyNo / _nphi;
  G4ThreeVector origin(0, 0, _zpos.at(iring));
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

void PHG4EtaPhiParameterization::ComputeDimensions(G4Tubs& ring, const G4int copyNo,
                                                   const G4VPhysicalVolume*) const
{
  int ieta = GetIEta(copyNo);
  int iphi = GetIPhi(copyNo);
  double phi = _phi0.at(iphi);
  double dphi = _phi1.at(iphi) - phi;
  ring.SetZHalfLength(_zhalf.at(ieta));
  ring.SetInnerRadius(_radiusIn);
  ring.SetOuterRadius(_radiusOut);
  ring.SetStartPhiAngle(phi);
  ring.SetDeltaPhiAngle(dphi);
}
