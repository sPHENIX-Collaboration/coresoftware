// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4ETAPARAMETERIZATION_H
#define G4MAIN_PHG4ETAPARAMETERIZATION_H

#include <Geant4/G4Types.hh>                // for G4int
#include <Geant4/G4VPVParameterisation.hh>

#include <iostream>                         // for cout, ostream
#include <vector>

class G4Tubs;
class G4VPhysicalVolume;

// Parameterization to define rings whose size changes with Z to 
// correspond to fixed width in eta.

class PHG4EtaParameterization : public G4VPVParameterisation
{ 
public:
  
  PHG4EtaParameterization(
			  unsigned int neta,       // Binning in eta
			  double minEta,  // "
			  double maxEta,  // "
			  double startPhi,
			  double deltaPhi,
			  double radiusIn,  // Radius of inner face of cylinder
			  double radiusOut,  // Radius of outer face of cylinder
			  double centerZ  // overall Z of center of rings
			  );
  
  ~PHG4EtaParameterization() override;
  
  virtual void Print(std::ostream& os = std::cout) const;

  void ComputeTransformation(const G4int copyNo,
			     G4VPhysicalVolume* physVol) const override;
  
  using G4VPVParameterisation::ComputeDimensions;// avoid warning for not implemented ComputeDimension methods
  void ComputeDimensions(G4Tubs& ring, const G4int copyNo,
			 const G4VPhysicalVolume* physVol) const override;

  int GetIEta(int copyNo) const { return _ieta.at(copyNo); }

private:  // Dummy declarations to get rid of warnings ...

//   void ComputeDimensions(G4Trd&,const G4int,const G4VPhysicalVolume*) const {}
//   void ComputeDimensions(G4Trap&,const G4int,const G4VPhysicalVolume*) const {}
//   void ComputeDimensions(G4Cons&,const G4int,const G4VPhysicalVolume*) const {}
//   void ComputeDimensions(G4Sphere&,const G4int,const G4VPhysicalVolume*) const {}
//   void ComputeDimensions(G4Orb&,const G4int,const G4VPhysicalVolume*) const {}
//   void ComputeDimensions(G4Torus&,const G4int,const G4VPhysicalVolume*) const {}
//   void ComputeDimensions(G4Para&,const G4int,const G4VPhysicalVolume*) const {}
//   void ComputeDimensions(G4Hype&,const G4int,const G4VPhysicalVolume*) const {}
//   void ComputeDimensions(G4Box&,const G4int,const G4VPhysicalVolume*) const {}
//   void ComputeDimensions(G4Polycone&,const G4int,const G4VPhysicalVolume*) const {}
//   void ComputeDimensions(G4Polyhedra&,const G4int,const G4VPhysicalVolume*) const {}

private:

  unsigned int _neta;   
  double _minEta;
  double _maxEta;
  double _startPhi;
  double _deltaPhi;
  double _radiusIn;
  double _radiusOut;
  double _centerZ;
  std::vector<double> _zpos; // Z positions of the rings
  std::vector<double> _zhalf; // Z half-widths of the rings
  std::vector<double> _phi0; // Lower edge of phi bins
  std::vector<double> _phi1; // Upper edge of phi bins
  std::vector<int> _ieta;
  std::vector<int> _iphi;
};


#endif
