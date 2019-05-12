// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4TPC_PHG4TPCDETECTOR_H
#define G4TPC_PHG4TPCDETECTOR_H

#include <g4main/PHG4Detector.h>

// cannot fwd declare G4RotationMatrix, it is a typedef pointing to clhep
#include <Geant4/G4RotationMatrix.hh>

#include <set>
#include <string>

class G4LogicalVolume;
class G4UserLimits;
class G4VPhysicalVolume;
class G4VSolid;
class PHG4TPCDisplayAction;
class PHG4TPCSubsystem;
class PHParameters;

class PHG4TPCDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4TPCDetector(PHG4TPCSubsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~PHG4TPCDetector(void)
  {
  }

  //! construct
  void Construct(G4LogicalVolume *world);

  int IsInTPC(G4VPhysicalVolume *) const;
  void SuperDetector(const std::string &name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }

 private:
  int ConstructTPCGasVolume(G4LogicalVolume *tpc_envelope);
  int ConstructTPCCageVolume(G4LogicalVolume *tpc_envelope);
  PHG4TPCDisplayAction *m_DisplayAction;
  PHParameters *params;
  G4UserLimits *g4userlimits;
  int active;
  int absorberactive;
  double inner_cage_radius;
  double outer_cage_radius;
  std::set<G4VPhysicalVolume *> absorbervols;
  std::set<G4VPhysicalVolume *> activevols;

  std::string superdetector;
};

#endif
