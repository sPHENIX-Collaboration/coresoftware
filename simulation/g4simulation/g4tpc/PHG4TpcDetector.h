// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4TPC_PHG4TPCDETECTOR_H
#define G4TPC_PHG4TPCDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <set>
#include <string>

class G4LogicalVolume;
class G4UserLimits;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4TpcDisplayAction;
class PHG4Subsystem;
class PHParameters;

class PHG4TpcDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4TpcDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  ~PHG4TpcDetector(void) override
  {
  }

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  int IsInTpc(G4VPhysicalVolume *) const;
  void SuperDetector(const std::string &name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }

 private:
  int ConstructTpcGasVolume(G4LogicalVolume *tpc_envelope);
  int ConstructTpcCageVolume(G4LogicalVolume *tpc_envelope);
  PHG4TpcDisplayAction *m_DisplayAction;
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
