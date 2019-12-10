// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CRYSTALCALORIMETERDISPLAYACTION_H
#define G4DETECTORS_PHG4CRYSTALCALORIMETERDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <map>
#include <string>
#include <vector>

class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;

class PHG4CrystalCalorimeterDisplayAction : public PHG4DisplayAction
{
 public:
  PHG4CrystalCalorimeterDisplayAction(const std::string &name);

  virtual ~PHG4CrystalCalorimeterDisplayAction();

  void ApplyDisplayAction(G4VPhysicalVolume *physvol);
  void AddVolume(G4LogicalVolume *logvol, const std::string &mat) { m_LogicalVolumeMap[logvol] = mat; }

 private:
  std::map<G4LogicalVolume *, std::string> m_LogicalVolumeMap;
  std::vector<G4VisAttributes *> m_VisAttVec;
};

#endif  // G4DETECTORS_PHG4CRYSTALCALORIMETERDISPLAYACTION_H
