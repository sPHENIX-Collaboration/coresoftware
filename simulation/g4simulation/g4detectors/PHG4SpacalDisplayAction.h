// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SPACALDISPLAYACTION_H
#define G4DETECTORS_PHG4SPACALDISPLAYACTION_H

#include "PHG4CylinderGeom_Spacalv3.h"

#include <g4main/PHG4DisplayAction.h>

#include <map>
#include <string>
#include <vector>

class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;

class PHG4SpacalDisplayAction : public PHG4DisplayAction
{
 public:
  PHG4SpacalDisplayAction(const std::string &name);

  virtual ~PHG4SpacalDisplayAction();

  void ApplyDisplayAction(G4VPhysicalVolume *physvol);
  void AddVolume(G4LogicalVolume *logvol, const std::string &mat) { m_LogicalVolumeMap[logvol] = mat; }
  void SetGeom(const PHG4CylinderGeom_Spacalv3 *geo) {m_GeomV3 = geo;}
  void AddMaterial(const std::string &name, const std::string &mat) {m_MaterialMap[name] = mat;}

 private:
  const PHG4CylinderGeom_Spacalv3 *m_GeomV3;
  std::map<G4LogicalVolume *, std::string> m_LogicalVolumeMap;
  std::vector<G4VisAttributes *> m_VisAttVec;
  std::map<std::string, std::string> m_MaterialMap;
};

#endif  // G4DETECTORS_PHG4TPCDISPLAYACTION_H
