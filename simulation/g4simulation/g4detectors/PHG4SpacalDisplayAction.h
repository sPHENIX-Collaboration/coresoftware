// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SPACALDISPLAYACTION_H
#define G4DETECTORS_PHG4SPACALDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <map>
#include <string>
#include <vector>

class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;
class PHG4CylinderGeom_Spacalv1;

class PHG4SpacalDisplayAction : public PHG4DisplayAction
{
 public:
  explicit PHG4SpacalDisplayAction(const std::string &name);

  ~PHG4SpacalDisplayAction() override;

  void ApplyDisplayAction(G4VPhysicalVolume *physvol) override;
  void AddVolume(G4LogicalVolume *logvol, const std::string &mat) { m_LogicalVolumeMap[logvol] = mat; }
  void SetGeom(const PHG4CylinderGeom_Spacalv1 *geo) { m_Geom = geo; }
  void AddMaterial(const std::string &name, const std::string &mat) { m_MaterialMap[name] = mat; }

 private:
  const PHG4CylinderGeom_Spacalv1 *m_Geom = nullptr;
  std::map<G4LogicalVolume *, std::string> m_LogicalVolumeMap;
  std::vector<G4VisAttributes *> m_VisAttVec;
  std::map<std::string, std::string> m_MaterialMap;
};

#endif  // G4DETECTORS_PHG4TPCDISPLAYACTION_H
