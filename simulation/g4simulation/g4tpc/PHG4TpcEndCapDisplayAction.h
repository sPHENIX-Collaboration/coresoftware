// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4TpcEndCapDisplayAction_H
#define G4DETECTORS_PHG4TpcEndCapDisplayAction_H

#include <g4main/PHG4DisplayAction.h>

#include <map>
#include <string>
#include <vector>

class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;

class PHG4TpcEndCapDisplayAction : public PHG4DisplayAction
{
 public:
  explicit PHG4TpcEndCapDisplayAction(const std::string &name);

  ~PHG4TpcEndCapDisplayAction() override;

  void ApplyDisplayAction(G4VPhysicalVolume *physvol) override;
  void AddVolume(G4LogicalVolume *logvol, const std::string &mat) { m_LogicalVolumeMap[logvol] = mat; }

 private:
  std::map<G4LogicalVolume *, std::string> m_LogicalVolumeMap;
  std::vector<G4VisAttributes *> m_VisAttVec;
};

#endif  // G4DETECTORS_PHG4TpcEndCapDisplayAction_H
