// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PHENIXDISPLAYACTION_H
#define G4MAIN_PHG4PHENIXDISPLAYACTION_H

#include "PHG4DisplayAction.h"

#include <map>
#include <string>
#include <vector>

class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;

class PHG4PhenixDisplayAction : public PHG4DisplayAction
{
 public:
  PHG4PhenixDisplayAction(const std::string &name);

  ~PHG4PhenixDisplayAction() override;

  void ApplyDisplayAction(G4VPhysicalVolume *physvol) override;
  void AddVolume(G4LogicalVolume *logvol, const std::string &mat) { m_LogicalVolumeMap[logvol] = mat; }

 private:
  std::map<G4LogicalVolume *, std::string> m_LogicalVolumeMap;
  std::vector<G4VisAttributes *> m_VisAttVec;
};

#endif  // G4MAIN_PHG4PHENIXDISPLAYACTION_H
