// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4EPD_PHG4EPDDISPLAYACTION_H
#define G4EPD_PHG4EPDDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <map>
#include <string>  // for string
#include <vector>

class G4LogicalVolume;
class G4VisAttributes;
class G4VPhysicalVolume;

class PHG4EPDDisplayAction : public PHG4DisplayAction
{
 public:
  PHG4EPDDisplayAction(const std::string &name);

  ~PHG4EPDDisplayAction() override;

  void ApplyDisplayAction(G4VPhysicalVolume *physvol) override;
  void AddVolume(G4LogicalVolume *logvol, const std::string &mat) { m_LogicalVolumeMap[logvol] = mat; }

 private:
  std::map<G4LogicalVolume *, std::string> m_LogicalVolumeMap;
  std::vector<G4VisAttributes *> m_VisAttVec;
};

#endif  // G4EPD_PHG4EPDDISPLAYACTION_H
