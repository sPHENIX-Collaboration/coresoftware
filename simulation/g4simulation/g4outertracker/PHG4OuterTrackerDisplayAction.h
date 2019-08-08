// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4MVTX_PHG4OUTERTRACKERDISPLAYACTION_H
#define G4MVTX_PHG4OUTERTRACKERDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <map>
#include <string>
#include <vector>

class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;

class PHG4OuterTrackerDisplayAction : public PHG4DisplayAction
{
 public:
  PHG4OuterTrackerDisplayAction(const std::string &name);

  virtual ~PHG4OuterTrackerDisplayAction();

  void ApplyDisplayAction(G4VPhysicalVolume *physvol);
  void AddVolume(G4LogicalVolume *logvol, const std::string &mat) { m_LogicalVolumeMap[logvol] = mat; }

 private:
  std::map<G4LogicalVolume *, std::string> m_LogicalVolumeMap;
  std::vector<G4VisAttributes *> m_VisAttVec;
};

#endif  // G4MVTX_PHG4MVTXDISPLAYACTION_H
