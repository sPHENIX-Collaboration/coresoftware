// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BEAMLINEMAGNETDISPLAYACTION_H
#define BEAMLINEMAGNETDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <map>
#include <string>
#include <vector>

class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;

class BeamLineMagnetDisplayAction : public PHG4DisplayAction
{
 public:
  BeamLineMagnetDisplayAction(const std::string &name);

  virtual ~BeamLineMagnetDisplayAction();

  void ApplyDisplayAction(G4VPhysicalVolume *physvol);
  void AddVolume(G4LogicalVolume *logvol, const std::string &mat) { m_LogicalVolumeMap[logvol] = mat; }

 private:
  std::map<G4LogicalVolume *, std::string> m_LogicalVolumeMap;
  std::vector<G4VisAttributes *> m_VisAttVec;
};

#endif  // BEAMLINEMAGNETDISPLAYACTION_H
