// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4JLEIC_G4JLEICBEAMLINEMAGNETDISPLAYACTION_H
#define G4JLEIC_G4JLEICBEAMLINEMAGNETDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <map>
#include <string>
#include <vector>

class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;
class PHParameters;

class G4JLeicBeamLineMagnetDisplayAction : public PHG4DisplayAction
{
 public:
  G4JLeicBeamLineMagnetDisplayAction(const std::string &name, PHParameters *parameters);

  virtual ~G4JLeicBeamLineMagnetDisplayAction();

  void ApplyDisplayAction(G4VPhysicalVolume *physvol);
  void AddVolume(G4LogicalVolume *logvol, const std::string &mat) { m_LogicalVolumeMap[logvol] = mat; }

 private:
  PHParameters *m_Params;
  std::map<G4LogicalVolume *, std::string> m_LogicalVolumeMap;
  std::vector<G4VisAttributes *> m_VisAttVec;
};

#endif  // G4JLEIC_G4JLEICBEAMLINEMAGNETDISPLAYACTION_H
