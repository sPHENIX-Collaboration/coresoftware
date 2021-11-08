// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4TPC_PHG4TPCDISPLAYACTION_H
#define G4TPC_PHG4TPCDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <map>
#include <string>
#include <vector>

class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;

class PHG4TpcDisplayAction : public PHG4DisplayAction
{
 public:
  PHG4TpcDisplayAction(const std::string &name);

  ~PHG4TpcDisplayAction() override;

  void ApplyDisplayAction(G4VPhysicalVolume *physvol) override;
  void AddVolume(G4LogicalVolume *logvol, const std::string &mat) { m_LogicalVolumeMap[logvol] = mat; }
  void AddTpcInnerLayer(G4LogicalVolume *logvol) { m_TpcInnerLayersVec.push_back(logvol); }
  void AddTpcOuterLayer(G4LogicalVolume *logvol) { m_TpcOuterLayersVec.push_back(logvol); }

 private:
  std::map<G4LogicalVolume *, std::string> m_LogicalVolumeMap;
  std::vector<G4VisAttributes *> m_VisAttVec;
  std::vector<G4LogicalVolume *> m_TpcInnerLayersVec;
  std::vector<G4LogicalVolume *> m_TpcOuterLayersVec;
};

#endif  // G4TPC_PHG4TPCDISPLAYACTION_H
