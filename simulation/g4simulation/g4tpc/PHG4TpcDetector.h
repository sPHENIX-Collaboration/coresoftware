// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4TPC_PHG4TPCDETECTOR_H
#define G4TPC_PHG4TPCDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <cmath>
#include <set>
#include <string>
#include <vector>

class G4LogicalVolume;
class G4UserLimits;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4TpcDisplayAction;
class PHG4Subsystem;
class PHParameters;

class PHG4TpcDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4TpcDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  ~PHG4TpcDetector(void) override
  {
  }

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  int IsInTpc(G4VPhysicalVolume *) const;
  void SuperDetector(const std::string &name) { m_SuperDetectorName = name; }
  const std::string SuperDetector() const { return m_SuperDetectorName; }

 private:
  int ConstructTpcGasVolume(G4LogicalVolume *tpc_envelope);
  int ConstructTpcCageVolume(G4LogicalVolume *tpc_envelope);
  int ConstructTpcExternalSupports(G4LogicalVolume *logicWorld);

  void CreateCompositeMaterial(std::string compositeName, std::vector<std::string> materialName, std::vector<double> thickness);

  PHG4TpcDisplayAction *m_DisplayAction = nullptr;
  PHParameters *m_Params = nullptr;
  G4UserLimits *m_G4UserLimits = nullptr;
  int m_ActiveFlag = 0;
  int m_AbsorberActiveFlag = 0;
  double m_InnerCageRadius = NAN;
  double m_OuterCageRadius = NAN;
  std::set<G4VPhysicalVolume *> m_AbsorberVolumeSet;
  std::set<G4VPhysicalVolume *> m_ActiveVolumeSet;

  std::string m_SuperDetectorName;
};

#endif
