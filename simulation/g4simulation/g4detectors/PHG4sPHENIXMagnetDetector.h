// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SPHENIXMAGNETDETECTOR_H
#define G4DETECTORS_PHG4SPHENIXMAGNETDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <set>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4sPHENIXMagnetDisplayAction;
class PHG4Subsystem;
class PHG4GDMLConfig;
class PHParameters;
class G4Box;
class G4Polycone;
class G4Tubs;


/**
 */

class PHG4sPHENIXMagnetDetector : public PHG4Detector
{
 public:
  //! constructor
  explicit PHG4sPHENIXMagnetDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int detid);

  //! destructor
  ~PHG4sPHENIXMagnetDetector() override = default;

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  //!@name volume accessors
  int IsInsPHENIXMagnet(G4VPhysicalVolume *) const;

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string &SuperDetector() const { return m_SuperDetector; }

  int get_Layer() const { return m_Layer; }

  PHG4sPHENIXMagnetDisplayAction *GetDisplayAction() { return m_DisplayAction; }

 private:
  PHParameters *GetParams() const { return m_Params; }

  PHG4sPHENIXMagnetDisplayAction *m_DisplayAction = nullptr;
  PHParameters *m_Params = nullptr;
  //! registry for volumes that should not be exported, i.e. fibers
  PHG4GDMLConfig *m_GdmlConfig = nullptr;

  G4Box* Block(G4int iblock) const;
  G4Tubs* CryoTubes(G4int itube) const;
  G4Polycone* SolenoidPolycones(G4int ipolycone) const;
  G4Tubs* SolenoidTubes(G4int itube) const;


  int m_ActiveFlag;
  int m_Layer;

  std::string m_SuperDetector;

  std::set<G4LogicalVolume *> m_LogicalVolSet;

};

#endif
