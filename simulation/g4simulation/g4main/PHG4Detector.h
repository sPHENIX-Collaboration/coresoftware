// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4DETECTOR_H
#define G4MAIN_PHG4DETECTOR_H

#include <Geant4/G4RotationMatrix.hh>

#include <iostream>
#include <string>

class G4LogicalVolume;
class G4Material;
class G4Element;
class G4UserSteppingAction;
class G4VSolid;
class PHCompositeNode;
class PHG4Subsystem;

//! base class for phenix detector creation
/*! derived classes must implement construct method, which takes the "world" logical volume as argument */
class PHG4Detector
{
 public:
  //! constructor
  // delete default ctor, nobody should use it
  PHG4Detector() = delete;
  // this is the ctor we use
  explicit PHG4Detector(PHG4Subsystem *subsys, PHCompositeNode *Node, const std::string &nam);

  //! destructor
  virtual ~PHG4Detector(void)
  {
  }

  //! construct method
  /*!
  construct all logical and physical volumes relevant for given detector and place them
  inside the world logical volume
  */
  virtual void Construct(G4LogicalVolume *world) final;

  virtual void ConstructMe(G4LogicalVolume *mothervolume) = 0;

  //! Optional PostConstruction call after all geometry is constructed
  virtual void PostConstruction() {};

  virtual void Verbosity(const int v) { m_Verbosity = v; }

  virtual int Verbosity() const { return m_Verbosity; }
  virtual G4UserSteppingAction *GetSteppingAction() { return nullptr; }
  virtual std::string GetName() const { return m_Name; }
  virtual void OverlapCheck(const bool chk) { m_OverlapCheck = chk; }
  virtual bool OverlapCheck() const { return m_OverlapCheck; }
  virtual void Print(const std::string &/*what*/ = "ALL") const
  {
    std::cout << GetName() << ": Print method not implemented" << std::endl;
  }
  virtual int DisplayVolume(G4VSolid *volume, G4LogicalVolume *logvol, G4RotationMatrix *rotm = nullptr);
  virtual int DisplayVolume(G4LogicalVolume *checksolid, G4LogicalVolume *logvol, G4RotationMatrix *rotm = nullptr);
  virtual PHCompositeNode *topNode() { return m_topNode; }
  virtual PHG4Subsystem *GetMySubsystem() {return m_MySubsystem;}
  static G4Material *GetDetectorMaterial(const std::string &name, const bool quit = true);
  static G4Element *GetDetectorElement(const std::string &name, const bool quit = true);

 private:
  PHCompositeNode *m_topNode = nullptr;
  PHG4Subsystem *m_MySubsystem = nullptr;
  int m_Verbosity = 0;
  bool m_OverlapCheck = false;
  int m_ColorIndex = 0;
  std::string m_Name;
};

#endif  // G4MAIN_PHG4DETECTOR_H
