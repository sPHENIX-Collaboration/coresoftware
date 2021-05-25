// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4HCALSUBSYSTEM_H
#define G4DETECTORS_PHG4HCALSUBSYSTEM_H

#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4String.hh>
#include <Geant4/G4Types.hh>

#include <string>  // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4HcalDetector;
class PHG4HcalSteppingAction;
class PHG4SteppingAction;

class PHG4HcalSubsystem : public PHG4Subsystem
{
 public:
  //! constructor
  PHG4HcalSubsystem(const std::string &name = "HCALCYLINDER", const int layer = 0);

  //! destructor
  ~PHG4HcalSubsystem(void) override
  {
  }

  //! init
  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRun(PHCompositeNode *) override;

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode *) override;

  //! accessors (reimplemented)
  PHG4Detector *GetDetector(void) const override;
  PHG4SteppingAction *GetSteppingAction(void) const override;

  void SetRadius(const G4double dbl) { radius = dbl; }
  void SetLength(const G4double dbl) { length = dbl; }
  void SetLengthViaRapidityCoverage(const G4bool bl) { lengthViaRapidityCoverage = bl; }
  void SetPosition(const G4double x, const G4double y, const G4double z)
  {
    xpos = x;
    ypos = y;
    zpos = z;
  }
  void SetTilt(const double tilt) { _sciTilt = tilt; }
  void SetTiltViaNcross(const int ncross);
  void SetScintWidth(const double wid) { _sciWidth = wid; }
  void SetNumScint(const int num) { _sciNum = num; }
  void SetScintPhi0(const G4double phi0) { _sciPhi0 = phi0; }  //  in units of sampling cells
  void SetThickness(const G4double dbl) { TrackerThickness = dbl; }
  void SetMaterial(const std::string &mat) { material = mat; }
  void SetActive(const int i = 1) { active = i; }
  void SetAbsorberActive(const int i = 1) { absorberactive = i; }
  void SuperDetector(const std::string &name) { superdetector = name; }
  const std::string SuperDetector() { return superdetector; }

  void SetLightCorrection(float inner_radius, float inner_corr,
                          float outer_radius, float outer_corr)
  {
    light_balance_ = true;
    light_balance_inner_radius_ = inner_radius;
    light_balance_inner_corr_ = inner_corr;
    light_balance_outer_radius_ = outer_radius;
    light_balance_outer_corr_ = outer_corr;
  }
  void SetLightScintModel(const bool b = true)
  {
    light_scint_model_ = b;
  }

  void Print(const std::string &what = "ALL") const override;

 private:
  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4HcalDetector *detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4HcalSteppingAction *steppingAction_;
  G4double radius;
  G4double length;
  G4double xpos, ypos, zpos;
  G4bool lengthViaRapidityCoverage;
  G4double TrackerThickness;
  G4String material;
  G4double _sciTilt;
  G4double _sciWidth;
  G4int _sciNum;
  G4double _sciPhi0;  //  in units of sampling cells
  int active;
  int absorberactive;
  int layer;
  std::string detector_type;
  std::string superdetector;

  bool light_scint_model_;
  bool light_balance_;
  float light_balance_inner_radius_;
  float light_balance_inner_corr_;
  float light_balance_outer_radius_;
  float light_balance_outer_corr_;
};

#endif
