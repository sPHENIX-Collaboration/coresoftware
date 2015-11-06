#ifndef PHG4InnerHcalSubsystem_h
#define PHG4InnerHcalSubsystem_h

#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

class PHG4InnerHcalDetector;
class PHG4InnerHcalParameters;
class PHG4InnerHcalSteppingAction;
class PHG4EventAction;

class PHG4InnerHcalSubsystem: public PHG4Subsystem
{

  public:

  //! constructor
  PHG4InnerHcalSubsystem( const std::string &name = "HCALIN", const int layer = 0 );

  //! destructor
  virtual ~PHG4InnerHcalSubsystem( void )
  {}

  //! init
  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRun(PHCompositeNode *);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode *);

  //! Print info (from SubsysReco)
  void Print(const std::string &what = "ALL") const;

  //! accessors (reimplemented)
  virtual PHG4Detector* GetDetector( void ) const;
  virtual PHG4SteppingAction* GetSteppingAction( void ) const;

  void SetPlaceZ(const G4double dbl);
  void SetPlace(const G4double place_x, const G4double place_y, const G4double place_z);

  void SetXRot(const G4double dbl);
  void SetYRot(const G4double dbl);
  void SetZRot(const G4double dbl);
  PHG4EventAction* GetEventAction() const {return eventAction_;}
  void SetActive(const int i = 1);
  void SetAbsorberActive(const int i = 1);
  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() {return superdetector;}

  void BlackHole(const int i=1);
  PHG4InnerHcalParameters *GetParameters();

  void SetTiltViaNcross(const int ncross);
  void SetInnerRadius(const double inner);
  double GetInnerRadius() const;
  void SetOuterRadius(const double outer);
  double GetOuterRadius() const;
  void SetStepLimits(const double slim);

  void SetLightCorrection(const float inner_radius, const float inner_corr,
			  const float outer_radius, const float outer_corr);
  void SetLightScintModel(const bool b = true);


  protected:

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4InnerHcalDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingAction */
  PHG4InnerHcalSteppingAction* steppingAction_;

  //! particle tracking "stepping" action
  /*! derives from PHG4EventAction */
  PHG4EventAction *eventAction_;

  PHG4InnerHcalParameters *params;

  int layer;
  std::string detector_type;
  std::string superdetector;
};

#endif
