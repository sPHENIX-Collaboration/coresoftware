#ifndef PHG4AccSubsystem_h
#define PHG4AccSubsystem_h

#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

class PHG4AccDetector;
class PHG4AccSteppingAction;
class PHG4EventAction;

class PHG4AccSubsystem: public PHG4Subsystem
{

  public:

  //! constructor
  PHG4AccSubsystem( const std::string &name = "ACCCYLINDER", const int layer = 0 );

  //! destructor
  virtual ~PHG4AccSubsystem( void )
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

  //! accessors (reimplemented)
  virtual PHG4Detector* GetDetector( void ) const;
  virtual PHG4SteppingAction* GetSteppingAction( void ) const;
  PHG4EventAction* GetEventAction() const {return eventAction_;}

  void SetRadius(const G4double dbl) {radius = dbl;}
  void SetLength(const G4double dbl) {length = dbl;}
  void SetLengthViaRapidityCoverage(const G4bool bl) {lengthViaRapidityCoverage = bl;}
  void SetPosition(const G4double x, const G4double y, const G4double z)
  {
    xpos = x;
    ypos = y;
    zpos = z;
  }
  void SetTilt(const double tilt) {_sciTilt = tilt;}
  void SetScintWidth(const double wid) {_sciWidth = wid;}
  void SetTungstenWidth(const double wid) {_tungstenWidth = wid;}
  void SetUndulMag(const double unmag) {_unmag = unmag;}
  void SetNumScint(const int num) {_sciNum = num;}
  void SetThickness(const G4double dbl) {TrackerThickness = dbl;}
  void SetMaterial(const std::string &mat) {material = mat;}
  void SetActive(const int i = 1) {active = i;}
  void SetAbsorberActive(const int i = 1) {absorberactive = i;}
  void SuperDetector(const std::string &name) {superdetector = name;}
  void UseDrawing(const int i=1) {usedrawing = i;}
  const std::string SuperDetector() {return superdetector;}

  void Print(const std::string &what = "ALL") const;

  private:

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4AccDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4AccSteppingAction* steppingAction_;
  PHG4EventAction *eventAction_;
  G4double radius;
  G4double length;
  G4double xpos,ypos,zpos;
  G4bool lengthViaRapidityCoverage;
  G4double TrackerThickness;
  G4String material;
  G4double _sciTilt;
  G4double _sciWidth;
  G4double _tungstenWidth;
  G4double _unmag;
  G4int    _sciNum;
  int active;
  int absorberactive;
  int layer;
  std::string detector_type;
  std::string superdetector;
  int usedrawing;
};

#endif
