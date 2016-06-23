#ifndef PHG4CylinderSubsystem_h
#define PHG4CylinderSubsystem_h

#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

#include <float.h>

class PHG4CylinderDetector;
class PHG4CylinderSteppingAction;
class PHG4EventAction;
class PHG4FlushStepTrackingAction;

class PHG4CylinderSubsystem: public PHG4Subsystem
{

  public:

  //! constructor
  PHG4CylinderSubsystem( const std::string &name = "CYLINDER", const int layer = 0 );

  //! destructor
  virtual ~PHG4CylinderSubsystem( void )
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
  PHG4Detector* GetDetector( void ) const;
  PHG4SteppingAction* GetSteppingAction( void ) const;
  PHG4TrackingAction* GetTrackingAction( void ) const;
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
  void SetThickness(const G4double dbl) {TrackerThickness = dbl;}
  void SetMaterial(const std::string &mat) {material = mat;}
  void SetActive(const int i = 1) {active = i;}
  void SetReducedTruthInfo(const bool b) {reduced = b;}
  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() {return superdetector;}

   void BlackHole(const int i=1,
		 const double tmin = DBL_MAX,
		 const double tmax = -1.0*DBL_MAX) {
    blackhole = i; blackhole_tmin = tmin; blackhole_tmax = tmax;
  }
  int IsBlackHole() const {return blackhole;}
  double GetBlackHoleTMin() {return blackhole_tmin;}
  double GetBlackHoleTMax() {return blackhole_tmax;}
  
  private:

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4CylinderDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4CylinderSteppingAction* steppingAction_;

  PHG4FlushStepTrackingAction *trackingAction_;

  PHG4EventAction *eventAction_;
  G4double radius;
  G4double length;
  G4double xpos,ypos,zpos;
  G4bool lengthViaRapidityCoverage;
  G4double TrackerThickness;
  G4String material;
  int active;
  bool reduced;
  int layer;
  int blackhole;
  double blackhole_tmin;
  double blackhole_tmax;
  std::string detector_type;
  std::string superdetector;
};

#endif
