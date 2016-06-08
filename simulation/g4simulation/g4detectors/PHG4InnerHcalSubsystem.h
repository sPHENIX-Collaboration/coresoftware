#ifndef PHG4InnerHcalSubsystem_h
#define PHG4InnerHcalSubsystem_h

#include "PHG4DetectorSubsystem.h"

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

#include <map>
#include <set>
#include <string>

class PHG4InnerHcalDetector;
class PHG4Parameters;
class PHG4InnerHcalSteppingAction;
class PHG4EventAction;
class PHG4FlushStepTrackingAction;

class PHG4InnerHcalSubsystem: public PHG4DetectorSubsystem
{

  public:

  //! constructor
  PHG4InnerHcalSubsystem( const std::string &name = "HCALIN", const int layer = 0 );

  //! destructor
  virtual ~PHG4InnerHcalSubsystem( void )
  {}

  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode *);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode *);

  //! Print info (from SubsysReco)
  void Print(const std::string &what = "ALL") const;

  //! accessors (reimplemented)
  PHG4Detector* GetDetector( void ) const;
  PHG4SteppingAction* GetSteppingAction( void ) const {return steppingAction_;}
  PHG4TrackingAction* GetTrackingAction( void ) const {return  trackingAction_;}

  PHG4EventAction* GetEventAction() const {return eventAction_;}

  void SetLightCorrection(const double inner_radius, const double inner_corr,const double outer_radius, const double outer_corr);


  private:

  void SetDefaultParameters();

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4InnerHcalDetector* detector_;

  //! detector "stepping" action, executes after every G4 step
  /*! derives from PHG4SteppingAction */
  PHG4SteppingAction* steppingAction_;

  //! detector "tracking" action, executes before/after every G4 track
  /*! derives from PHG4TrackingAction */
  PHG4TrackingAction *trackingAction_;

  //! detector event action executes before/after every event
  /*! derives from PHG4EventAction */
  PHG4EventAction *eventAction_;
};

#endif
