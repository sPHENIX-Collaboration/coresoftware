#ifndef G4DETECTORS_PHG4FPBSCSUBSYSTEM_H
#define G4DETECTORS_PHG4FPBSCSUBSYSTEM_H

#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

class PHG4FPbScDetector;
class PHG4FPbScSteppingAction;
class PHG4EventAction;

class PHG4FPbScSubsystem: public PHG4Subsystem
{
  
  public:
    
    //! constructor
  PHG4FPbScSubsystem( const std::string &name = "FPBSC" );
    
    //! destructor
    virtual ~PHG4FPbScSubsystem( void )
    {}
    
    //! init
    /*!
    creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
    reates the stepping action and place it on the node tree, under "ACTIONS" node
    creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
    */
    int Init(PHCompositeNode *);
    
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

    void set_Place(double x, double y, double z)
    {
      x_position = x;
      y_position = y;
      z_position = z;
    }
   
  private:
    
    //! detector geometry
    /*! defives from PHG4Detector */
    PHG4FPbScDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
    PHG4FPbScSteppingAction* steppingAction_;
    PHG4EventAction *eventAction_;
    
    double x_position;
    double y_position;
    double z_position;
    
};

#endif
