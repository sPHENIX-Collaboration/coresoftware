#ifndef PHG4SectorSubsystem_h
#define PHG4SectorSubsystem_h


#include "PHG4SectorConstructor.h"

#include <g4main/PHG4Subsystem.h>

class PHG4SectorDetector;
class PHG4SectorSteppingAction;
class PHG4EventAction;

class PHG4SectorSubsystem : public PHG4Subsystem
{

public:

  //! constructor
  PHG4SectorSubsystem(const std::string &name = "Sector");

  //! destructor
  virtual
  ~PHG4SectorSubsystem(void)
  {
  }

  //! init
  /*!
   creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
   reates the stepping action and place it on the node tree, under "ACTIONS" node
   creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
   */
  int
  Init(PHCompositeNode *);

  //! event processing
  /*!
   get all relevant nodes from top nodes (namely hit list)
   and pass that to the stepping action
   */
  int
  process_event(PHCompositeNode *);

  //! accessors (reimplemented)
  virtual PHG4Detector*
  GetDetector(void) const;
  virtual PHG4SteppingAction*
  GetSteppingAction(void) const;
  PHG4EventAction* GetEventAction() const {return eventAction_;}

  void
  SuperDetector(const std::string &name)
  {
    superdetector = name;
  }

  //! geometry manager PHG4Sector::Sector_Geometry
  PHG4Sector::Sector_Geometry &
  get_geometry()
  {
    return geom;
  }

  //! geometry manager PHG4Sector::Sector_Geometry
  void
  set_geometry(const PHG4Sector::Sector_Geometry &g)
  {
    geom = g;
  }

private:

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4SectorDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SectorSteppingAction* steppingAction_;

  PHG4EventAction *eventAction_;

  std::string superdetector;

  PHG4Sector::Sector_Geometry geom;
};

#endif
