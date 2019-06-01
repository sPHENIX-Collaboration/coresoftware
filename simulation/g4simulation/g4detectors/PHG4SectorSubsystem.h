// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SECTORSUBSYSTEM_H
#define G4DETECTORS_PHG4SECTORSUBSYSTEM_H

#include "PHG4SectorConstructor.h"

#include <g4main/PHG4Subsystem.h>

#include <string>                   // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4SectorDetector;
class PHG4SteppingAction;

class PHG4SectorSubsystem : public PHG4Subsystem
{
 public:
  //! constructor
  PHG4SectorSubsystem(const std::string& name = "Sector");

  //! destructor
  virtual ~PHG4SectorSubsystem();

  //! init
  /*!
   creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
   reates the stepping action and place it on the node tree, under "ACTIONS" node
   creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
   */
  int Init(PHCompositeNode*);

  //! event processing
  /*!
   get all relevant nodes from top nodes (namely hit list)
   and pass that to the stepping action
   */
  int process_event(PHCompositeNode*);

  //! accessors (reimplemented)
  virtual PHG4Detector*
  GetDetector(void) const;
  virtual PHG4SteppingAction* GetSteppingAction(void) const { return m_SteppingAction; }

  PHG4DisplayAction* GetDisplayAction() const { return m_DisplayAction; }

  void
  SuperDetector(const std::string& name)
  {
    superdetector = name;
  }

  //! geometry manager PHG4Sector::Sector_Geometry
  PHG4Sector::Sector_Geometry&
  get_geometry()
  {
    return geom;
  }

  //! geometry manager PHG4Sector::Sector_Geometry
  void
  set_geometry(const PHG4Sector::Sector_Geometry& g)
  {
    geom = g;
  }

 private:
  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4SectorDetector* m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;

  std::string superdetector;

  PHG4Sector::Sector_Geometry geom;
};

#endif
