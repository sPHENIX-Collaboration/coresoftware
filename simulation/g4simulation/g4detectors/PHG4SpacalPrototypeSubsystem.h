// $$Id: PHG4SpacalPrototypeSubsystem.h,v 1.2 2014/08/12 03:49:12 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2014/08/12 03:49:12 $$
 */

#ifndef PHG4SpacalPrototypeSubsystem_h
#define PHG4SpacalPrototypeSubsystem_h

#include "g4main/PHG4Subsystem.h"
#include "PHG4Parameters.h"

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

class PHG4SpacalPrototypeDetector;
class PHG4SpacalPrototypeSteppingAction;
class PHG4EventAction;

class PHG4SpacalPrototypeSubsystem : public PHG4Subsystem
{

public:

  //! constructor
  PHG4SpacalPrototypeSubsystem(const std::string &name =
      "CEMC");

  //! destructor
  virtual
  ~PHG4SpacalPrototypeSubsystem(void)
  {
  }

  //! init
  /*!
   creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
   reates the stepping action and place it on the node tree, under "ACTIONS" node
   creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
   */
  int
  InitRun(PHCompositeNode *);

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
  PHG4EventAction*
  GetEventAction() const
  {
    return eventAction_;
  }

  void
  SetActive(const int i = 1)
  {
    active = i;
  }
  void
  SetAbsorberActive(const int i = 1)
  {
    absorberactive = i;
  }
  void
  SuperDetector(const std::string &name)
  {
    superdetector = name;
    Params.set_name(superdetector);
  }
  const std::string
  SuperDetector()
  {
    return superdetector;
  }

  void
  Print(const std::string &what = "ALL") const;

  //! load the default parameter to param
  void
  SetDefaultParameters(PHG4Parameters * param);

  //! Get the parameters for readonly
  const PHG4Parameters &
  GetParameters() const
  {
    return Params;
  }

  //! Get the parameters for update. Useful fields are listed in SetDefaultParameters();
  PHG4Parameters &
  GetParameters()
  {
    return Params;
  }

  //! Overwrite the parameter. Useful fields are listed in SetDefaultParameters();
  void
  SetParameters(const PHG4Parameters & geom)
  {
    Params = geom;
  }

  //! use database?
  void
  UseDB(const bool b = true)
  {
    useDB = b;
  }

private:

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4SpacalPrototypeDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SpacalPrototypeSteppingAction* steppingAction_;
  PHG4EventAction *eventAction_;

  int active;
  int absorberactive;
  std::string detector_type;
  std::string superdetector;

  int useDB;
  PHG4Parameters Params;
};

#endif
