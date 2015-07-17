// $$Id: PHG4SpacalSubsystem.h,v 1.2 2014/08/12 03:49:12 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2014/08/12 03:49:12 $$
 */

#ifndef PHG4SpacalSubsystem_h
#define PHG4SpacalSubsystem_h

#include "g4main/PHG4Subsystem.h"
#include "PHG4CylinderGeom_Spacalv2.h"

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

class PHG4SpacalDetector;
class PHG4SpacalSteppingAction;
class PHG4EventAction;

class PHG4SpacalSubsystem : public PHG4Subsystem
{

public:

  //! constructor
  PHG4SpacalSubsystem(const std::string &name = "PHG4SpacalSubsystem",
      const int layer = 0);

  //! destructor
  virtual
  ~PHG4SpacalSubsystem(void)
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
  SetLengthViaRapidityCoverage(const G4bool bl)
  {
    lengthViaRapidityCoverage = bl;
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
  }
  const std::string
  SuperDetector()
  {
    return superdetector;
  }

  void
  Print(const std::string &what = "ALL") const;

  typedef PHG4CylinderGeom_Spacalv2 SpacalGeom_t;

  const SpacalGeom_t &
  get_geom() const
  {
    return _geom;
  }

SpacalGeom_t &
  get_geom()
  {
    return _geom;
  }

  void
  set_geom(const SpacalGeom_t & geom)
  {
    _geom = geom;
  }

private:
  SpacalGeom_t _geom;

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4SpacalDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SpacalSteppingAction* steppingAction_;
  PHG4EventAction *eventAction_;

  int active;
  int absorberactive;
  int layer;
  G4bool lengthViaRapidityCoverage;
  std::string detector_type;
  std::string superdetector;
};

#endif
