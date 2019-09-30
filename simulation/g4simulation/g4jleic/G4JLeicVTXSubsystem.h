// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4JLEIC_G4JLEICVTXSUBSYSTEM_H
#define G4JLEIC_G4JLEICVTXSUBSYSTEM_H

#include <g4detectors/PHG4DetectorGroupSubsystem.h>

#include <string>

class PHCompositeNode;
class PHG4Detector;
class G4JLeicVTXDetector;
class PHG4SteppingAction;

/**
   * \brief Fun4All module to simulate the Barrel PSTOF detector.
   *
   * The detector is constructed and registered via G4JLeicVTXDetector
   *
   * The PHG4SteppingAction needs to be written, but will provide the info for the hit time
   *
   * \see G4JLeicVTXDetector
   * \see G4JLeicVTXSubsystem
   *
   */
class G4JLeicVTXSubsystem : public PHG4DetectorGroupSubsystem
{
 public:
  //! constructor
  G4JLeicVTXSubsystem(const std::string& name = "PSTOF");

  //! destructor
  virtual ~G4JLeicVTXSubsystem(void)
  {
  }

  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  virtual int InitRunSubsystem(PHCompositeNode*);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  virtual int process_event(PHCompositeNode*);

  //! accessors (reimplemented)
  virtual PHG4Detector* GetDetector(void) const;
  virtual PHG4SteppingAction* GetSteppingAction(void) const;
  //! Print info (from SubsysReco)
  virtual void Print(const std::string& what = "ALL") const;

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  G4JLeicVTXDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* steppingAction_;
};

#endif  // G4JLEIC_G4JLEICVTXSUBSYSTEM_H
