// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id: $

/*!
 * \file PHG4GDMLSubsystem.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4DETECTORS_PHG4GDMLSUBSYSTEM_H
#define G4DETECTORS_PHG4GDMLSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

#include <string>

class PHCompositeNode;
class PHG4GDMLDetector;
class PHG4Detector;
class PHG4SteppingAction;

/*!
 * \brief PHG4GDMLSubsystem is a generic detector built from a GDML import
 */
class PHG4GDMLSubsystem : public PHG4DetectorSubsystem
{
 public:
  explicit PHG4GDMLSubsystem(const std::string& name);
  ~PHG4GDMLSubsystem() override;

  /*!
  creates the m_Detector object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode*) override;

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode*) override;

  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const override;

  //! accessors (reimplemented)
  PHG4Detector* GetDetector(void) const override;
  PHG4SteppingAction* GetSteppingAction(void) const override { return nullptr; }

 private:
  void SetDefaultParameters() override;

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4GDMLDetector* m_Detector = nullptr;
};

#endif /* PHG4GDMLSUBSYSTEM_H_ */
