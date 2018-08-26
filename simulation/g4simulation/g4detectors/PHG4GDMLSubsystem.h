// $Id: $

/*!
 * \file PHG4GDMLSubsystem.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHG4GDMLSUBSYSTEM_H_
#define PHG4GDMLSUBSYSTEM_H_

#include "PHG4DetectorSubsystem.h"
#include <string>

class PHG4GDMLDetector;

/*!
 * \brief PHG4GDMLSubsystem is a generic detector built from a GDML import
 */
class PHG4GDMLSubsystem : public PHG4Subsystem
{
 public:
  PHG4GDMLSubsystem(const std::string &name);
  virtual ~PHG4GDMLSubsystem();


  /*!
  creates the m_Detector object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode*);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode*);

  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const;

  //! accessors (reimplemented)
  PHG4Detector* GetDetector(void) const;
  PHG4SteppingAction* GetSteppingAction(void) const { return nullptr; }

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4GDMLDetector* m_Detector;

};

#endif /* PHG4GDMLSUBSYSTEM_H_ */
