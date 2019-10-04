// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKSUBSYSTEM_H
#define G4DETECTORS_PHG4BLOCKSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

#if !defined(__CINT__) || defined(__CLING__)
#include <array>   // for array
#endif

#include <string>  // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4BlockDetector;
class PHG4DisplayAction;
class PHG4SteppingAction;

class PHG4BlockSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4BlockSubsystem(const std::string& name = "BLOCK", const int layer = 0);

  //! destructor
  virtual ~PHG4BlockSubsystem();

  //! InitRunSubsystem
  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
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

  //! accessors (reimplemented)
  PHG4Detector* GetDetector() const;

  PHG4SteppingAction* GetSteppingAction() const { return m_SteppingAction; }

  PHG4DisplayAction* GetDisplayAction() const { return m_DisplayAction; }

  void set_color(const double red, const double green, const double blue, const double alpha = 1.)
  {
    m_ColorArray[0] = red;
    m_ColorArray[1] = green;
    m_ColorArray[2] = blue;
    m_ColorArray[3] = alpha;
  }
// this method is used to check if it can be used as mothervolume
// Subsystems which can be mothervolume need to implement this 
// and return true
  virtual bool CanBeMotherSubsystem() const {return true;}

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defines from PHG4Detector */
  PHG4BlockDetector* m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;
  //! Color setting if we want to override the default
#if !defined(__CINT__) || defined(__CLING__)
  std::array<double, 4> m_ColorArray;
#else
  double m_ColorArray[4];
#endif
};

#endif
