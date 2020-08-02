// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERSUBSYSTEM_H
#define G4DETECTORS_PHG4CYLINDERSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

#include <array>  // for array
#include <string>  // for string

class PHCompositeNode;
class PHG4CylinderDetector;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4SteppingAction;

class PHG4CylinderSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4CylinderSubsystem(const std::string& name = "CYLINDER", const int layer = 0);

  //! destructor
  virtual ~PHG4CylinderSubsystem(void);

  //! init runwise stuff
  /*!
  creates the m_Detector object
  creates the stepping action
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
  PHG4SteppingAction* GetSteppingAction(void) const { return m_SteppingAction; }

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
  virtual bool CanBeMotherSubsystem() const { return true; }

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4CylinderDetector* m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;

  //! Color setting if we want to override the default
  std::array<double, 4> m_ColorArray;
};

#endif  // G4DETECTORS_PHG4CYLINDERSUBSYSTEM_H
