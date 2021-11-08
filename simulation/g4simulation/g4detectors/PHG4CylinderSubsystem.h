// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERSUBSYSTEM_H
#define G4DETECTORS_PHG4CYLINDERSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

#include <array>   // for array
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
  ~PHG4CylinderSubsystem(void) override;

  //! init runwise stuff
  /*!
  creates the m_Detector object
  creates the stepping action
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
  PHG4SteppingAction* GetSteppingAction(void) const override { return m_SteppingAction; }

  PHG4DisplayAction* GetDisplayAction() const override { return m_DisplayAction; }
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
  bool CanBeMotherSubsystem() const override { return true; }

  // this is just needed for use as reference plane for projections
  // this is the only detector using this - there is no need to add
  // this to our parameters
  void SaveAllHits(bool i = true) { m_SaveAllHitsFlag = i; }

 private:
  void SetDefaultParameters() override;

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4CylinderDetector* m_Detector = nullptr;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction = nullptr;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction = nullptr;

  bool m_SaveAllHitsFlag = false;
  //! Color setting if we want to override the default
  std::array<double, 4> m_ColorArray;
};

#endif  // G4DETECTORS_PHG4CYLINDERSUBSYSTEM_H
