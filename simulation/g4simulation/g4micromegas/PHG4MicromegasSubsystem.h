// Tell emacs that this is a C++ source
//  -*- C++ -*-.

#ifndef G4MICROMEGAS_PHG4MICROMEGASSUBSYSTEM_H
#define G4MICROMEGAS_PHG4MICROMEGASSUBSYSTEM_H

/*!
 * \file PHG4MicromegasSubsystem.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <string>                               // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4MicromegasDetector;
class PHG4DisplayAction;
class PHG4MicromegasSteppingAction;
class PHG4SteppingAction;

/*!
   * \brief Detector Subsystem module
   * The detector is constructed and registered via PHG4MicromegasDetector
   * \see PHG4MicromegasDetector
   * \see PHG4MicromegasSubsystem
   */
class PHG4MicromegasSubsystem : public PHG4DetectorSubsystem
{

  public:
  //! constructor
  PHG4MicromegasSubsystem(const std::string& name = "MICROMEGAS", int layer = 0);

  ~PHG4MicromegasSubsystem() override;

  /*!
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode*) override;

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode*) override;

  //!@name accessors (reimplemented)
  //@{
  PHG4Detector* GetDetector() const override;
  PHG4SteppingAction* GetSteppingAction() const override;
  //@}

  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const override;

  //! get the display action if display is started
  PHG4DisplayAction* GetDisplayAction() const override { return m_DisplayAction; }

  private:
  
  // \brief Set default parameter values
  void SetDefaultParameters() override;

  private:

  //! detector construction
  /*! derives from PHG4Detector */
  PHG4MicromegasDetector *m_Detector = nullptr;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4MicromegasSteppingAction *m_SteppingAction = nullptr;

  //! display attribute setting
  PHG4DisplayAction* m_DisplayAction = nullptr;
};

#endif // G4MICROMEGAS_PHG4MICROMEGASSUBSYSTEM_H
