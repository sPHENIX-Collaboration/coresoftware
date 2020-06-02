// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MICROMEGASSUBSYSTEM_H
#define MICROMEGASSUBSYSTEM_H

/*!
 * \file MicromegasSubsystem.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <g4detectors/PHG4DetectorSubsystem.h>

class PHCompositeNode;
class PHG4Detector;
class MicromegasDetector;
class MicromegasSteppingAction;
class PHG4SteppingAction;

/**
   * \brief Detector Subsystem module
   * The detector is constructed and registered via MicromegasDetector
   * \see MicromegasDetector
   * \see MicromegasSubsystem
   */
class MicromegasSubsystem : public PHG4DetectorSubsystem
{
  
  public:
  //! constructor
  MicromegasSubsystem(const std::string& name = "Micromegas");

  //! destructor
  virtual ~MicromegasSubsystem() = default;

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

  protected:
  // \brief Set default parameter values
  void SetDefaultParameters() override;

  private:
  //! detector construction
  /*! derives from PHG4Detector */
  MicromegasDetector *m_Detector = nullptr;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  MicromegasSteppingAction *m_SteppingAction = nullptr;
};

#endif // MICROMEGASSUBSYSTEM_H
