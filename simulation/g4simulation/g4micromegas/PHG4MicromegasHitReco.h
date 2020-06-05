// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4MICROMEGAS_PHG4MICROMEGASHITRECO_H
#define G4MICROMEGAS_PHG4MICROMEGASHITRECO_H

/*!
 * \file PHG4MicromegasHitReco.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>
#include <utility>  // for pair

class PHCompositeNode;

class PHG4MicromegasHitReco : public SubsysReco, public PHParameterInterface
{

  public:
  explicit PHG4MicromegasHitReco(
    const std::string &name = "PHG4MicromegasHitReco",
    const std::string &detector = "MICROMEGAS");

  //! module initialization
  int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode*) override;

  //!@name modifiers
  //@{

  //! parameters
  void SetDefaultParameters() override;

  //@}

  private:

  //! setup tiles definition in CylinderGeom
  void setup_tiles(PHCompositeNode*);

  //! detector name
  std::string m_detector;

  //! timing window (ns)
  double m_tmin = 0;
  double m_tmax = 0;
};

#endif
