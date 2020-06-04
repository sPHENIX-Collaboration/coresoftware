#ifndef G4MVTX_PHG4MicromegasHitReco_H
#define G4MVTX_PHG4MicromegasHitReco_H

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
  int InitRun(PHCompositeNode*);

  //! event processing
  int process_event(PHCompositeNode*);

  //!@name modifiers
  //@{

  //! parameters
  void SetDefaultParameters();

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
