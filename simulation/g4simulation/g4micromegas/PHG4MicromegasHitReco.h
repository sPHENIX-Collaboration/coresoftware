// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4MICROMEGAS_PHG4MICROMEGASHITRECO_H
#define G4MICROMEGAS_PHG4MICROMEGASHITRECO_H

/*!
 * \file PHG4MicromegasHitReco.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <micromegas/MicromegasTile.h>
#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>
#include <utility>  // for pair

class CylinderGeomMicromegas;
class PHCompositeNode;
class TVector3;

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

  //! set micromegas tiles
  void set_tiles( const MicromegasTile::List& tiles )
  { m_tiles = tiles; }

  //@}

  private:

  //! setup tiles definition in CylinderGeom
  void setup_tiles(PHCompositeNode*);

  //! stores strip number and corresponding charge fraciton
  using charge_pair_t = std::pair<int, float>;
    
  //! list of charge fractions
  using charge_list_t = std::vector<charge_pair_t>;
  
  //! tile and list of charge fractions
  using charge_info_t = std::pair<int, charge_list_t>;

  //! distribute a charge across adjacent strips
  charge_info_t distribute_charge( CylinderGeomMicromegas*, const TVector3& position, double sigma ) const;
  
  //! detector name
  std::string m_detector;

  //! timing window (ns)
  double m_tmin = 0;
  double m_tmax = 0;

  //! micromegas tiles
  MicromegasTile::List m_tiles;

};

#endif
