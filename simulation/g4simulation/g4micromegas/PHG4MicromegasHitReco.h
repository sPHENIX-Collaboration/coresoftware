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

#include <gsl/gsl_rng.h>
#include <map>
#include <memory>
#include <string>
#include <utility>  // for pair

class CylinderGeomMicromegas;
class PHCompositeNode;
class PHG4Hit;
class TVector3;

class PHG4MicromegasHitReco : public SubsysReco, public PHParameterInterface
{

  public:
  explicit PHG4MicromegasHitReco(
    const std::string &name = "PHG4MicromegasHitReco",
    const std::string &detector = "MICROMEGAS");

  //! run initialization
  int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode*) override;

  //! parameters
  void SetDefaultParameters() override;

  //! set micromegas tiles
  void set_tiles( const MicromegasTile::List& tiles )
  { m_tiles = tiles; }

  private:

  //! setup tiles definition in CylinderGeom
  void setup_tiles(PHCompositeNode*);

  //! stores strip number and corresponding charge fraciton
  using charge_pair_t = std::pair<int, double>;

  //! list of charge fractions
  using charge_list_t = std::vector<charge_pair_t>;

  //! tile and list of charge fractions
  using charge_info_t = std::pair<int, charge_list_t>;

  //! get total number of electrons collected for a give g4hit
  /*! this accounts for the number of primary electrons, the detector gain, and fluctuations */
  uint get_electrons( PHG4Hit* ) const;
  
  //! distribute a Gaussian charge across adjacent strips
  charge_info_t distribute_charge( CylinderGeomMicromegas*, const TVector3& position, double sigma ) const;

  //! detector name
  std::string m_detector;

  //! timing window (ns)
  double m_tmin = 0;

  //! timing window (ns)
  double m_tmax = 0;

  //! number of primary electrons per GeV
  double m_electrons_per_gev = 0;
  
  //! min gain
  double m_gain = 0;
  
  //! electron cloud sigma (cm) after avalanche
  double m_cloud_sigma = 0.04;

  //! use zig zag pads
  bool m_zigzag_strips = true;
  
  //! micromegas tiles
  MicromegasTile::List m_tiles;

  //! rng de-allocator
  class Deleter
  {
    public:
    //! deletion operator
    void operator() (gsl_rng* rng) const { gsl_rng_free(rng); }
  };

  //! random generator that conform with sPHENIX standard
  /*! using a unique_ptr with custom Deleter ensures that the structure is properly freed when parent object is destroyed */
  std::unique_ptr<gsl_rng, Deleter> m_rng;

};

#endif
