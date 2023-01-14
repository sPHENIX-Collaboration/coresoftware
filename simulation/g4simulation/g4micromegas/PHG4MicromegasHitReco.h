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

#include <gsl/gsl_rng.h>
#include <memory>
#include <string>
#include <utility>
#include <vector>

class ActsGeometry;

class CylinderGeomMicromegas;
class PHCompositeNode;
class PHG4Hit;
class TVector2;

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

  private:

  //! return full geo node name, that also contains tile information
  std::string full_geonodename() const
  { return "CYLINDERGEOM_" + m_detector + "_FULL"; }

  //! get total number of electrons collected for a give g4hit
  /*! this accounts for the number of primary electrons, the detector gain, and fluctuations */
  uint get_primary_electrons( PHG4Hit* ) const;

  //! get single electron amplification
  uint get_single_electron_amplification() const;

  //! stores strip number and corresponding charge fraction
  using charge_pair_t = std::pair<int, double>;

  //! map strip number to charge fraction
  using charge_list_t = std::vector<charge_pair_t>;

  //! distribute a Gaussian charge across adjacent strips
  charge_list_t distribute_charge( CylinderGeomMicromegas*, uint tileid, const TVector2& local_position, double sigma ) const;

  //! acts geometry
  ActsGeometry* m_acts_geometry = nullptr;
  
  //! detector name
  std::string m_detector;

  //! timing window (ns)
  double m_tmin = -20;

  //! timing window (ns)
  double m_tmax = 800;

  //! number of primary electrons per GeV
  double m_electrons_per_gev = 0;

  //! min gain
  double m_gain = 0;

  //! electron cloud sigma (cm) after avalanche
  double m_cloud_sigma = 0.04;

  //! electron transverse diffusion (cm/sqrt(cm))
  double m_diffusion_trans = 0.03;

  //! additional smearing of primary electrons (cm)
  /** it is used to adjust the Micromegas resolution to actual measurements */
  double m_added_smear_sigma_z = 0;
  double m_added_smear_sigma_rphi = 0;

  
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
