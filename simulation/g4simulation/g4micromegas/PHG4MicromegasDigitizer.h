// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4MICROMEGAS_PHG4MICROMEGASDIGITIZER_H
#define G4MICROMEGAS_PHG4MICROMEGASDIGITIZER_H

/*!
 * \file PHG4MicromegasDigitizer.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/SubsysReco.h>

#include <phparameter/PHParameterInterface.h>

#include <gsl/gsl_rng.h>

#include <memory>
#include <string>                              // for string

class PHCompositeNode;

class PHG4MicromegasDigitizer : public SubsysReco, public PHParameterInterface
{

  public:
  PHG4MicromegasDigitizer(const std::string &name = "PHG4MicromegasDigitizer");

  //! run initialization
  int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  //! parameters
  void SetDefaultParameters() override;

  private:

  //! add noise to a measurement
  double add_noise() const;

  //! threshold (electrons)
  double m_adc_threshold = 2700;

  //! noise (electrons)
  double m_enc = 670;
  
  //! pedestal (electrons)
  double m_pedestal = 50000;
  
  //! conversion factor mv/fc
  double m_volts_per_charge = 20;

  //! conversion factor (mv/electron)
  double m_volt_per_electron_signal = 0;

  //! conversion factor (mv/electron)
  double m_volt_per_electron_noise = 0;

  //! conversion factor (adc/mv)
  /*! this is a fixed parameter, from SAMPA */
  static constexpr double m_adc_per_volt = 1024./2200;
  
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
