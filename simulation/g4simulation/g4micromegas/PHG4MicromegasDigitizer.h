// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4MICROMEGAS_PHG4MICROMEGASDIGITIZER_H
#define G4MICROMEGAS_PHG4MICROMEGASDIGITIZER_H

/*!
 * \file PHG4MicromegasDigitizer.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>
#include <memory>

class PHCompositeNode;

class PHG4MicromegasDigitizer : public SubsysReco
{

  public:
  PHG4MicromegasDigitizer(const std::string &name = "PHG4MicromegasDigitizer");

  //! run initialization
  int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  //! adc parameters
  void set_adc_scale(unsigned int max_adc, double energy_scale)
  {
    m_max_adc  = max_adc;
    m_energy_scale = energy_scale;
  }

  //! threshold
  void set_energy_threshold(double value)
  { m_energy_threshold = value; }

  private:

  // settings
  unsigned int m_max_adc = 0;
  double m_energy_scale = 1;
  double m_energy_threshold = 0;

  //! rng de-allocator
  class Deleter
  {
    public:
    //! deletion operator
    void operator() (gsl_rng* rng) const
    { gsl_rng_free(rng); }
  };

  //! random generator that conform with sPHENIX standard
  /*! using a unique_ptr with custom Deleter ensures that the structure is properly freed when parent object is destroyed */
  std::unique_ptr<gsl_rng, Deleter> m_rng;

};

#endif
