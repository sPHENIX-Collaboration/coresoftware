#ifndef G4MVTX_PHG4MVTXHITRECO_H
#define G4MVTX_PHG4MVTXHITRECO_H

#include <phparameter/PHParameterInterface.h>
#include <trackbase/TrkrDefs.h>

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <map>
#include <memory>  // for unique_ptr
#include <string>

class PHCompositeNode;

class PHG4MvtxHitReco : public SubsysReco, public PHParameterInterface
{
 public:
  explicit PHG4MvtxHitReco(
      const std::string &name = "PHG4MvtxHitReco",
      const std::string &detector = "MVTX");

  ~PHG4MvtxHitReco() override {}

  //! module initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  void Detector(const std::string &d) { m_detector = d; }

  //! TODO keep it for backward compatibily. remove after PR merged
  //In the future use the relevant set parameter function
  void set_timing_window(const int detid, const double tmin, const double tmax);

  //! parameters
  void SetDefaultParameters() override;

 private:
  std::pair<double, double> generate_alpide_pulse(const double energy_deposited);

  double generate_strobe_zero_tm_start();

  int get_strobe_frame(double alpide_time, double strobe_zero_tm_start);

  TrkrDefs::hitsetkey zero_strobe_bits(TrkrDefs::hitsetkey hitsetkey);

  std::string m_detector;

  double m_tmin;
  double m_tmax;
  double m_strobe_width;
  double m_strobe_separation;
  //double crossing_period = 106.0;
  double m_extended_readout_time = 0.0;

  bool m_in_sphenix_srdo = false;

  class Deleter
  {
   public:
    //! delection operation
    void operator()(gsl_rng *rng) const { gsl_rng_free(rng); }
  };

  std::unique_ptr<gsl_rng, Deleter> m_rng;
};

#endif
