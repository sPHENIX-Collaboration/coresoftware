#ifndef G4MVTX_PHG4MVTXHITRECO_H
#define G4MVTX_PHG4MVTXHITRECO_H

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <map>
#include <string>
#include <memory>  // for unique_ptr

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

  std::string m_detector;

  double m_tmin;
  double m_tmax;
  double crossing_period = 106.0;

  class Deleter
  {
    public:
      //! delection operation
      void operator() (gsl_rng* rng) const { gsl_rng_free(rng); }
  };

  std::unique_ptr<gsl_rng, Deleter> m_rng;
};

#endif
