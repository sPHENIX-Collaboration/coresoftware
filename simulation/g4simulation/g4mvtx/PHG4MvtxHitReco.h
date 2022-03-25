#ifndef G4MVTX_PHG4MVTXHITRECO_H
#define G4MVTX_PHG4MVTXHITRECO_H

#include <phparameter/PHParameterContainerInterface.h>

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <map>
#include <string>
#include <utility>  // for pair

class PHCompositeNode;

class PHG4MvtxHitReco : public SubsysReco, public PHParameterContainerInterface
{
 public:
  explicit PHG4MvtxHitReco(const std::string &name = "PHG4MvtxRECO");

  ~PHG4MvtxHitReco();

  //! module initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  void Detector(const std::string &d) { detector = d; }

  double get_timing_window_min(const int i) { return tmin_max[i].first; }
  double get_timing_window_max(const int i) { return tmin_max[i].second; }
  std::pair<double, double> generate_alpide_pulse(const double energy_deposited);
  double generate_strobe();
  void set_timing_window(const int detid, const double tmin, const double tmax);
  void use_strobed_mode(bool use) {m_use_strobe = use;}
  void set_strobe_properties(double width, double separation) {m_strobe_width = width; m_strobe_separation = separation;}

  void SetDefaultParameters() override;

 protected:
  std::string detector;
  std::string hitnodename;
  std::string geonodename;
  std::map<int, std::pair<double, double> > tmin_max;
  bool m_use_strobe = false;
  double m_strobe_width = 5.;
  double m_strobe_separation = 0.;

private:
  gsl_rng *m_RandomGenerator = nullptr;
};

#endif
