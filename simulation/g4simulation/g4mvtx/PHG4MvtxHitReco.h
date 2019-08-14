#ifndef G4MVTX_PHG4MVTXHITRECO_H
#define G4MVTX_PHG4MVTXHITRECO_H

#include <phparameter/PHParameterContainerInterface.h>

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>
#include <utility>                                      // for pair

class PHCompositeNode;

class PHG4MvtxHitReco : public SubsysReco, public PHParameterContainerInterface
{
 public:
  explicit PHG4MvtxHitReco(const std::string &name = "PHG4MvtxRECO");

  virtual ~PHG4MvtxHitReco() {}

  //! module initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  void Detector(const std::string &d) { detector = d; }

  double get_timing_window_min(const int i) { return tmin_max[i].first; }
  double get_timing_window_max(const int i) { return tmin_max[i].second; }
  void set_timing_window(const int detid, const double tmin, const double tmax);

  void SetDefaultParameters();

 protected:

  std::string detector;
  std::string hitnodename;
  std::string geonodename;
  std::map<int, std::pair<double, double> > tmin_max;
};

#endif
