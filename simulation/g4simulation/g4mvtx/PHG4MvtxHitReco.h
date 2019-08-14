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

  double circle_rectangle_intersection(double x1, double y1, double x2, double y2, double mx, double my, double r);
  double sA(double r, double x, double y);

  std::map<int, int> binning;
  std::map<int, std::pair<double, double> > zmin_max;  // zmin/zmax for each layer for faster lookup
  std::map<int, double> etastep;
  std::string detector;
  std::string hitnodename;
  std::string cellnodename;
  std::string geonodename;
  std::string seggeonodename;
  int nbins[2];
  int chkenergyconservation;
  std::map<int, std::pair<double, double> > tmin_max;
  //std::map<unsigned long long, PHG4Cell *> celllist;  // This map holds the hit cells
};

#endif
