#ifndef PHG4SILICONTRACKERCELLRECO_H
#define PHG4SILICONTRACKERCELLRECO_H

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

#ifndef __CINT__
#include <gsl/gsl_vector.h>
#endif

#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class PHG4Cell;

class PHG4SiliconTrackerCellReco : public SubsysReco
{
 public:
  PHG4SiliconTrackerCellReco(const std::string &name = "SILICON_TRACKER");

  virtual ~PHG4SiliconTrackerCellReco();
  //! module initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  //! end of process
  int End(PHCompositeNode *topNode);

  void Detector(const std::string &d) { detector = d; }
  void checkenergy(const int i = 1) { chkenergyconservation = i; }
  double get_timing_window_min(const int i) { return tmin_max[i].first; }
  double get_timing_window_max(const int i) { return tmin_max[i].second; }
  void set_timing_window(const int i, const double tmin, const double tmax)
  {
    tmin_max[i] = std::make_pair(tmin, tmax);
  }
  void set_timing_window_defaults(const double tmin, const double tmax)
  {
    tmin_default = tmin;
    tmax_default = tmax;
  }

  double circle_rectangle_intersection(double x1, double y1,  double x2,  double y2,  double mx, double my,  double r);
  double sA(double r, double x, double y) ;
  
 protected:
  int CheckEnergy(PHCompositeNode *topNode);
  std::map<int, int> binning;
  std::string detector;
  std::string hitnodename;
  std::string cellnodename;
  std::string geonodename;
  PHTimeServer::timer _timer;
  int nbins[2];
  int chkenergyconservation;
  std::map<std::string, PHG4Cell *> celllist;  // This map holds the hit cells

  double tmin_default;
  double tmax_default;
  std::map<int, std::pair<double, double> > tmin_max;
#ifndef __CINT__
  gsl_vector *m_LocalOutVec;
  gsl_vector *m_PathVec;
  gsl_vector *m_SegmentVec;
#endif
};

#endif
