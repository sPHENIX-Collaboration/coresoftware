#ifndef PHG4CYLINDERCELLRECO_H
#define PHG4CYLINDERCELLRECO_H

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <string>
#include <map>

class PHCompositeNode;
class PHG4CylinderCell;

class PHG4CylinderCellReco : public SubsysReco
{
 public:

  PHG4CylinderCellReco(const std::string &name = "CYLINDERRECO");

  virtual ~PHG4CylinderCellReco(){}
  
  //! module initialization
  int InitRun(PHCompositeNode *topNode);
  
  //! run initialization
  int Init(PHCompositeNode *topNode) {return 0;}
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode);
  
  void Detector(const std::string &d);
  void cellsize(const int i, const double sr, const double sz);
  void etaphisize(const int i, const double deltaeta, const double deltaphi);
  void checkenergy(const int i=1) {chkenergyconservation = i;}
  void OutputDetector(const std::string &d) {outdetector = d;}

  double get_timing_window_min(const int i) {return tmin_max[i].first;}
  double get_timing_window_max(const int i) {return tmin_max[i].second;}
  void   set_timing_window(const int i, const double tmin, const double tmax) {
    tmin_max[i] = std::make_pair(tmin,tmax);
  }
  void   set_timing_window_defaults(const double tmin, const double tmax) {
    tmin_default = tmin; tmax_default = tmax;
  }

 protected:
  void set_size(const int i, const double sizeA, const double sizeB, const int what);
  int CheckEnergy(PHCompositeNode *topNode);
  static std::pair<double, double> get_etaphi(const double x, const double y, const double z);
  static double get_eta(const double radius, const double z);
  bool lines_intersect( double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy, double* rx, double* ry);
  bool line_and_rectangle_intersect( double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy, double* rr);

  std::map<int, int>  binning;
  std::map<int, std::pair <double,double> > cell_size; // cell size in phi/z
  std::map<int, std::pair <double,double> > zmin_max; // zmin/zmax for each layer for faster lookup
  std::map<int, double> phistep;
  std::map<int, double> etastep;
  std::string detector;
  std::string outdetector;
  std::string hitnodename;
  std::string cellnodename;
  std::string geonodename;
  std::string seggeonodename;
  std::map<int, std::pair<int, int> > n_phi_z_bins;
  std::map<std::string, PHG4CylinderCell*> cellptmap;  // This map holds the hit cells
  std::map<std::string, PHG4CylinderCell*>::iterator it;

  PHTimeServer::timer _timer;
  int nbins[2];
  int chkenergyconservation;

  //! timing window size in ns. This is for a simple simulation of the ADC integration window starting from 0ns to this value. Default to infinity, i.e. include all hits
  double tmin_default;
  double tmax_default;
  std::map<int, std::pair<double,double> > tmin_max;
};

#endif
