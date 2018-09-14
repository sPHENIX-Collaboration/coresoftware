#ifndef PHG4BLOCKCELLRECO_H
#define PHG4BLOCKCELLRECO_H

#include <phparameter/PHParameterContainerInterface.h>

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

#include <map>
#include <set>
#include <string>

class PHCompositeNode;
class PHG4BlockCell;

class PHG4BlockCellReco : public SubsysReco, public PHParameterContainerInterface
{
 public:

  PHG4BlockCellReco(const std::string &name = "BLOCKRECO");

  virtual ~PHG4BlockCellReco(){}
  
  //! module initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  int ResetEvent(PHCompositeNode *topNode);

  void SetDefaultParameters();

  void Detector(const std::string &d) {detector = d;}
  void etaxsize(const int i, const double deltaeta, const double deltax);
  void checkenergy(const int i=1) {chkenergyconservation = i;}
  
  void   set_timing_window(const int detid, const double tmin, const double tmax);
  
 protected:
  void set_size(const int i, const double sizeA, const double sizeB, const int what);
  int CheckEnergy(PHCompositeNode *topNode);
  static std::pair<double, double> get_etaphi(const double x, const double y, const double z);
  static double get_eta(const double radius, const double z);
  bool lines_intersect( double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy, double* rx, double* ry);
  bool line_and_rectangle_intersect( double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy, double* rr);

  double sum_energy_g4hit;
  std::map<int, int>  binning;
  std::map<int, std::pair <double,double> > cell_size; // cell size in x/z
  std::map<int, std::pair <double,double> > zmin_max; // zmin/zmax for each layer for faster lookup
  std::map<int, double> xstep;
  std::map<int, double> etastep;
  std::map<int, std::pair<double,double> > tmin_max;
  std::set<int> implemented_detid;
  std::string detector;
  std::string hitnodename;
  std::string cellnodename;
  std::string geonodename;
  std::string seggeonodename;
  std::map<int, std::pair<int, int> > n_x_z_bins;
  PHTimeServer::timer _timer;
  int chkenergyconservation;

};

#endif
