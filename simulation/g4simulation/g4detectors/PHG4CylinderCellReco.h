#ifndef G4DETECTORS_PHG4CYLINDERCELLRECO_H
#define G4DETECTORS_PHG4CYLINDERCELLRECO_H


#include <phparameter/PHParameterContainerInterface.h>

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

#include <map>
#include <set>
#include <string>

class PHCompositeNode;
class PHG4Cell;

class PHG4CylinderCellReco : public SubsysReco, public PHParameterContainerInterface
{
 public:

  explicit PHG4CylinderCellReco(const std::string &name = "CYLINDERRECO");

  virtual ~PHG4CylinderCellReco(){}
  
  //! module initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);

  int ResetEvent(PHCompositeNode *topNode);

  void SetDefaultParameters();

  void Detector(const std::string &d);
  void cellsize(const int i, const double sr, const double sz);
  void etaphisize(const int i, const double deltaeta, const double deltaphi);
  void checkenergy(const int i=1) {chkenergyconservation = i;}
  void OutputDetector(const std::string &d) {outdetector = d;}

  double get_timing_window_min(const int i) {return tmin_max[i].first;}
  double get_timing_window_max(const int i) {return tmin_max[i].second;}
  void   set_timing_window(const int detid, const double tmin, const double tmax);

 protected:
  void set_size(const int i, const double sizeA, const double sizeB);
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
  std::set<int> implemented_detid;
  std::string detector;
  std::string outdetector;
  std::string hitnodename;
  std::string cellnodename;
  std::string geonodename;
  std::string seggeonodename;
  std::map<int, std::pair<int, int> > n_phi_z_bins;
  std::map<unsigned long long, PHG4Cell*> cellptmap;  // This map holds the hit cells
  std::map<unsigned long long, PHG4Cell*>::iterator it;
  std::map<int, std::pair<double,double> > tmin_max;

  PHTimeServer::timer _timer;
  int nbins[2];
  int chkenergyconservation;

  double sum_energy_before_cuts;
  double sum_energy_g4hit;
};

#endif
