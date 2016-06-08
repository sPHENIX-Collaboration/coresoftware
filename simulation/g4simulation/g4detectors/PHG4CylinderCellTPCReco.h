#ifndef PHG4CYLINDERCELLTPCRECO_H
#define PHG4CYLINDERCELLTPCRECO_H

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <string>
#include <map>
#include "TRandom3.h"

class PHCompositeNode;
class PHG4CylinderCell;

class PHG4CylinderCellTPCReco : public SubsysReco
{
public:
  
  PHG4CylinderCellTPCReco(const std::string &name = "CYLINDERTPCRECO");
  
  virtual ~PHG4CylinderCellTPCReco(){}
  
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
//   void etaphisize(const int i, const double deltaeta, const double deltaphi);
  void checkenergy(const int i=1) {chkenergyconservation = i;}
  void OutputDetector(const std::string &d) {outdetector = d;}

  void setDiffusion( double diff ){diffusion = diff;}
  void setElectronsPerKeV( double epk ){elec_per_kev = epk;}

  double get_timing_window_min() {return timing_min;}
  double get_timing_window_max() {return timing_max;}
  void set_timing_window(const double tmin, const double tmax) {
    timing_min = tmin; timing_max = tmax;
  }
  
  //! get timing window size in ns.
  double get_timing_window_size() const {return timing_max - timing_min;}
  //! set timing window size in ns. This is for a simple simulation of the ADC integration window starting from 0ns to this value. Default to infinity, i.e. include all hits
  void set_timing_window_size(const double s) {set_timing_window(0.0,s);}
  
protected:
//   void set_size(const int i, const double sizeA, const double sizeB, const int what);
//   int CheckEnergy(PHCompositeNode *topNode);
//   static std::pair<double, double> get_etaphi(const double x, const double y, const double z);
//   static double get_eta(const double radius, const double z);
//   bool lines_intersect( double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy, double* rx, double* ry);
//   bool line_and_rectangle_intersect( double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy, double* rr);
  
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
  
  int nbins[2];
  int chkenergyconservation;
  
  TRandom3 rand;

  double diffusion;
  double elec_per_kev;

  double timing_min;
  double timing_max;
};

#endif
