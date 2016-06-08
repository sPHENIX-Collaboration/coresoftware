#ifndef PHG4HCALCELLRECO_H
#define PHG4HCALCELLRECO_H

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <string>
#include <map>
#include <vector>

class PHCompositeNode;
class PHG4CylinderCell;

class PHG4HcalCellReco : public SubsysReco
{
 public:

  PHG4HcalCellReco(const std::string &name = "HCALCELLRECO");

  virtual ~PHG4HcalCellReco(){}
  
  //! module initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode);
  
  void Detector(const std::string &d) {detector = d;}
  void etasize_nslat(const int i, const double deltaeta, const int nslat);
  void checkenergy(const int i=1) {chkenergyconservation = i;}
  void set_etabins(const int nbins=24) {netabins = nbins;}

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
  void set_size(const int i, const double sizeA, const int sizeB, const int what);
  int CheckEnergy(PHCompositeNode *topNode);
  double get_phi_slat_zero_low(const double radius, const double thickness, const double tiltangle);
  double get_phi_slat_zero_up(const double radius, const double thickness, const double tiltangle);
  static std::pair<double, double> get_etaphi(const double x, const double y, const double z);
  static double get_eta(const double radius, const double z);
  std::map<int, int>  binning;
  std::map<int, std::pair <double,int> > cell_size; // cell size in eta/nslats
  std::map<int, std::pair <double,double> > zmin_max; // zmin/zmax for each layer for faster lookup
  std::map<int, double> etastep;
  std::string detector;
  std::string hitnodename;
  std::string cellnodename;
  std::string geonodename;
  std::string seggeonodename;
  PHTimeServer::timer _timer;
  int nbins[2];
  int nslatscombined;
  int netabins;
  int chkenergyconservation;
  std::map<unsigned int, PHG4CylinderCell *> celllist;

  //! timing window size in ns. This is for a simple simulation of the ADC integration window starting from 0ns to this value. Default to infinity, i.e. include all hits
  double timing_min;
  double timing_max;

};

#endif
