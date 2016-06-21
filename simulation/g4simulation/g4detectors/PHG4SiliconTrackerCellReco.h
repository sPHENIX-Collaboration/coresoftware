#ifndef PHG4SILICONTRACKERCELLRECO_H
#define PHG4SILICONTRACKERCELLRECO_H

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <string>
#include <map>
#include <vector>

class PHCompositeNode;
class PHG4CylinderCell;

class PHG4SiliconTrackerCellReco : public SubsysReco
{
 public:

  PHG4SiliconTrackerCellReco(const std::string &name);

  virtual ~PHG4SiliconTrackerCellReco(){}
  
  //! module initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode);
  
  void Detector(const std::string &d) {detector = d;}
  void checkenergy(const int i=1) {chkenergyconservation = i;}

  double get_timing_window_min(const int i) {return tmin_max[i].first;}
  double get_timing_window_max(const int i) {return tmin_max[i].second;}
  void   set_timing_window(const int i, const double tmin, const double tmax) {
    tmin_max[i] = std::make_pair(tmin,tmax);
  }
  void   set_timing_window_defaults(const double tmin, const double tmax) {
    tmin_default = tmin; tmax_default = tmax;
  }
  
 protected:
  //void set_size(const int i, const double sizeA, const int sizeB, const int what);
  int CheckEnergy(PHCompositeNode *topNode);
  //double get_phi_slat_zero_low(const double radius, const double thickness, const double tiltangle);
  //double get_phi_slat_zero_up(const double radius, const double thickness, const double tiltangle);
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
  int chkenergyconservation;
  int layer;
  //std::map<unsigned int, PHG4CylinderCell *> celllist;
  std::map<std::string, PHG4CylinderCell*> celllist;  // This map holds the hit cells

  double tmin_default;
  double tmax_default;
  std::map<int, std::pair<double,double> > tmin_max;
};

#endif
