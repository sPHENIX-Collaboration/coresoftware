#ifndef PHG4Prototype2HcalCELLRECO_H
#define PHG4Prototype2HcalCELLRECO_H

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <string>
#include <map>
#include <vector>

class PHCompositeNode;
class PHG4CylinderCell;

class PHG4Prototype2HcalCellReco : public SubsysReco
{
 public:

  PHG4Prototype2HcalCellReco(const std::string &name = "Prototype2HcalCELLRECO");

  virtual ~PHG4Prototype2HcalCellReco(){}
  
  //! module initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode);
  
  void Detector(const std::string &d) {detector = d;}
  void checkenergy(const int i=1) {chkenergyconservation = i;}

  //! get timing window size in ns.
  double get_timing_window_size() const {return timing_window_size;}
  //! set timing window size in ns. This is for a simple simulation of the ADC integration window starting from 0ns to this value. Default to infinity, i.e. include all hits
  void set_timing_window_size(const double s) {timing_window_size = s;}

 protected:
  int CheckEnergy(PHCompositeNode *topNode);
  std::map<int, int>  binning;
  std::map<int, std::pair <double,int> > cell_size; // cell size in eta/nslats
  std::map<int, std::pair <double,double> > zmin_max; // zmin/zmax for each layer for faster lookup
  std::map<int, double> etastep;
  std::string detector;
  std::string hitnodename;
  std::string cellnodename;
  std::string geonodename;
  std::string seggeonodename;
  std::map<int, std::pair<int, int> > n_phi_z_bins;
  PHTimeServer::timer _timer;
  int nbins[2];
  int nslatscombined;
  int chkenergyconservation;

  //! timing window size in ns. This is for a simple simulation of the ADC integration window starting from 0ns to this value. Default to infinity, i.e. include all hits
  double timing_window_size;

};

#endif
