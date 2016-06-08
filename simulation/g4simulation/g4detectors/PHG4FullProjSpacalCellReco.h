#ifndef PHG4FullProjSpacalCellReco_H
#define PHG4FullProjSpacalCellReco_H

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <string>
#include <map>
#include <vector>

class PHCompositeNode;
class PHG4CylinderCell;

class PHG4FullProjSpacalCellReco : public SubsysReco
{
 public:

  PHG4FullProjSpacalCellReco(const std::string &name = "HCALCELLRECO");

  virtual ~PHG4FullProjSpacalCellReco(){}
  
  //! module initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode);
  
  void Detector(const std::string &d) {detector = d;}

  void checkenergy(const int i=1) {chkenergyconservation = i;}

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

  int CheckEnergy(PHCompositeNode *topNode);


  std::string detector;
  std::string hitnodename;
  std::string cellnodename;
  std::string geonodename;
  std::string seggeonodename;

  PHTimeServer::timer _timer;
  int chkenergyconservation;
  std::map<unsigned int, PHG4CylinderCell *> celllist;

  //! timing window size in ns. This is for a simple simulation of the ADC integration window starting from 0ns to this value. Default to infinity, i.e. include all hits
  double timing_min;
  double timing_max;

};

#endif
