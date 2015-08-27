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
  void etasize_nslat(const int i, const double deltaeta, const int nslat);
  void checkenergy(const int i=1) {chkenergyconservation = i;}

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
  int chkenergyconservation;
  std::map<unsigned int, PHG4CylinderCell *> celllist;
};

#endif
