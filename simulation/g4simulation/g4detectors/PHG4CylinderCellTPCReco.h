#ifndef PHG4CYLINDERCELLTPCRECO_H
#define PHG4CYLINDERCELLTPCRECO_H

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

#include <TRandom3.h>

#include <string>
#include <map>

class PHCompositeNode;
class PHG4TPCDistortion;
class TH1F;
class TProfile2D;
class TStopwatch;

class PHG4CylinderCellTPCReco : public SubsysReco
{
public:
  
  PHG4CylinderCellTPCReco( int n_pixel=2, const std::string &name = "CYLINDERTPCRECO");
  
  virtual ~PHG4CylinderCellTPCReco();
  
  //! module initialization
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  
  //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode);
  
  void Detector(const std::string &d);
  void cellsize(const int i, const double sr, const double sz);
//   void etaphisize(const int i, const double deltaeta, const double deltaphi);
  void OutputDetector(const std::string &d) {outdetector = d;}

  void setHalfLength( double hz ){fHalfLength = hz;}
  void setDiffusionL( double diff ){fDiffusionL = diff;}
  void setDiffusionT( double diff ){fDiffusionT = diff;}
  void setDiffusion( double diff ){setDiffusionL(diff); setDiffusionT(diff);} //deprecated
  void setElectronsPerKeV( double epk ){elec_per_kev = epk;}
  void set_drift_velocity( const double cm_per_ns) { driftv = cm_per_ns;}
  
  double get_timing_window_min(const int i) {return tmin_max[i].first;}
  double get_timing_window_max(const int i) {return tmin_max[i].second;}
  void   set_timing_window(const int i, const double tmin, const double tmax) {
    tmin_max[i] = std::make_pair(tmin,tmax);
  }
  void   set_timing_window_defaults(const double tmin, const double tmax) {
    tmin_default = tmin; tmax_default = tmax;
  }

  //! distortion to the primary ionization
  void setDistortion (PHG4TPCDistortion * d) {distortion = d;}

protected:
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
  PHTimeServer::timer _timer;
  int nbins[2];
  
  TRandom3 rand;

  double fHalfLength;
  double fDiffusionT;
  double fDiffusionL;
  double elec_per_kev;
  double driftv;

  int num_pixel_layers;

  double tmin_default;
  double tmax_default;
  std::map<int, std::pair<double,double> > tmin_max;
  
  //! distortion to the primary ionization if not NULL
  PHG4TPCDistortion * distortion;
  TH1F *fHElectrons;
  TProfile2D *fHWindowP;
  TProfile2D *fHWindowZ;
  TProfile2D *fHMeanEDepPerCell;
  TProfile2D *fHMeanElectronsPerCell;
  TProfile2D *fHErrorRPhi;
  TProfile2D *fHErrorZ;
  TStopwatch *fSW;
  TH1F *fHTime;
};

#endif
