#ifndef PHG4CYLINDERCELLTPCRECO_H
#define PHG4CYLINDERCELLTPCRECO_H

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

// rootcint barfs with this header so we need to hide it
#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif


#include <string>
#include <map>

class PHCompositeNode;
class PHG4TPCDistortion;
class TH1;
class TF1;
class TProfile2D;
class PHG4CylinderCellGeom;

class PHG4CylinderCellTPCReco : public SubsysReco
{
public:
  
  PHG4CylinderCellTPCReco( const int n_pixel=2, const std::string &name = "CYLINDERTPCRECO");
  
  virtual ~PHG4CylinderCellTPCReco();
  
  //! module initialization
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  
  //! event processing
  int process_event(PHCompositeNode *topNode);
  
  void Detector(const std::string &d);
  void cellsize(const int i, const double sr, const double sz);
  void setZigzags(const bool zzpads){zigzag_pads = zzpads;}
  void populate_rectangular_phibins(PHG4CylinderCellGeom *geo, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_zigzag_phibins(PHG4CylinderCellGeom *geo, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_zbins(PHG4CylinderCellGeom *geo, const double z,  const double cloud_sig_zz[2], std::vector<int> &pad_zbin, std::vector<double> &pad_zbin_share);
  void OutputDetector(const std::string &d) {outdetector = d;}

  void setHalfLength(const double hz){fHalfLength = hz;}
  void setDiffusionL(const double diff){fDiffusionL = diff;}
  void setDiffusionT(const double diff){fDiffusionT = diff;}
  void setSigmaT    (const double gem) {sigmaT = gem;}  //  avalanche-induced charge spread...
  void setElectronsPerKeV(const double epk){elec_per_gev = epk*1e6;}
  void set_drift_velocity(const double cm_per_ns) { driftv = cm_per_ns;}
 void setEffElectronsThreshold(const double elec_threshold) { neffelectrons_threshold = elec_threshold;}

  void setSmearRPhi(const double v) {fFractRPsm=v;}
  void setSmearZ(const double v) {fFractZZsm=v;}

  void setShapingRMSLead(const double v) {fShapingLead=v;}
  void setShapingRMSTail(const double v) {fShapingTail=v;}
  
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
  std::map<int, int> binning;
  std::map<int, std::pair<double,double>> cell_size; // cell size in phi/z
  std::map<int, double> phistep;
  std::map<int, double> etastep;
  std::vector<int> adc_zbin;
  std::vector<int> pad_phibin;
  std::vector<double> pad_phibin_share;
  std::vector<double> adc_zbin_share;
  std::string detector;
  std::string outdetector;
  std::string hitnodename;
  std::string cellnodename;
  std::string geonodename;
  std::string seggeonodename;
  std::map<int, std::pair<int, int>> n_phi_z_bins;
  PHTimeServer::timer _timer;
  int nbins[2];
  
  double fHalfLength;
  double fDiffusionT;
  double fDiffusionL;
  double sigmaT;
  double elec_per_gev;
  double driftv;
  double neffelectrons_threshold;

  TF1 *fpad[10];
  TF1 *fcharge;

  int num_pixel_layers;

  double tmin_default;
  double tmax_default;
  std::map<int,std::pair<double,double>> tmin_max;
  
  //! distortion to the primary ionization if not NULL
  PHG4TPCDistortion * distortion;
  TH1 *fHElectrons;
  TProfile2D *fHWindowP;
  TProfile2D *fHWindowZ;
  TProfile2D *fHMeanEDepPerCell;
  TProfile2D *fHMeanElectronsPerCell;
  TProfile2D *fHErrorRPhi;
  TProfile2D *fHErrorZ;
  double fFractRPsm;
  double fFractZZsm;
  double fShapingLead;
  double fShapingTail;
  bool zigzag_pads;

#ifndef __CINT__
  //! random generator that conform with sPHENIX standard
  gsl_rng *RandomGenerator;
#endif

};

#endif
