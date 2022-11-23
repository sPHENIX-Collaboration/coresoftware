// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FILLSPACECHARGEMAPS_FILLSPACECHARGEMAPS_H
#define FILLSPACECHARGEMAPS_FILLSPACECHARGEMAPS_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>
#include <vector>

#include <cmath>      // for sin, asin, cos, floor, M_PI

// Forward declerations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TH1;
class TH2;
class TH3;
class TTree;

class fillSpaceChargeMaps : public SubsysReco
{
 public:
  fillSpaceChargeMaps(const std::string &name = "fillSpaceChargeMaps", const std::string &filename = "Hist.root");

  virtual ~fillSpaceChargeMaps();

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode * /*topNode*/) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode * /*topNode*/) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode * /*topNode*/) override;

  void SetFrequency(int freq);
  void SetBeamXing(const std::vector<int> &beamXs);
  void SetEvtStart(int newEvtStart);
  void SetUseIBFMap(bool useIBFMap = true);
  void SetGain(float ampGain = 2e3);
  void SetIBF(float ampIBFfrac = 0.004);
  void SetCollSyst(int coll_syst = 0);
  void SetAvg(int fAvg = 0);
  void UseSliming(int fSliming = 0);
  void UseFieldMaps(int shiftElectrons = 0);

 private:
  std::vector<double> getNewWeights(TH3 *_h_SC_ibf, TH2 *_h_modules_anode, TH2 *_h_modules_measuredibf, double _hit_r, double _hit_phi, double dr_bin, double dphi_bin, bool _fUseIBFMap);
  bool IsOverFrame(double r, double phi);
  std::vector<double> putOnPlane(double r, double phi);

  Fun4AllHistoManager *hm = nullptr;
  std::string _filename;
  std::set<std::string> _node_postfix;
  std::map<int, int> _timestamps;
  std::vector<int> _keys;
  TFile *outfile = nullptr;
  float _ampGain = 2e3;
  float _ampIBFfrac = 0.02;
  int _collSyst = 0;
  int _shiftElectrons = 0;

  double _freqKhz = 22;
  //int _beamxing = 0;
  std::vector<int> _beamxing;
  //std::vector<int> _beamxing_end;

  int _evtstart = 0;
  int _fAvg = 0;
  int _fSliming = 0;
  TTree *_rawHits = nullptr;
  int _isOnPlane = 0;
  float _hit_z = 0;
  float _hit_r = 0;
  float _hit_phi = 0;
  float _hit_eion = 0;
  float _ibf_vol = 0;
  float _amp_ele_vol = 0;
  float _event_timestamp = 0;
  float _event_bunchXing = 0;

  bool _fUseIBFMap = false;
  TH2 *_h_modules_anode = nullptr;
  TH2 *_h_modules_measuredibf = nullptr;
  TH1 *_h_hits = nullptr;
  TH1 *_h_R = nullptr;
  TH2 *_h_DC_E = nullptr;
  static const int nFrames = 30;
  TH3 *_h_SC_prim[nFrames] = {nullptr};
  TH3 *_h_SC_ibf[nFrames] = {nullptr};

  float f = 0.5;                    //for now, just pick the middle of the hit.  Do better later.
  float ns = 1e-9, s = 1.0;           // us=1e-6,ms=1e-3,
  float mm = 1.0, cm = 10.0, m = 1000.0;  //um=1e-3,//changed to make 'cm' 1.0, for convenience.
  float kHz = 1e3, MHz = 1e6;       //Hz=1,
  float V = 1.0;
  //used two ways:  1) to apply units to variables when defined
  //                2) to divide by certain units so that those variables are expressed in those units.

  //float ionMobility=3.37*cm*cm/V/s;
  float ionMobility = 1.65 * cm * cm / V / s;
  float vIon = ionMobility * 400 * V / cm;
  //float vIon=16.0*um/us;
  float mbRate = _freqKhz * kHz;
  float xingRate = 9.383 * MHz;
  //float mean = mbRate/xingRate;
  //float z_rdo=105.5*cm;
  //float rmin=20*cm;
  //float rmax=78*cm;

  double Ne_dEdx = 1.56 / cm;    // keV/cm
  double CF4_dEdx = 7.00 / cm;   // keV/cm
  double Ne_NTotal = 43 / cm;    // Number/cm
  double CF4_NTotal = 100 / cm;  // Number/cm
  //double Tpc_NTot = 0.90 * Ne_NTotal + 0.10 * CF4_NTotal;
  //double Tpc_dEdx = 0.90 * Ne_dEdx + 0.10 * CF4_dEdx;
  double Tpc_NTot = 0.50 * Ne_NTotal + 0.50 * CF4_NTotal;
  double Tpc_dEdx = 0.50 * Ne_dEdx + 0.50 * CF4_dEdx;

  //double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
  double Tpc_ElectronsPerGeV = Tpc_NTot / Tpc_dEdx * 1e6;  //electrons per gev.
  double phi_dead_bins[24] ={ 6.5314-2 * M_PI, 6.545-2 * M_PI, 
                                  0.7718, 0.7854, 
                                  1.2954, 1.309, 
                                  1.819, 1.8326, 
                                  2.3426, 2.3562, 
                                  2.8662, 2.8798, 
                                  3.3898, 3.4034, 
                                  3.9134, 3.927, 
                                  4.437, 4.4506, 
                                  4.9606, 4.9742, 
                                  5.4842, 5.4978, 
                                  6.0078, 6.0214};
  //int nr=159;
  //int nphi=360;
  //int nz=62*2;

  //double hrstep=(rmax-rmin)/nr;
  //double hphistep=2*pi/nphi;
  //double hzstep=z_rdo/nz;

  //int nBeams = z_rdo/(vIon/xingRate); //numaber of beamcrossings to fill TPC

  float _mbRate = 0;
  float _xingRate = 0;
  float _mean = 0;
};

#endif  // FILLSPACECHARGEMAPS_FILLSPACECHARGEMAPS_H
