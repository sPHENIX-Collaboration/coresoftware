// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FILLSPACECHARGEMAPS_H
#define FILLSPACECHARGEMAPS_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>
#include <vector>

#include <math.h>

// Forward declerations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TH1;
class TH2;
class TH3;
class TNtuple;
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
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  //int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  void SetFrequency(int freq);
  //void SetBeamXing(int newBeamXing);
  void SetBeamXing(std::vector<int> beamXs);
  //void SetBeamXingEnd(std::vector<int> beamXs_end);
  void SetEvtStart(int newEvtStart);
  void SetUseIBFMap(bool useIBFMap = true);
  void SetGain(float ampGain=2e3);
  void SetIBF(float ampIBFfrac=0.004);
  void SetCollSyst(int coll_syst=0);
  void SetAvg(int fAvg=0);
  void UseSliming(int fSliming=0);
  void UseFieldMaps(int shiftElectrons=0);

  double pi = 2 * acos(0.0);

 protected:
   Fun4AllHistoManager *hm;
   std::string _filename;
   std::set<std::string> _node_postfix;
   std::map<int,int> _timestamps;
   std::vector<int> _keys;
   TFile *outfile;
   float _ampGain;
   float _ampIBFfrac;
   int _collSyst;
   int _shiftElectrons;
 private:
   double _freqKhz = 22;
   //int _beamxing = 0;
   std::vector<int> _beamxing;
   //std::vector<int> _beamxing_end;

   int _evtstart = 0;
   int _fAvg = 0;
   int _fSliming = 0;
   TTree *_rawHits;
   int _isOnPlane;
   float _hit_z  ;
   float _hit_r  ;
   float _hit_phi;
   float _hit_eion;
   float _ibf_vol;   
   float _amp_ele_vol;
   float _event_timestamp;
   float _event_bunchXing;
   
   bool _fUseIBFMap;
   TH2*   _h_modules_anode;
   TH2*   _h_modules_measuredibf;
   TH1*   _h_hits;
   TH1*   _h_R;
   TH2*   _h_DC_E;
   TH3*   _h_SC_prim[30];
   TH3*   _h_SC_ibf[30];


   float f=0.5;//for now, just pick the middle of the hit.  Do better later.
   float ns=1e-9,us=1e-6,ms=1e-3,s=1;
   float um=1e-3, mm=1, cm=10,m=1000; //changed to make 'cm' 1.0, for convenience.
   float Hz=1, kHz=1e3, MHz=1e6;
   float V=1;
   //used two ways:  1) to apply units to variables when defined
   //                2) to divide by certain units so that those variables are expressed in those units.

   //float ionMobility=3.37*cm*cm/V/s;
   float ionMobility=1.65*cm*cm/V/s;
   float vIon=ionMobility*400*V/cm;
   //float vIon=16.0*um/us;
   float mbRate=_freqKhz*kHz;
   float xingRate = 9.383*MHz;
   float mean = mbRate/xingRate;
   float z_rdo=105.5*cm;
   float rmin=20*cm;
   float rmax=78*cm;

   double Ne_dEdx = 1.56/cm;   // keV/cm
   double CF4_dEdx = 7.00/cm;  // keV/cm
   double Ne_NTotal = 43/cm;    // Number/cm
   double CF4_NTotal = 100/cm;  // Number/cm
   //double Tpc_NTot = 0.90 * Ne_NTotal + 0.10 * CF4_NTotal;
   //double Tpc_dEdx = 0.90 * Ne_dEdx + 0.10 * CF4_dEdx;
   double Tpc_NTot = 0.50 * Ne_NTotal + 0.50 * CF4_NTotal;
   double Tpc_dEdx = 0.50 * Ne_dEdx   + 0.50 * CF4_dEdx;

   double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
   double Tpc_ElectronsPerGeV = Tpc_NTot / Tpc_dEdx*1e6; //electrons per gev.

  int nr=159;
  int nphi=360;
  int nz=62*2;
  
  double hrstep=(rmax-rmin)/nr;
  double hphistep=2*pi/nphi;
  double hzstep=z_rdo/nz;

   int nBeams = z_rdo/(vIon/xingRate); //numaber of beamcrossings to fill TPC


   float _mbRate   ;
   float _xingRate ;
   float _mean     ;
};

#endif // fillSpaceChargeMaps_H
