// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef READDIGITALCURRENTS_READDIGITALCURRENTS_H
#define READDIGITALCURRENTS_READDIGITALCURRENTS_H

#include <fun4all/SubsysReco.h>

#include <fstream>
#include <map>
#include <string>
#include <vector>

class Fun4AllHistoManager;
class PHCompositeNode;

// class PHG4CylinderCellGeom;

// class TFile;
class TH1;
class TH2;
class TH3;

class readDigitalCurrents : public SubsysReco
{
 public:
  readDigitalCurrents(const std::string &name = "readDigitalCurrents", const std::string &filename = "DC_Hist_OnPlane_WIBF.root");

  virtual ~readDigitalCurrents();

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
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  void SetBeamXing(const std::vector<int> &beamXs);
  void SetEvtStart(int newEvtStart);
  void FillCSV(int fillCSVFile);
  void SetCollSyst(int coll_syst = 0);
  void SetIBF(double ampIBFfrac = 0.004);
  void SetCCGC(double f_ccgc = 0);

  // double pi = 3.14159265358979323846;//2 * acos(0.0);

 protected:
  Fun4AllHistoManager *hm{nullptr};
  std::string _filename;
  // TFile *outfile {nullptr};
  std::map<int, int> _timestamps;
  std::vector<int> _keys;
  double _ampIBFfrac{0.02};
  int _collSyst{0};
  std::ofstream myCSVFile;

 private:
  std::vector<int> _beamxing;
  int _evtstart{0};

  int _fillCSVFile{0};

  int _f_ccgc{0};

  TH2 *_h_modules_measuredibf{nullptr};

  TH1 *_h_R{nullptr};
  TH1 *_h_hits{nullptr};
  TH3 *_h_DC_SC{nullptr};
  TH2 *_h_DC_SC_XY{nullptr};
  TH2 *_h_hit_XY{nullptr};
  TH2 *_h_DC_E{nullptr};
  // TH3*   _h_SC_ibf {nullptr};
  static const int nFrames{30};
  TH3 *_h_SC_ibf[nFrames]{};

  float _event_timestamp{0};
  float _event_bunchXing{0};

  // double pi = 2 * acos(0.0);
  double adc_pedestal{0.};  // 74.4;
  double cm{1e1};
  double m{1e3};
  double mm{1};  // changed to make 'm' 1.0, for convenience.
  //    float ns=1e-9,us=1e-6,ms=1e-3,s=1;
  double ns{1e-9};
  double s{1};
  double V{1};
  // double ionMobility=3.37*cm*cm/V/s;
  double ionMobility{1.65 * cm * cm / V / s};
  double vIon{ionMobility * 400 * V / cm};

  //    float f=0.5;//for now, just pick the middle of the hit.  Do better later.
};

#endif  // READDIGITALCURRENTS_H
