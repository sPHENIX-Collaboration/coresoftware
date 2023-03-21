#ifndef G4BBC_BBCSIMRECO_H
#define G4BBC_BBCSIMRECO_H

#include <fun4all/SubsysReco.h>

#include <Rtypes.h>

#include <gsl/gsl_rng.h>

#include <cmath>
#include <map>
#include <string>

// Forward declarations
class PHCompositeNode;
class PHG4HitContainer;
class PHG4TruthInfoContainer;
class EventHeader;
class BbcOut;
class BbcPmtContainer;
class TDatabasePDG;
class TRandom3;
class TH1;
class TH2;
class TF1;

class BbcSimReco : public SubsysReco
{
 public:
  // Default constructor
  BbcSimReco(const std::string &name = "BbcSimReco");

  ~BbcSimReco() override;

  //! Initialization, called for at overall initialization
  int Init(PHCompositeNode *) override;

  //! Initialization at start of every run
  int InitRun(PHCompositeNode *) override;

  //! Process Event, called for each event
  int process_event(PHCompositeNode *) override;

  //! Reset after every event
  // int ResetEvent(PHCompositeNode * /*topNode*/) override;

  //! End, write and close files
  int End(PHCompositeNode *) override { return 0; }

  //! Set time resolution (each channel has same time resol)
  void set_tres(const Float_t tr) { _tres = tr; }

 private:
  void CreateNodes(PHCompositeNode *topNode);  // Create all the nodes
  void GetNodes(PHCompositeNode *);            // Get all the needed nodes

  Int_t f_evt = 0;
  Float_t f_vx = NAN;
  Float_t f_vy = NAN;
  Float_t f_vz = NAN;
  Float_t f_vt = NAN;
  Float_t f_pmtq[128]{};   // npe in each arm
  Float_t f_pmtt0[128]{};  // time in each arm
  Float_t f_pmtt1[128]{};  // time in each arm
  Short_t f_bbcn[2]{};     // num hits for each arm (north and south)
  Float_t f_bbcq[2]{};     // total charge (currently npe) in each arm
  Float_t f_bbct[2]{};     // time in arm
  Float_t f_bbcte[2]{};    // earliest hit time in arm
  Float_t f_bbcz = NAN;    // z-vertex
  Float_t f_bbct0 = NAN;   // start time

  TH1 *hevt_bbct[2]{};  // time in each bbc, per event
  TF1 *gaussian = nullptr;

  //
  TDatabasePDG *_pdg = nullptr;
  gsl_rng *m_RandomGenerator = nullptr;
  unsigned int m_Seed = 0;
  Float_t _tres = NAN;  // time resolution of one channel

  std::map<int, int> _pids;  // PIDs of tracks in the BBC

  // Input Objects from DST
  PHG4TruthInfoContainer *_truth_container = nullptr;
  PHG4HitContainer *_bbchits = nullptr;
  EventHeader *_evtheader = nullptr;

  // Output to DST
  BbcOut *_bbcout = nullptr;
  BbcPmtContainer *_bbcpmts = nullptr;
};

#endif  //* __BBCSIMRECO_H__ *//
