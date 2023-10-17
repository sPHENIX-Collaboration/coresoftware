#ifndef G4BBC_BBCDIGITIZATION_H
#define G4BBC_BBCDIGITIZATION_H

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
class BbcPmtInfoContainerV1;
class TDatabasePDG;
class TRandom3;
class TH1;
class TH2;
class TF1;

class BbcDigitization : public SubsysReco
{
 public:
  // Default constructor
  BbcDigitization(const std::string &name = "BbcDigitization");

  ~BbcDigitization() override;

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

  Float_t f_vx = NAN;
  Float_t f_vy = NAN;
  Float_t f_vz = NAN;
  Float_t f_vt = NAN;
  Float_t f_pmtq[128]{};   // npe in each arm
  Float_t f_pmtt0[128]{};  // time in each arm
  Float_t f_pmtt1[128]{};  // time in each arm

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

  // Output to DST
  BbcPmtInfoContainerV1 *_bbcpmts = nullptr;
};

#endif  //* __BBCDIGITIZATION_H__ *//
