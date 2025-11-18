#ifndef G4MBD_MBDDIGITIZATION_H
#define G4MBD_MBDDIGITIZATION_H

#include <mbd/MbdDefs.h>

#include <fun4all/SubsysReco.h>

#include <Rtypes.h>

#include <gsl/gsl_rng.h>

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <string>

// Forward declarations
class PHCompositeNode;
class PHG4HitContainer;
class PHG4TruthInfoContainer;
class EventHeader;
class MbdPmtContainer;
class TDatabasePDG;
class TRandom3;
class TH1;
class TH2;
class TF1;

class MbdDigitization : public SubsysReco
{
 public:
  // Default constructor
  MbdDigitization(const std::string &name = "MbdDigitization");

  ~MbdDigitization() override;

  //! Initialization, called for at overall initialization
  int Init(PHCompositeNode *) override;

  //! Initialization at start of every run
  int InitRun(PHCompositeNode *) override;

  //! Process Event, called for each event
  int process_event(PHCompositeNode *) override;

  //! Set time resolution (each channel has same time resol)
  void set_tres(const Float_t tr) { _tres = tr; }

 private:
  void CreateNodes(PHCompositeNode *topNode);  // Create all the nodes
  void GetNodes(PHCompositeNode *);            // Get all the needed nodes

  gsl_rng *m_RandomGenerator{nullptr};

  // Output to DST
  MbdPmtContainer *_bbcpmts{nullptr};
  // Input Objects from DST
  PHG4TruthInfoContainer *_truth_container{nullptr};
  PHG4HitContainer *_bbchits{nullptr};

  TDatabasePDG *_pdg{nullptr};
  TF1 *gaussian{nullptr};

  Float_t f_vx{std::numeric_limits<Float_t>::quiet_NaN()};
  Float_t f_vy{std::numeric_limits<Float_t>::quiet_NaN()};
  Float_t f_vz{std::numeric_limits<Float_t>::quiet_NaN()};
  Float_t f_vt{std::numeric_limits<Float_t>::quiet_NaN()};
  Float_t _tres{0.05};  // time resolution of one channel

  std::array<Float_t, MbdDefs::MBD_N_PMT> f_pmtq{};    // equiv. nch in each pmt
  std::array<Float_t, MbdDefs::MBD_N_PMT> f_pmtt0{};   // time in each pmt
  std::array<Float_t, MbdDefs::MBD_N_PMT> f_pmtt1{};   // time in each pmt
  std::array<Float_t, MbdDefs::MBD_N_PMT> f_pmtnpe{};  // npe in each pmt

  // gains
  std::array<Float_t, MbdDefs::MBD_N_PMT> _gains = {};

  unsigned int m_Seed{0};

  std::map<int, int> _pids;  // PIDs of tracks in the BBC
};

#endif  //* G4MBD_MBDDIGITIZATION_H *//
