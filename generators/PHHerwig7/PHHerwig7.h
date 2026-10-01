// PHHerwig7 -- Herwig 7 as an in-process Fun4All generator, the way PHPythia8
// wraps PYTHIA 8.  One Fun4All process, no FIFO, no HepMC on disk; seed from
// PHRandomSeed(); generated-event count and cross section in RUN/PHGenIntegral;
// truth triggers applied inside the generation loop so only triggered events
// reach the node tree.
//
// Input is a Herwig .run file (the output of `Herwig read card.in`, which
// hwrun.sh produces from a config).  Reading the card in-process is deliberately
// not done here: `read` builds the whole repository and is best kept as the
// separate, cached step Herwig designed it to be.
#ifndef PHHERWIG7_H
#define PHHERWIG7_H

#include <fun4all/SubsysReco.h>
#include <phhepmc/PHHepMCGenHelper.h>

#include <hepmctruthtrigger/HepMCTruthTriggerCore.h>

#include <memory>
#include <string>
#include <vector>

class PHCompositeNode;
class PHGenIntegral;
class PHHerwig7Core;

class PHHerwig7 : public SubsysReco, public PHHepMCGenHelper
{
 public:
  explicit PHHerwig7(const std::string &name = "PHHerwig7");
  ~PHHerwig7() override;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void Print(const std::string &what = "ALL") const override;

  // ---- what to run --------------------------------------------------------
  void set_run_file(const std::string &f) { m_RunFile = f; }
  // optional Herwig setup file applied on top of the .run (Herwig --setupfile)
  void set_setup_file(const std::string &f) { m_SetupFile = f; }
  // directories searched for files referenced by the setup file
  void add_read_directory(const std::string &d) { m_ReadDirs.push_back(d); }
  // 0 (default): seed from PHRandomSeed(); >0: use this seed
  void set_seed(int s) { m_Seed = s; }

  // ---- truth trigger, same definitions as hwfilter / HepMCTruthTrigger -----
  void require_jet(bool b = true) { m_Trig.requireJet = b; }
  void set_jet_R(double r) { m_Trig.jetR = r; }
  void set_jet_ptmin(double pt) { m_Trig.jetPtMin = pt; }
  void set_jet_etamax(double eta) { m_Trig.jetEtaMax = eta; }
  void set_leadjet_window(double lo, double hi)
  {
    m_Trig.windowLo = lo;
    m_Trig.windowHi = hi;
  }
  void require_photon(bool b = true) { m_Trig.requirePhoton = b; }
  void set_photon_ptmin(double pt) { m_Trig.photonPtMin = pt; }
  void set_photon_etamax(double eta) { m_Trig.photonEtaMax = eta; }
  void set_photon_isolation(bool on, double r = 0.4, double frac = 0.1)
  {
    m_Trig.photonIso = on;
    m_Trig.photonIsoR = r;
    m_Trig.photonIsoFrac = frac;
  }
  void set_photon_prompt(bool b = true) { m_Trig.photonPrompt = b; }
  void set_and_mode(bool b = true) { m_Trig.andMode = b; }
  // give up after this many consecutive failed attempts (0 = never)
  void set_max_trigger_attempts(long n) { m_MaxAttempts = n; }

  // ---- bookkeeping ---------------------------------------------------------
  void save_integrated_luminosity(bool b) { m_SaveIntegratedLuminosityFlag = b; }
  void save_event_weight(bool b) { m_SaveEventWeightFlag = b; }

  unsigned long get_n_generated() const { return m_NGenerated; }
  unsigned long get_n_passed() const { return m_NPassed; }

 private:
  int create_node_tree(PHCompositeNode *topNode) override;

  std::string m_RunFile;
  std::string m_SetupFile;
  std::vector<std::string> m_ReadDirs;
  int m_Seed = 0;

  HepMCTruthTriggerCore::Config m_Trig;
  long m_MaxAttempts = 0;

  bool m_SaveIntegratedLuminosityFlag = true;
  bool m_SaveEventWeightFlag = true;

  std::unique_ptr<PHHerwig7Core> m_Core;
  PHGenIntegral *m_IntegralNode = nullptr;

  unsigned long m_NGenerated = 0;
  unsigned long m_NPassed = 0;
  double m_SumWAll = 0.0;
  double m_SumWPass = 0.0;
};

#endif
