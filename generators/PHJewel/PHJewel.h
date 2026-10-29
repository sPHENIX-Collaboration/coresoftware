// PHJewel -- JEWEL 2.6.0 as an in-process Fun4All generator, the way
// PHPythia8 wraps PYTHIA 8 and PHHerwig7 wraps Herwig 7.  One Fun4All process,
// no FIFO, no HepMC on disk; seed from PHRandomSeed(); generated-event count,
// weight sum and cross section in RUN/PHGenIntegral; truth triggers applied
// inside the generation loop so only triggered events reach the node tree.
//
// Input is a JEWEL card pair (main card + medium card) plus the medium type,
// which selects the core library (libjewel_simple / _brick / _vac).  Cards are
// read strictly (JewelCard.h); what JEWEL is actually given is a canonical
// copy written next to the log, and its hash goes into PHGenIntegral.
//
// One PHJewel per process: JEWEL keeps all its state in Fortran common blocks.
// Do not combine with PHPythia6 in the same process (shared Fortran runtime
// unit numbers), although the PYTHIA 6 symbols themselves are isolated.
#ifndef PHJEWEL_H
#define PHJEWEL_H

#include <fun4all/SubsysReco.h>
#include <phhepmc/PHHepMCGenHelper.h>

#include <memory>
#include <string>
#include <vector>

class PHCompositeNode;
class PHGenIntegral;
class PHJewelGenTrigger;
class JewelCore;

class PHJewel : public SubsysReco, public PHHepMCGenHelper
{
 public:
  explicit PHJewel(const std::string &name = "PHJewel");
  ~PHJewel() override;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void Print(const std::string &what = "ALL") const override;

  // ---- what to run --------------------------------------------------------
  // main JEWEL card; relative paths inside it are relative to its directory
  void set_card(const std::string &f) { m_Card = f; }
  // simple | brick | vac -- must match the medium card's keys
  void set_medium(const std::string &m) { m_Medium = m; }
  // 0 (default): seed from PHRandomSeed(); >0: this PYTHIA seed (MRPY(1)).
  // A standalone run with NJOB = k reproduces set_seed(1000*k) exactly
  // (fixed ISOCHANNEL).
  void set_seed(int s) { m_Seed = s; }
  // output stem for <stem>.log, <stem>.card.dat, <stem>.medium.dat;
  // default "PHJewel_<seed>" in the working directory
  void set_output_stem(const std::string &s) { m_Stem = s; }
  // override the core library (default libjewel_<medium>.so on LD_LIBRARY_PATH)
  void set_core_library(const std::string &l) { m_CoreLib = l; }
  // validation: JEWEL's own Fortran HepMC writer writes every JEWEL event here
  void set_fortran_hepmc_file(const std::string &f) { m_FortranHepMC = f; }
  // JEWEL's EVINFOFILE (one line per JEWEL event; event number = HepMC event number)
  void set_evinfo_file(const std::string &f) { m_EvInfo = f; }
  // allow JEWEL to integrate missing tables (slow; never in parallel jobs)
  void allow_table_generation(bool b = true) { m_AllowTables = b; }
  // ISOCHANNEL XX: JEWEL events per isospin-channel draw (PYTHIA re-init on change)
  void set_isospin_block(long n) { m_IsospinBlock = n; }

  // ---- truth triggers, as PHPythia8 / PHHerwig7 -----------------------------
  // PHJewelJetTrigger, PHJewelPhotonTrigger; applied inside the generation
  // loop, so only passing events reach the node tree.  Not owned by PHJewel.
  // Note: with KEEPRECOILS T the jet trigger sees UNSUBTRACTED jets (see
  // PHJewelJetTrigger.h).
  void register_trigger(PHJewelGenTrigger *trigger);
  void set_trigger_OR() { m_TriggersOR = true; }  // default
  void set_trigger_AND() { m_TriggersOR = false; }
  // give up after this many consecutive failed attempts (0 = never)
  void set_max_trigger_attempts(long n) { m_MaxAttempts = n; }

  // ---- bookkeeping ---------------------------------------------------------
  void save_integrated_luminosity(bool b) { m_SaveIntegratedLuminosityFlag = b; }

  unsigned long get_n_generated() const { return m_NGenerated; }
  unsigned long get_n_passed() const { return m_NPassed; }

 private:
  int create_node_tree(PHCompositeNode *topNode) override;

  std::string m_Card;
  std::string m_Medium;
  int m_Seed = 0;
  std::string m_Stem;
  std::string m_CoreLib;
  std::string m_FortranHepMC;
  std::string m_EvInfo;
  bool m_AllowTables = false;
  long m_IsospinBlock = 1000;
  std::string m_Hash;

  std::vector<PHJewelGenTrigger *> m_RegisteredTriggers;
  bool m_TriggersOR = true;
  long m_MaxAttempts = 0;

  bool m_SaveIntegratedLuminosityFlag = true;

  std::unique_ptr<JewelCore> m_Core;
  PHGenIntegral *m_IntegralNode = nullptr;

  unsigned long m_NGenerated = 0;  // events JEWEL wrote (before the trigger)
  unsigned long m_NPassed = 0;
  double m_SumWAll = 0.0;
  double m_SumWPass = 0.0;

  static bool s_Instance;
};

#endif
