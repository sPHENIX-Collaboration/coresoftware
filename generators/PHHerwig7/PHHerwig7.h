#ifndef PHHERWIG7_H
#define PHHERWIG7_H

#include <fun4all/SubsysReco.h>
#include <phhepmc/PHHepMCGenHelper.h>

#include <memory>
#include <string>
#include <vector>

class PHCompositeNode;
class PHGenIntegral;
class PHHerwig7Core;
class PHHerwig7GenTrigger;

class PHHerwig7 : public SubsysReco, public PHHepMCGenHelper
{
 public:
  explicit PHHerwig7(const std::string &name = "PHHerwig7");
  ~PHHerwig7() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void Print(const std::string &what = "ALL") const override;

  // set the Herwig config file (this is the .in file read by `Herwig read`)
  void set_config_file(const std::string &f) { m_ConfigFile = f; }

  // (optional) corresponds to "Herwig run --setupfile"; example use: one base card plus a per-job setting file
  // This allows changing settings on top of a .run without regenerating it (.run is not produced within Fun4All)
  void set_setup_file(const std::string &f) { m_SetupFile = f; }

  // "Herwig read -I <dir>", adds a directory to the search path for read commands inside a card
  void add_read_directory(const std::string &d) { m_ReadDirs.push_back(d); }

  // set the seed. default is 0, i.e seed from PHRandomSeed(). For reproducibility, set a positive integer
  void set_seed(int s) { m_Seed = s; }

  // give up after this many consecutive failed attempts (0 = never)
  void set_max_trigger_attempts(long n) { m_MaxAttempts = n; }

  // bookkeeping
  void save_integrated_luminosity(bool b) { m_SaveIntegratedLuminosityFlag = b; }
  void save_event_weight(bool b) { m_SaveEventWeightFlag = b; }

  unsigned long get_n_generated() const { return m_NGenerated; }
  unsigned long get_n_passed() const { return m_NPassed; }

  // public setters to register trigger and logic as PHPythia8
  void register_trigger(PHHerwig7GenTrigger *trigger);
  void set_trigger_OR()
  {
    set_trigger_logic(true);
  }
  void set_trigger_AND()
  {
    set_trigger_logic(false);
  }

 private:
  void set_trigger_logic(bool triggerOR)
  {
    m_TriggersOR = triggerOR;
    m_TriggersAND = !triggerOR;
  }

  int create_node_tree(PHCompositeNode *topNode) override;

  std::string m_ConfigFile;
  std::string m_SetupFile;
  std::vector<std::string> m_ReadDirs;
  int m_Seed = 0;

  long m_MaxAttempts = 0;

  bool m_SaveIntegratedLuminosityFlag = true;
  bool m_SaveEventWeightFlag = true;

  std::unique_ptr<PHHerwig7Core> m_Core;
  PHGenIntegral *m_IntegralNode = nullptr;

  unsigned long m_NGenerated = 0;
  unsigned long m_NPassed = 0;
  double m_SumWAll = 0.0;
  double m_SumWPass = 0.0;

  // event selection
  std::vector<PHHerwig7GenTrigger *> m_RegisteredTriggers;
  bool m_TriggersOR{true};
  bool m_TriggersAND{false};
};

#endif
