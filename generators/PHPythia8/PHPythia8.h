#ifndef PHPYTHIA8_PHPYTHIA8_H
#define PHPYTHIA8_PHPYTHIA8_H

#include <fun4all/SubsysReco.h>

#include <phhepmc/PHHepMCGenHelper.h>

#include <algorithm> // for max
#include <cmath>
#include <string>
#include <vector>

class PHCompositeNode;
class PHGenIntegral;
class PHPy8GenTrigger;

namespace HepMC
{
  class Pythia8ToHepMC;
}  // namespace HepMC

namespace Pythia8
{
  class Pythia;
}

class PHPythia8 : public SubsysReco, public PHHepMCGenHelper
{
 public:
  PHPythia8(const std::string &name = "PHPythia8");
  ~PHPythia8() override;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_config_file(const std::string &cfg_file)
  {
    m_ConfigFileName = cfg_file;
  }

  void print_config() const;

  /// set event selection criteria
  void register_trigger(PHPy8GenTrigger *theTrigger);
  void set_trigger_OR()
  {
    m_TriggersOR = true;
    m_TriggersAND = false;
  }  // default true
  void set_trigger_AND()
  {
    m_TriggersAND = true;
    m_TriggersOR = false;
  }

  /// pass commands directly to PYTHIA8
  void process_string(const std::string &s) { m_Commands.push_back(s); }
  void beam_vertex_parameters(double beamX,
                              double beamY,
                              double beamZ,
                              double beamXsigma,
                              double beamYsigma,
                              double beamZsigma)
  {
    set_vertex_distribution_mean(beamX, beamY, beamZ, 0);
    set_vertex_distribution_width(beamXsigma, beamYsigma, beamZsigma, 0);
  }

  void save_event_weight(const bool b) { m_SaveEventWeightFlag = b; }
  void save_integrated_luminosity(const bool b) { m_SaveIntegratedLuminosityFlag = b; }

 private:
  int read_config(const std::string &cfg_file);
  int create_node_tree(PHCompositeNode *topNode) final;
  double percent_diff(const double a, const double b) { return fabs((a - b) / a); }
  int m_EventCount;

  // event selection
  std::vector<PHPy8GenTrigger *> m_RegisteredTriggers;
  bool m_TriggersOR;
  bool m_TriggersAND;

  // PYTHIA
  Pythia8::Pythia *m_Pythia8;

  std::string m_ConfigFileName;
  std::vector<std::string> m_Commands;

  // HepMC
  HepMC::Pythia8ToHepMC *m_Pythia8ToHepMC;

  //! whether to store the overall event weight into the HepMC weights
  bool m_SaveEventWeightFlag;

  //! whether to store the integrated luminosity and other event statistics to the TOP/RUN/PHGenIntegral node
  bool m_SaveIntegratedLuminosityFlag;

  //! pointer to data node saving the integrated luminosity
  PHGenIntegral *m_IntegralNode;
};

#endif /* PHPYTHIA8_PHPYTHIA8_H */
