#ifndef PHPYTHIA8_PHPYTHIA8_H
#define PHPYTHIA8_PHPYTHIA8_H

#include <fun4all/SubsysReco.h>

#include <phhepmc/PHHepMCGenHelper.h>

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

class PHPythia8 : public SubsysReco
{
 public:
  PHPythia8(const std::string &name = "PHPythia8");
  virtual ~PHPythia8();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

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

  //! toss a new vertex according to a Uniform or Gaus distribution
  void set_vertex_distribution_function(PHHepMCGenHelper::VTXFUNC x, PHHepMCGenHelper::VTXFUNC y, PHHepMCGenHelper::VTXFUNC z, PHHepMCGenHelper::VTXFUNC t)
  {
    m_HepMC_Helper.set_vertex_distribution_function(x, y, z, t);
  }

  //! set the mean value of the vertex distribution, use PHENIX units of cm, ns
  void set_vertex_distribution_mean(const double x, const double y, const double z, const double t)
  {
    m_HepMC_Helper.set_vertex_distribution_mean(x, y, z, t);
  }

  //! set the width of the vertex distribution function about the mean, use PHENIX units of cm, ns
  void set_vertex_distribution_width(const double x, const double y, const double z, const double t)
  {
    m_HepMC_Helper.set_vertex_distribution_width(x, y, z, t);
  }
  //
  //! reuse vertex from another PHHepMCGenEvent with embedding_id = src_embedding_id Additional smearing and shift possible with set_vertex_distribution_*()
  void set_reuse_vertex(int src_embedding_id)
  {
    m_HepMC_Helper.set_reuse_vertex(src_embedding_id);
  }

  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int get_embedding_id() const { return m_HepMC_Helper.get_embedding_id(); }
  //
  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  void set_embedding_id(int id) { m_HepMC_Helper.set_embedding_id(id); }
  //! whether to store the integrated luminosity and other event statistics to the TOP/RUN/PHGenIntegral node
  void save_integrated_luminosity(const bool b) { m_SaveIntegratedLuminosityFlag = b; }

 private:
  int read_config(const std::string &cfg_file);
  int create_node_tree(PHCompositeNode *topNode);
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

  //! helper for insert HepMC event to DST node and add vertex smearing
  PHHepMCGenHelper m_HepMC_Helper;

  //! whether to store the integrated luminosity and other event statistics to the TOP/RUN/PHGenIntegral node
  bool m_SaveIntegratedLuminosityFlag;

  //! pointer to data node saving the integrated luminosity
  PHGenIntegral *m_IntegralNode;
};

#endif /* PHPYTHIA8_PHPYTHIA8_H */
