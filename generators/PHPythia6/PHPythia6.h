#ifndef PHPYTHIA6_PHPYTHIA6_H
#define PHPYTHIA6_PHPYTHIA6_H

#include <fun4all/SubsysReco.h>

#include <phhepmc/PHHepMCGenHelper.h>

#include <string>
#include <vector>

class PHCompositeNode;
class PHPy6GenTrigger;

class PHPythia6 : public SubsysReco, public PHHepMCGenHelper
{
 public:
  PHPythia6(const std::string &name = "PHPythia6");
  ~PHPythia6() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  void set_config_file(const std::string &cfg_file) { _configFile = cfg_file; }

  void print_config() const;

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

  void save_ascii(const std::string &fname = "pythia_hepmc.dat")
  {
    _save_ascii = true;
    _filename_ascii = fname;
  }

  /// set event selection criteria
  void register_trigger(PHPy6GenTrigger *theTrigger);
  void set_trigger_OR() { _triggersOR = true; }  // default true
  void set_trigger_AND() { _triggersAND = true; }

 private:
  int ReadConfig(const std::string &cfg_file = "");

  /** Certain Pythia switches and parameters only accept integer values
   * This function checks if input values are integers and
   * warns the user if they are not
   */
  void IntegerTest(double number);

  int _eventcount;
  int _geneventcount;

  // Pythia6 configuration file
  std::string _configFile;

  /**
   * Save HepMC event to ASCII file?
   */
  bool _save_ascii;

  /**
   * ASCII file name to save HepMC event to
   */
  std::string _filename_ascii;

  // event selection
  std::vector<PHPy6GenTrigger *> _registeredTriggers;
  bool _triggersOR;
  bool _triggersAND;

  /**
   * definition needed to use pythia wrapper headers from HepMC
   */
  void initPythia();
};

#endif /* PHPYTHIA6_PHPYTHIA6_H */
