// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4MVTX_PHG4MVTXDIGITIZER_H
#define G4MVTX_PHG4MVTXDIGITIZER_H

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <map>
#include <string>   // for string
#include <utility>  // for pair, make_pair
#include <vector>

class PHCompositeNode;

class PHG4MvtxDigitizer : public SubsysReco
{
 public:
  PHG4MvtxDigitizer(const std::string &name = "PHG4MvtxDigitizer");
  ~PHG4MvtxDigitizer() override;

  //! module initialization
  int Init(PHCompositeNode * /*topNode*/) override { return 0; }

  //! run initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  //! end of process
  int End(PHCompositeNode * /*topNode*/) override { return 0; };

  void set_adc_scale(const int layer, const unsigned short max_adc, const float energy_per_adc)
  {
    _max_adc.insert(std::make_pair(layer, max_adc));
    _energy_scale.insert(std::make_pair(layer, energy_per_adc));
  }

  void set_energy_threshold(const float threshold)
  {
    _energy_threshold = threshold;
  }

  float get_energy_threshold() { return _energy_threshold; }

 private:
  void CalculateMvtxLadderCellADCScale(PHCompositeNode *topNode);
  void DigitizeMvtxLadderCells(PHCompositeNode *topNode);

  std::vector<float> adc_input;
  std::vector<int> is_populated;

  // settings
  std::map<int, unsigned short> _max_adc;
  std::map<int, float> _energy_scale;
  float _energy_threshold;

  //! random generator that conform with sPHENIX standard
  gsl_rng *RandomGenerator;
};

#endif
