// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4MVTX_PHG4MVTXDIGITIZER_H
#define G4MVTX_PHG4MVTXDIGITIZER_H

#include <fun4all/SubsysReco.h>

// rootcint barfs with this header so we need to hide it
#if !defined(__CINT__) || defined(__CLING__)
#include <gsl/gsl_rng.h>
#endif

#include <map>
#include <string>                // for string
#include <utility>               // for pair, make_pair
#include <vector>


class PHCompositeNode;

class PHG4MvtxDigitizer : public SubsysReco
{
 public:
  PHG4MvtxDigitizer(const std::string &name = "PHG4MvtxDigitizer");
  virtual ~PHG4MvtxDigitizer();

  //! module initialization
  int Init(PHCompositeNode *topNode) { return 0; }

  //! run initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  //! end of process
  int End(PHCompositeNode *topNode) { return 0; };

  void set_adc_scale(const int layer, const unsigned int max_adc, const float energy_per_adc)
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
  std::map<int, unsigned int> _max_adc;
  std::map<int, float> _energy_scale;
  float _energy_threshold;

#if !defined(__CINT__) || defined(__CLING__)
  //! random generator that conform with sPHENIX standard
  gsl_rng *RandomGenerator;
#endif
};

#endif
