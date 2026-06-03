// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4TPC_PHG4TPCDIGITIZER_H
#define G4TPC_PHG4TPCDIGITIZER_H

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>

#include <gsl/gsl_rng.h>

#include <map>
#include <string>   // for string
#include <utility>  // for pair, make_pair
#include <vector>

class PHCompositeNode;
class TrkrHit;

class PHG4TpcDigitizer : public SubsysReco
{
 public:
  PHG4TpcDigitizer(const std::string &name = "PHG4TpcDigitizer");
  ~PHG4TpcDigitizer() override;

  //! run initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  void set_adc_scale(const int layer, const unsigned int max_adc, const float energy_per_adc)
  {
    _max_adc.insert(std::make_pair(layer, max_adc));
    _energy_scale.insert(std::make_pair(layer, energy_per_adc));
  }

  void SetTpcMinLayer(const int minlayer) { TpcMinLayer = minlayer; };
  void SetADCThreshold(const float thresh) { ADCThreshold = thresh; };
  void SetENC(const float enc) { TpcEnc = enc; };
  void set_skip_noise_flag(const bool skip) { skip_noise = skip; }

 private:
  void CalculateCylinderCellADCScale(PHCompositeNode *topNode);
  void DigitizeCylinderCells(PHCompositeNode *topNode);
  float added_noise();
  float add_noise_to_bin(float signal);

  unsigned int TpcMinLayer {7};
  unsigned int TpcNLayers {48};
  float ADCThreshold {2700};
  float ADCThreshold_mV {0};
  float TpcEnc {670};
  float Pedestal {50000};
  float ChargeToPeakVolts {20};
  float ADCSignalConversionGain {std::numeric_limits<float>::quiet_NaN()};
  float ADCNoiseConversionGain {std::numeric_limits<float>::quiet_NaN()};

  bool skip_noise {false};

  std::vector<std::vector<TrkrHitSet::ConstIterator> > phi_sorted_hits;
  std::vector<float> adc_input;
  std::vector<TrkrHit *> signal_hit_by_tbin;

  // settings
  std::map<int, unsigned int> _max_adc;
  std::map<int, float> _energy_scale;

  //! random generator that conform with sPHENIX standard
  gsl_rng *RandomGenerator {nullptr};
};

#endif
