// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4TPC_PHG4TPCDIGITIZER_H
#define G4TPC_PHG4TPCDIGITIZER_H

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>

#include <map>
#include <string>                  // for string
#include <utility>                 // for pair, make_pair
#include <vector>

// rootcint barfs with this header so we need to hide it
#if !defined(__CINT__) || defined(__CLING__)
#include <gsl/gsl_rng.h>
#endif

class PHCompositeNode;
class SvtxHitMap;

class PHG4TpcDigitizer : public SubsysReco
{
 public:
  PHG4TpcDigitizer(const std::string &name = "PHG4TpcDigitizer");
  virtual ~PHG4TpcDigitizer();

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

  void SetTpcMinLayer(const int minlayer) { TpcMinLayer = minlayer; };
  void SetADCThreshold(const float thresh) { ADCThreshold = thresh; };
  void SetENC(const float enc) { TpcEnc = enc; };

 private:
  void CalculateCylinderCellADCScale(PHCompositeNode *topNode);
  void DigitizeCylinderCells(PHCompositeNode *topNode);
  float added_noise();

  unsigned int TpcMinLayer;
  float ADCThreshold;
  float TpcEnc;
  float Pedestal;
  float ChargeToPeakVolts;

  float ADCSignalConversionGain;
  float ADCNoiseConversionGain;

  std::vector<std::vector<TrkrHitSet::ConstIterator> > phi_sorted_hits;
  std::vector<std::vector<TrkrHitSet::ConstIterator> > z_sorted_hits;

  std::vector<float> adc_input;
  std::vector<TrkrDefs::hitkey> adc_hitid;
  std::vector<int> is_populated;

  // settings
  std::map<int, unsigned int> _max_adc;
  std::map<int, float> _energy_scale;

#if !defined(__CINT__) || defined(__CLING__)
  //! random generator that conform with sPHENIX standard
  gsl_rng *RandomGenerator;
#endif
};

#endif
