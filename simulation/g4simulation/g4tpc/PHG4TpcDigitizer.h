#ifndef G4TPC_PHG4TpcDIGITIZER_H
#define G4TPC_PHG4TpcDIGITIZER_H

#include <fun4all/SubsysReco.h>
#include <g4detectors/PHG4CellDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>

#include <map>
#include <vector>

// rootcint barfs with this header so we need to hide it
#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

class SvtxHitMap;
class PHG4Cell;
class TrkrHit;

class PHG4TpcDigitizer : public SubsysReco
{
 public:
  PHG4TpcDigitizer(const std::string &name = "PHG4TpcDigitizer");
  virtual ~PHG4TpcDigitizer() {}

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

  void SetTPCMinLayer(const int minlayer) { TPCMinLayer = minlayer; };
  void SetADCThreshold(const float thresh) { ADCThreshold = thresh; };
  void SetENC(const float enc) { TPCEnc = enc; };

 private:
  void CalculateCylinderCellADCScale(PHCompositeNode *topNode);
  void DigitizeCylinderCells(PHCompositeNode *topNode);
  void PrintHits(PHCompositeNode *topNode);
  float added_noise();

  unsigned int TPCMinLayer;
  float ADCThreshold;
  float TPCEnc;
  float Pedestal;
  float ChargeToPeakVolts;

  float ADCSignalConversionGain;
  float ADCNoiseConversionGain;

  /*
  std::vector<std::vector<const PHG4Cell *> > layer_sorted_cells;
  std::vector<std::vector<const PHG4Cell *> > phi_sorted_cells;
  std::vector<std::vector<const PHG4Cell *> > z_sorted_cells;
  */
  std::vector<std::vector<TrkrHitSet::ConstIterator > > layer_sorted_hits;
  std::vector<std::vector<TrkrHitSet::ConstIterator> > phi_sorted_hits;
  std::vector<std::vector<TrkrHitSet::ConstIterator> > z_sorted_hits;

  std::vector<float> adc_input;
  //std::vector<PHG4CellDefs::keytype> adc_cellid;
  std::vector<TrkrDefs::hitkey> adc_hitid;
  std::vector<int> is_populated;

  // settings
  std::map<int, unsigned int> _max_adc;
  std::map<int, float> _energy_scale;

  // storage
  SvtxHitMap *_hitmap;

#ifndef __CINT__
  //! random generator that conform with sPHENIX standard
  gsl_rng *RandomGenerator;
#endif
};

#endif
