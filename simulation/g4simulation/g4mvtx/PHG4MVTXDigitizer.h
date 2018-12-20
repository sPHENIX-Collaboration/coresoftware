#ifndef __PHG4MVTXDIGITIZER_H__
#define __PHG4MVTXDIGITIZER_H__

#include <fun4all/SubsysReco.h>
#include <g4detectors/PHG4CellDefs.h>
#include <phool/PHTimeServer.h>

#include <map>
#include <vector>

// rootcint barfs with this header so we need to hide it
#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

class SvtxHitMap;
class PHG4Cell;

class PHG4MVTXDigitizer : public SubsysReco
{
 public:
  PHG4MVTXDigitizer(const std::string &name = "PHG4MVTXDigitizer");
  virtual ~PHG4MVTXDigitizer() {}

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

 private:
  void CalculateMVTXLadderCellADCScale(PHCompositeNode *topNode);
  void DigitizeMVTXLadderCells(PHCompositeNode *topNode);
  void PrintHits(PHCompositeNode *topNode);

  std::vector<std::vector<const PHG4Cell *> > layer_sorted_cells;
  std::vector<std::vector<const PHG4Cell *> > phi_sorted_cells;
  std::vector<std::vector<const PHG4Cell *> > z_sorted_cells;
  std::vector<float> adc_input;
  std::vector<PHG4CellDefs::keytype> adc_cellid;
  std::vector<int> is_populated;

  // settings
  std::map<int, unsigned int> _max_adc;
  std::map<int, float> _energy_scale;

  // storage
  SvtxHitMap *_hitmap;

  PHTimeServer::timer _timer;  ///< Timer

#ifndef __CINT__
  //! random generator that conform with sPHENIX standard
  gsl_rng *RandomGenerator;
#endif
};

#endif
