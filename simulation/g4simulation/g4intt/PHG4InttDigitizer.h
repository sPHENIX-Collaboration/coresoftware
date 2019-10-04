#ifndef G4INTT_PHG4INTTDIGITIZER_H
#define G4INTT_PHG4INTTDIGITIZER_H

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

// rootcint barfs with this header so we need to hide it
#if !defined(__CINT__) || defined(__CLING__)
#include <gsl/gsl_rng.h>
#endif

#include <map>
#include <string>                              // for string
#include <utility>                             // for pair
#include <vector>

class PHCompositeNode;

class PHG4InttDigitizer : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4InttDigitizer(const std::string &name = "PHG4InttDigitizer");
  virtual ~PHG4InttDigitizer();

  //! run initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  //! end of process
  int End(PHCompositeNode *topNode);

  void SetDefaultParameters();

  void Detector(const std::string &d) { detector = d; }

  void set_adc_scale(const int &layer, const std::vector<double> &userrange);

 private:
  void CalculateLadderCellADCScale(PHCompositeNode *topNode);

  void DigitizeLadderCells(PHCompositeNode *topNode);

  std::string detector;
  // noise electrons
  float added_noise();

  float mNoiseMean;      // Mean of noise electron distribution
  float mNoiseSigma;     // Sigma of noise electron distribution
  float mEnergyPerPair;  // GeV/e-h pair

  // settings
  std::map<int, unsigned int> _max_adc;
  std::map<int, float> _energy_scale;

  // storage
  //SvtxHitMap *_hitmap;

  const unsigned int nadcbins = 8;
  std::map<int, std::vector<std::pair<double, double> > > _max_fphx_adc;

  unsigned int m_nCells;
  unsigned int m_nDeadCells;

#if !defined(__CINT__) || defined(__CLING__)
  //! random generator that conform with sPHENIX standard
  gsl_rng *RandomGenerator;
#endif
};

#endif
