#ifndef G4INTT_PHG4INTTDIGITIZER_H
#define G4INTT_PHG4INTTDIGITIZER_H

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <map>
#include <string>   // for string
#include <utility>  // for pair
#include <vector>

class PHCompositeNode;

class PHG4InttDigitizer : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4InttDigitizer(const std::string &name = "PHG4InttDigitizer");
  ~PHG4InttDigitizer() override;

  //! run initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  //! end of process
  int End(PHCompositeNode *topNode) override;

  void SetDefaultParameters() override;

  void Detector(const std::string &d) { detector = d; }

  void set_adc_scale(const int &layer, const std::vector<double> &userrange);

 private:
  void CalculateLadderCellADCScale(PHCompositeNode *topNode);

  void DigitizeLadderCells(PHCompositeNode *topNode);

  // noise electrons
  float added_noise();

  std::string detector = "INTT";

  float mNoiseMean = 457.2;        // Mean of noise electron distribution
  float mNoiseSigma = 166.6;       // Sigma of noise electron distribution
  float mEnergyPerPair = 3.62e-9;  // GeV/e-h pair

  // settings
  std::map<int, unsigned int> _max_adc;
  std::map<int, float> _energy_scale;

  // storage
  //SvtxHitMap *_hitmap;

  const unsigned int nadcbins = 8;
  std::map<int, std::vector<std::pair<double, double> > > _max_fphx_adc;

  unsigned int m_nCells = 0;
  unsigned int m_nDeadCells = 0;

  //! random generator that conform with sPHENIX standard
  gsl_rng *RandomGenerator = nullptr;
};

#endif
