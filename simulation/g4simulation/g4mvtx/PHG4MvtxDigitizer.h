/**
 * @file g4mvtx/PHG4MvtxDigitizer.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Digitization class for simulated mvtx hits
 */
#ifndef G4MVTX_PHG4MVTXDIGITIZER_H
#define G4MVTX_PHG4MVTXDIGITIZER_H

#include <fun4all/SubsysReco.h>

#include <map>

class TrkrHitSetContainer;

/**
 * @brief Digitization class for simulated mvtx hits
 */
class PHG4MvtxDigitizer : public SubsysReco
{
 public:
  PHG4MvtxDigitizer(const std::string &name = "PHG4MvtxDigitizer");
  virtual ~PHG4MvtxDigitizer() {}

  //! module initialization
  int Init(PHCompositeNode *topNode) { return 0; }

  //! run initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  //! end of process
  int End(PHCompositeNode *topNode) { return 0; }

  void set_adc_scale(const int layer, const unsigned int max_adc, const float energy_per_adc)
  {
    m_maxAdc.insert(std::make_pair(layer, max_adc));
    m_energyScale.insert(std::make_pair(layer, energy_per_adc));
  }

 private:
  void CalculateADCScale(PHCompositeNode *topNode);

  void DigitizeCells(PHCompositeNode *topNode);
  void PrintHits(PHCompositeNode *topNode);

  // settings
  std::map<int, unsigned int> m_maxAdc;
  std::map<int, float> m_energyScale;

  // storage
  TrkrHitSetContainer *m_hitsets;
};

#endif  //G4MVTX_PHG4MVTXDIGITIZER_H
