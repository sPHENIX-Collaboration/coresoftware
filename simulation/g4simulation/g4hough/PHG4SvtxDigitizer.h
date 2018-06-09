#ifndef __PHG4SVTXDIGITIZER_H__
#define __PHG4SVTXDIGITIZER_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

#include <vector>

// rootcint barfs with this header so we need to hide it
#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

class SvtxHitMap;

class PHG4SvtxDigitizer : public SubsysReco
{
 public:

  PHG4SvtxDigitizer(const std::string &name = "PHG4SvtxDigitizer");
  virtual ~PHG4SvtxDigitizer(){}
  
  //! module initialization
  int Init(PHCompositeNode *topNode){return 0;}
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode) {return 0;};
  
  void set_adc_scale(const int layer, const unsigned int max_adc, const float energy_per_adc) {
    _max_adc.insert(std::make_pair(layer,max_adc));
    _energy_scale.insert(std::make_pair(layer,energy_per_adc));
  }
    
  void SetTPCMinLayer(const int minlayer) {TPCMinLayer = minlayer;};
  void SetADCThreshold(const float thresh) {ADCThreshold = thresh;};
  void SetENC(const float enc){TPCEnc = enc;};
    
 private:

  void CalculateCylinderCellADCScale(PHCompositeNode *topNode);
  void CalculateLadderCellADCScale(PHCompositeNode *topNode);
  void CalculateMapsLadderCellADCScale(PHCompositeNode *topNode);

  void DigitizeCylinderCells(PHCompositeNode *topNode);
  void DigitizeLadderCells(PHCompositeNode *topNode);
  void DigitizeMapsLadderCells(PHCompositeNode *topNode);
  void PrintHits(PHCompositeNode *topNode);
  float added_noise();

  unsigned int TPCMinLayer;  
  float ADCThreshold;
  float TPCEnc;

  // settings
  std::map<int,unsigned int> _max_adc;
  std::map<int,float> _energy_scale;

  // storage
  SvtxHitMap* _hitmap;
  
  PHTimeServer::timer _timer;   ///< Timer


#ifndef __CINT__
  //! random generator that conform with sPHENIX standard
  gsl_rng *RandomGenerator;
#endif

};

#endif
