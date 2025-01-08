#ifndef CALORECO_CALOWAVEFORMFITTING_H
#define CALORECO_CALOWAVEFORMFITTING_H

#include <string>
#include <vector>

class TProfile;

class CaloWaveformFitting
{
 public:
  CaloWaveformFitting() = default;
  ~CaloWaveformFitting();

  void set_template_file(const std::string &template_input_file)
  {
    m_template_input_file = template_input_file;
    return;
  }

  void set_nthreads(int nthreads)
  {
    _nthreads = nthreads;
    return;
  }

  void set_softwarezerosuppression(bool usezerosuppression, int softwarezerosuppression)
  {
    _nsoftwarezerosuppression = softwarezerosuppression;
    _bdosoftwarezerosuppression = usezerosuppression;
  }
  void set_maxsoftwarezerosuppression(bool usezerosuppression, int softwarezerosuppression)
  {
    _nsoftwarezerosuppression = softwarezerosuppression;
    _maxsoftwarezerosuppression = usezerosuppression;
  }

  int get_nthreads()
  {
    return _nthreads;
  }
  void set_timeFitLim(float low,float high)
  {
    m_setTimeLim = true;
    m_timeLim_low = low;
    m_timeLim_high = high;
    return;
  }

  void set_bitFlipRecovery(bool dobitfliprecovery) {
    _dobitfliprecovery = dobitfliprecovery;
  }

  std::vector<std::vector<float>> process_waveform(std::vector<std::vector<float>> waveformvector);
  std::vector<std::vector<float>> calo_processing_templatefit(std::vector<std::vector<float>> chnlvector);
  std::vector<std::vector<float>> calo_processing_fast(std::vector<std::vector<float>> chnlvector);
  std::vector<std::vector<float>> calo_processing_nyquist(std::vector<std::vector<float>> chnlvector);

  void initialize_processing(const std::string &templatefile);

 private:
  void FastMax(float x0, float x1, float x2, float y0, float y1, float y2, float &xmax, float &ymax);
  std::vector<float> NyquistInterpolation(std::vector<float> &vec_signal_samples);
  double Dkernelodd(double x, int N);
  double Dkernel(double x, int N);

  float stablepsinc(float t, std::vector<float> &vec_signal_samples);

  float psinc(float t, std::vector<float> &vec_signal_samples);
  double template_function(double *x, double *par);

  TProfile *h_template {nullptr};
  double m_peakTimeTemp {0};
  int _nthreads{1};
  int _nzerosuppresssamples{2};
  int _nsoftwarezerosuppression{40};
//  float _stepsize{0.001};
  float m_timeLim_low{-3.0};
  float m_timeLim_high{4.0};
  float _chi2threshold {100000};
  float _chi2lowthreshold {10000};
  float _bfr_lowpedestalthreshold {1200};
  float _bfr_highpedestalthreshold {4000};
  bool _bdosoftwarezerosuppression{false};
  bool _maxsoftwarezerosuppression{false};
  bool m_setTimeLim{false};
  bool _dobitfliprecovery {false};

  std::string m_template_input_file;
  std::string url_template;
  std::string url_onnx;
  std::string m_model_name;
};
#endif
