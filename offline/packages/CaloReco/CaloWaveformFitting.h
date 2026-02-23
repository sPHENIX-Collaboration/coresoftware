#ifndef CALORECO_CALOWAVEFORMFITTING_H
#define CALORECO_CALOWAVEFORMFITTING_H

#include <string>
#include <vector>

class TProfile;

class CaloWaveformFitting
{
 public:
  enum FuncFitType
  {
    POWERLAWEXP = 0,
    POWERLAWDOUBLEEXP = 1,
    FERMIEXP = 2,
  };

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
  void set_timeFitLim(float low, float high)
  {
    m_setTimeLim = true;
    m_timeLim_low = low;
    m_timeLim_high = high;
    return;
  }

  void set_bitFlipRecovery(bool dobitfliprecovery)
  {
    _dobitfliprecovery = dobitfliprecovery;
  }

  void set_handleSaturation(bool handleSaturation = true)
  {
    _handleSaturation = handleSaturation;
  }

  std::vector<std::vector<float>> process_waveform(std::vector<std::vector<float>> waveformvector);
  std::vector<std::vector<float>> calo_processing_templatefit(std::vector<std::vector<float>> chnlvector);
  static std::vector<std::vector<float>> calo_processing_fast(const std::vector<std::vector<float>> &chnlvector);
  std::vector<std::vector<float>> calo_processing_nyquist(const std::vector<std::vector<float>> &chnlvector);
  std::vector<std::vector<float>> calo_processing_funcfit(const std::vector<std::vector<float>> &chnlvector);

  void initialize_processing(const std::string &templatefile);

  // Power-law fit function: amplitude * (x-t0)^power * exp(-(x-t0)*decay) + pedestal
  static double SignalShape_PowerLawExp(double *x, double *par);
  // Double exponential power-law fit function
  static double SignalShape_PowerLawDoubleExp(double *x, double *par);
  static double SignalShape_FermiExp(double *x, double *par);

  void set_funcfit_type(FuncFitType type)
  {
    m_funcfit_type = type;
  }

  void set_powerlaw_params(double power, double decay)
  {
    m_powerlaw_power = power;
    m_powerlaw_decay = decay;
  }

  void set_doubleexp_params(double power, double peaktime1, double peaktime2, double ratio)
  {
    m_doubleexp_power = power;
    m_doubleexp_peaktime1 = peaktime1;
    m_doubleexp_peaktime2 = peaktime2;
    m_doubleexp_ratio = ratio;
  }

 private:
  static void FastMax(float x0, float x1, float x2, float y0, float y1, float y2, float &xmax, float &ymax);
  std::vector<float> NyquistInterpolation(std::vector<float> &vec_signal_samples);
  static double Dkernelodd(double x, int N);
  static double Dkernel(double x, int N);

  static float stablepsinc(float t, std::vector<float> &vec_signal_samples);

  static float psinc(float t, std::vector<float> &vec_signal_samples);
  double template_function(double *x, double *par);

  TProfile *h_template{nullptr};
  double m_peakTimeTemp{0};
  int _nthreads{1};
  int _nzerosuppresssamples{2};
  int _nsoftwarezerosuppression{40};
  //  float _stepsize{0.001};
  float m_timeLim_low{-3.0};
  float m_timeLim_high{4.0};
  float _chi2threshold{100000};
  float _chi2lowthreshold{10000};
  float _bfr_lowpedestalthreshold{1200};
  float _bfr_highpedestalthreshold{4000};
  bool _bdosoftwarezerosuppression{false};
  bool _maxsoftwarezerosuppression{false};
  bool m_setTimeLim{false};
  bool _dobitfliprecovery{false};
  bool _handleSaturation{true};

  std::string m_template_input_file;
  std::string url_template;
  std::string url_onnx;
  std::string m_model_name;

  // Functional fit type selector
  FuncFitType m_funcfit_type{POWERLAWDOUBLEEXP};

  // Power-law fit parameters
  double m_powerlaw_power{4.0};
  double m_powerlaw_decay{1.5};

  // Double exponential fit parameters
  double m_doubleexp_power{2.0};
  double m_doubleexp_peaktime1{5.0};
  double m_doubleexp_peaktime2{5.0};
  double m_doubleexp_ratio{0.3};
};
#endif
