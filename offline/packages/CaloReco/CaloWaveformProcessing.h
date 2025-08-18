#ifndef CALORECO_CALOWAVEFORMPROCESSING_H
#define CALORECO_CALOWAVEFORMPROCESSING_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>
#include <vector>

class CaloWaveformFitting;

class CaloWaveformProcessing : public SubsysReco
{
 public:
  enum process
  {
    NONE = 0,
    TEMPLATE = 1,
    ONNX = 2,
    FAST = 3,
    NYQUIST = 4,
    TEMPLATE_NOSAT = 5,
  };

  CaloWaveformProcessing() = default;
  ~CaloWaveformProcessing() override;

  void set_processing_type(CaloWaveformProcessing::process modelno)
  {
    m_processingtype = modelno;
    return;
  }

  CaloWaveformProcessing::process get_processing_type()
  {
    return m_processingtype;
  }

  void set_template_file(const std::string &template_input_file)
  {
    m_template_input_file = template_input_file;
    return;
  }

  void set_template_name(const std::string &template_name)
  {
    m_template_name = template_name;
    return;
  }
  void set_model_file(const std::string &model_name)
  {
    m_model_name = model_name;
    return;
  }

  void set_nthreads(int nthreads);

  int get_nthreads();

  void set_softwarezerosuppression(bool usezerosuppression, int softwarezerosuppression)
  {
    _nsoftwarezerosuppression = softwarezerosuppression;
    _bdosoftwarezerosuppression = usezerosuppression;
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

  std::vector<std::vector<float>> process_waveform(std::vector<std::vector<float>> waveformvector);
  std::vector<std::vector<float>> calo_processing_ONNX(const std::vector<std::vector<float>> &chnlvector);

  void initialize_processing();

  // onnx options
  void set_onnx_factor(const int i, const double val) { m_Onnx_factor.at(i) = val; }
  void set_onnx_offset(const int i, const double val) { m_Onnx_offset.at(i) = val; }

 private:
  CaloWaveformFitting *m_Fitter{nullptr};

  CaloWaveformProcessing::process m_processingtype{CaloWaveformProcessing::TEMPLATE};
  int _nthreads{1};
  int _nzerosuppresssamples{2};

  int _nsoftwarezerosuppression{40};
  bool _bdosoftwarezerosuppression{false};
  bool _dobitfliprecovery{false};

  std::string m_template_input_file;
  std::string url_template;
  std::string m_template_name{"NONE"};

  bool m_setTimeLim{false};
  bool _maxsoftwarezerosuppression{false};
  float m_timeLim_low{-3.0};
  float m_timeLim_high{4.0};

  std::string url_onnx;
  std::string m_model_name{"CEMC_ONNX"};
  std::array<double, 3> m_Onnx_factor{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
  std::array<double, 3> m_Onnx_offset{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
};

#endif
