#ifndef CALORECO_CALOWAVEFORMPROCESSING_H
#define CALORECO_CALOWAVEFORMPROCESSING_H

#include <fun4all/SubsysReco.h>

#include <TProfile.h>

#include <string>

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
  };

  CaloWaveformProcessing()
    : m_processingtype(CaloWaveformProcessing::TEMPLATE)
    , m_template_input_file("CEMC_TEMPLATE")
    , m_model_name("CEMC_ONNX"){};
  ~CaloWaveformProcessing() override {}

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
  void set_model_file(const std::string &model_name)
  {
    m_model_name = model_name;
    return;
  }

  void set_nthreads(int nthreads);

  int get_nthreads();

  void set_softwarezerosuppression(bool usezerosuppression,int softwarezerosuppression)
  {
    _nsoftwarezerosuppression = softwarezerosuppression;
    _bdosoftwarezerosuppression = usezerosuppression;
  }


  std::vector<std::vector<float>> process_waveform(std::vector<std::vector<float>> waveformvector);
  std::vector<std::vector<float>> calo_processing_ONNX(std::vector<std::vector<float>> chnlvector);

  void initialize_processing();

 private:
  CaloWaveformFitting *m_Fitter = nullptr;

  CaloWaveformProcessing::process m_processingtype = CaloWaveformProcessing::NONE;
  int _nthreads = 1;
  int _nsoftwarezerosuppression = 40;
  bool _bdosoftwarezerosuppression = false;

  std::string m_template_input_file;
  std::string url_template;

  std::string url_onnx;
  std::string m_model_name;
};
#endif
