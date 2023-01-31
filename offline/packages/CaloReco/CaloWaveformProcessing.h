#ifndef CALORECO_CALOWAVEFORMPROCESSING_H
#define CALORECO_CALOWAVEFORMPROCESSING_H

#include <fun4all/SubsysReco.h>

#include<TProfile.h>

#include <string>

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
    :m_processingtype(CaloWaveformProcessing::TEMPLATE)
    , m_template_input_file("CEMC_TEMPLATE")
    , m_model_name("CEMC_ONNX")
{};
  ~CaloWaveformProcessing() override  { }

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

  void set_nthreads (int nthreads)
  {
    _nthreads = nthreads;
    return;
  }

  int get_nthreads ()
  {
    return _nthreads;
  }

  std::vector<std::vector<float>>  process_waveform(std::vector<std::vector<float>> waveformvector);
  std::vector<std::vector<float>>  calo_processing_ONNX(std::vector<std::vector<float>> chnlvector);
  std::vector<std::vector<float>>  calo_processing_templatefit(std::vector<std::vector<float>> chnlvector);
  std::vector<std::vector<float>>  calo_processing_fast(std::vector<std::vector<float>> chnlvector);


  void initialize_processing(); 

 private:

  static TProfile* h_template; 
  static double template_function(double *x, double *par);

  CaloWaveformProcessing::process m_processingtype = CaloWaveformProcessing::NONE; 
  int _nthreads = 1;

  std::string m_template_input_file;
  std::string url_template;

  std::string url_onnx;
  std::string m_model_name;
  void FastMax(float x0, float x1, float x2, float y0, float y1, float y2, float & xmax, float & ymax);
};
#endif
