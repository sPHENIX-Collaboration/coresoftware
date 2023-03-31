#ifndef CALORECO_CALOWAVEFORMFITTING_H
#define CALORECO_CALOWAVEFORMFITTING_H

#include <string>
#include <vector>

class TProfile;

class CaloWaveformFitting
{
 public:
  CaloWaveformFitting() = default;
  ~CaloWaveformFitting() = default;

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

  int get_nthreads()
  {
    return _nthreads;
  }

  std::vector<std::vector<float>> process_waveform(std::vector<std::vector<float>> waveformvector);
  std::vector<std::vector<float>> calo_processing_templatefit(std::vector<std::vector<float>> chnlvector);
  std::vector<std::vector<float>> calo_processing_fast(std::vector<std::vector<float>> chnlvector);

  void initialize_processing(const std::string &templatefile);

 private:
  void FastMax(float x0, float x1, float x2, float y0, float y1, float y2, float &xmax, float &ymax);
  static TProfile *h_template;
  static double template_function(double *x, double *par);
  int _nthreads = 1;
  std::string m_template_input_file;
  std::string url_template;

  std::string url_onnx;
  std::string m_model_name;
};
#endif
