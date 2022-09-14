#ifndef CALORECO_CALOWAVEFORMPROCESSING_H
#define CALORECO_CALOWAVEFORMPROCESSING_H

#include <fun4all/SubsysReco.h>
#include <string>
#include<TProfile.h>


class CaloWaveformProcessing : public SubsysReco
{
 public:
  CaloWaveformProcessing(){};
  ~CaloWaveformProcessing() override  { }

  void set_processing_type(int modelno)
  {
    m_modeltype = modelno;
    return;
  }

  void set_template_file(std::string template_input_file)
  {
    m_template_input_file = template_input_file;
    return;
  }
  void set_model_file(std::string model_name)
  {
    m_model_name = model_name;
    return;
  }

  int get_model_type()
  {
    return m_modeltype;
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

  std::vector<std::vector<float>> process_waveform(std::vector<std::vector<float>> waveformvector);
  std::vector<std::vector<float>>  calo_processing_ONNX(std::vector<std::vector<float>> chnlvector);
  std::vector<std::vector<float>>  calo_processing_templatefit(std::vector<std::vector<float>> chnlvector);


  void initialize_processing(); 

 private:

  int m_modeltype;
  int _nthreads = 1;
  static TProfile* h_template; 
  static double template_function(double *x, double *par);
  


  std::string m_template_input_file = "testbeam_cemc_template.root";
  std::string url_template;



  std::string url_onnx;
  std::string m_model_name = "testbeamtrained_cemc.onnx";

};
#endif
