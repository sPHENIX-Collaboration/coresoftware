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

  void set_model_type(int modelno)
  {
    m_modeltype = modelno;
    return;
  }

  void set_template_file(std::string template_input_file)
  {
    m_template_input_file = template_input_file;
    return;
  }

  int get_model_type()
  {
    return m_modeltype;
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


  std::string m_template_input_file = std::string(getenv("CALIBRATIONROOT")) + "/WaveformProcessing/templates/testbeam_cemc_template.root";

  std::string m_model_name = std::string(getenv("CALIBRATIONROOT")) + "/WaveformProcessing/models/testbeamtrained_cemc.onnx";

};
#endif
