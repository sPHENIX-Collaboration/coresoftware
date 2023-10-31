#include "CaloWaveformProcessing.h"
#include "CaloWaveformFitting.h"

#include <ffamodules/CDBInterface.h>

#include <phool/onnxlib.h>

#include <algorithm>                  // for max
#include <cassert>
#include <cstdlib>                   // for getenv
#include <iostream>
#include <memory>                     // for allocator_traits<>::value_type
#include <string>

Ort::Session *onnxmodule;

CaloWaveformProcessing::~CaloWaveformProcessing()
{
  delete m_Fitter;
}

void CaloWaveformProcessing::initialize_processing()
{
  char *calibrationsroot = getenv("CALIBRATIONROOT");
  assert(calibrationsroot);
  if (m_processingtype == CaloWaveformProcessing::TEMPLATE)
  {
    std::string calibrations_repo_template = std::string(calibrationsroot) + "/WaveformProcessing/templates/" + m_template_input_file;
    url_template = CDBInterface::instance()->getUrl(m_template_name, calibrations_repo_template);
    m_Fitter = new CaloWaveformFitting();
    m_Fitter->initialize_processing(url_template);
    m_Fitter->set_nthreads(get_nthreads());

    if (_bdosoftwarezerosuppression == true)
      {
	m_Fitter->set_softwarezerosuppression(_bdosoftwarezerosuppression,_nsoftwarezerosuppression);
      }
  }
  else if (m_processingtype == CaloWaveformProcessing::ONNX)
  {
    std::string calibrations_repo_model = std::string(calibrationsroot) + "/WaveformProcessing/models/" + m_model_name;
    url_onnx = CDBInterface::instance()->getUrl(m_model_name, calibrations_repo_model);
    onnxmodule = onnxSession(url_onnx);
  }
}

std::vector<std::vector<float>> CaloWaveformProcessing::process_waveform(std::vector<std::vector<float>> waveformvector)
{
  int size1 = waveformvector.size();
  std::vector<std::vector<float>> fitresults;
  if (m_processingtype == CaloWaveformProcessing::TEMPLATE)
  {
    for (int i = 0; i < size1; i++)
    {
      waveformvector.at(i).push_back(i);
    }
    fitresults = m_Fitter->calo_processing_templatefit(waveformvector);
  }
  if (m_processingtype == CaloWaveformProcessing::ONNX)
  {
    fitresults = CaloWaveformProcessing::calo_processing_ONNX(waveformvector);
  }
  if (m_processingtype == CaloWaveformProcessing::FAST)
  {
    fitresults = m_Fitter->calo_processing_fast(waveformvector);
  }
  return fitresults;
}

std::vector<std::vector<float>> CaloWaveformProcessing::calo_processing_ONNX(std::vector<std::vector<float>> chnlvector)
{
  std::vector<std::vector<float>> fit_values;
  int nchnls = chnlvector.size();
  for (int m = 0; m < nchnls; m++)
  {
    std::vector<float> v = chnlvector.at(m);
    int nsamples = v.size() - 1;
    std::vector<float> vtmp;
    vtmp.reserve(nsamples);
    for (int k = 0; k < nsamples; k++)
    {
      vtmp.push_back(v.at(k) / 1000.0);
    }
    std::vector<float> val = onnxInference(onnxmodule, vtmp, 1, 31, 3);
    int nvals = val.size();
    for (int i = 0; i < nvals; i++)
    {
      if (i == 0 || i == 2)
      {
        val.at(i) = val.at(i) * 1000;
      }
    }
    fit_values.push_back(val);
    val.clear();
  }
  return fit_values;
}

int CaloWaveformProcessing::get_nthreads()
{
  if (m_Fitter)
  {
    return m_Fitter->get_nthreads();
  }
  return _nthreads;
}
void CaloWaveformProcessing::set_nthreads(int nthreads)
{
  _nthreads = nthreads;
  if (m_Fitter)
  {
    m_Fitter->set_nthreads(nthreads);
  }
  return;
}
