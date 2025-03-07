#include "CaloWaveformProcessing.h"
#include "CaloWaveformFitting.h"

#include <ffamodules/CDBInterface.h>

#include <phool/onnxlib.h>

#include <algorithm>  // for max
#include <cassert>
#include <cstdlib>  // for getenv
#include <iostream>
#include <memory>  // for allocator_traits<>::value_type
#include <string>

namespace
{
  Ort::Session *onnxmodule;
}

CaloWaveformProcessing::~CaloWaveformProcessing()
{
  delete m_Fitter;
}

void CaloWaveformProcessing::initialize_processing()
{
  char *calibrationsroot = getenv("CALIBRATIONROOT");
  assert(calibrationsroot);
  if (m_processingtype == CaloWaveformProcessing::TEMPLATE || m_processingtype == CaloWaveformProcessing::TEMPLATE_NOSAT)
  {
    std::string calibrations_repo_template = std::string(calibrationsroot) + "/WaveformProcessing/templates/" + m_template_input_file;
    url_template = CDBInterface::instance()->getUrl(m_template_name, calibrations_repo_template);
    m_Fitter = new CaloWaveformFitting();
    m_Fitter->initialize_processing(url_template);
    if(m_processingtype == CaloWaveformProcessing::TEMPLATE_NOSAT)
    {
      m_Fitter->set_handleSaturation(false);
    }
    m_Fitter->set_nthreads(get_nthreads());
    if (m_setTimeLim)
    {
      m_Fitter->set_timeFitLim(m_timeLim_low, m_timeLim_high);
    }

    if (_bdosoftwarezerosuppression)
    {
      m_Fitter->set_softwarezerosuppression(_bdosoftwarezerosuppression, _nsoftwarezerosuppression);
    }
    if (_dobitfliprecovery)
    {
      m_Fitter->set_bitFlipRecovery(_dobitfliprecovery);
    }
  }
  else if (m_processingtype == CaloWaveformProcessing::ONNX)
  {
    std::string calibrations_repo_model = std::string(calibrationsroot) + "/WaveformProcessing/models/" + m_model_name;
    url_onnx = CDBInterface::instance()->getUrl(m_model_name, calibrations_repo_model);
    onnxmodule = onnxSession(url_onnx);
  }
  else if (m_processingtype == CaloWaveformProcessing::NYQUIST)
  {
    std::string calibrations_repo_template = std::string(calibrationsroot) + "/WaveformProcessing/templates/" + m_template_input_file;
    url_template = CDBInterface::instance()->getUrl(m_template_name, calibrations_repo_template);
    m_Fitter = new CaloWaveformFitting();
    m_Fitter->initialize_processing(url_template);
  }
}

std::vector<std::vector<float>> CaloWaveformProcessing::process_waveform(std::vector<std::vector<float>> waveformvector)
{
  unsigned int size1 = waveformvector.size();
  std::vector<std::vector<float>> fitresults;
  if (m_processingtype == CaloWaveformProcessing::TEMPLATE)
  {
    for (unsigned int i = 0; i < size1; i++)
    {
      waveformvector.at(i).push_back((float) i);
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
  if (m_processingtype == CaloWaveformProcessing::NYQUIST)
  {
    fitresults = m_Fitter->calo_processing_nyquist(waveformvector);
  }
  return fitresults;
}

std::vector<std::vector<float>> CaloWaveformProcessing::calo_processing_ONNX(const std::vector<std::vector<float>> &chnlvector)
{
  std::vector<std::vector<float>> fit_values;
  unsigned int nchnls = chnlvector.size();
  for (unsigned int m = 0; m < nchnls; m++)
  {
    const std::vector<float> &v = chnlvector.at(m);
    unsigned int nsamples = v.size() - 1;
    std::vector<float> vtmp;
    vtmp.reserve(nsamples);
    for (unsigned int k = 0; k < nsamples; k++)
    {
      vtmp.push_back((float) (v.at(k) / 1000.0));
    }
    std::vector<float> val = onnxInference(onnxmodule, vtmp, 1, 31, 3);
    unsigned int nvals = val.size();
    for (unsigned int i = 0; i < nvals; i++)
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
