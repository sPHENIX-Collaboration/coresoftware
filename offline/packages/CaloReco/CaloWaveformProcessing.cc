#include "CaloWaveformProcessing.h"
#include "CaloWaveformFitting.h"

#include <ffamodules/CDBInterface.h>

#include <phool/onnxlib.h>

#include <algorithm>  // for max
#include <cassert>
#include <cstdlib>  // for getenv
#include <iostream>
#include <limits>
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
    if (m_processingtype == CaloWaveformProcessing::TEMPLATE_NOSAT)
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
    // std::string calibrations_repo_model = m_model_name;
    // url_onnx = CDBInterface::instance()->getUrl("CEMC_ONNX", m_model_name);
    onnxmodule = onnxSession(m_model_name);
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
  if (m_processingtype == CaloWaveformProcessing::TEMPLATE || m_processingtype == CaloWaveformProcessing::TEMPLATE_NOSAT)
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
    fitresults = CaloWaveformFitting::calo_processing_fast(waveformvector);
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
  std::vector<float> val;  // single row to return
  unsigned int nchnls = chnlvector.size();
  for (unsigned int m = 0; m < nchnls; m++)
  {
    val.clear();
    const std::vector<float> &v = chnlvector.at(m);
    int size1 = v.size();
    if (size1 == _nzerosuppresssamples)
    {
      val.push_back(v.at(1) - v.at(0));
      val.push_back(std::numeric_limits<float>::quiet_NaN());
      val.push_back(v.at(0));
      if (v.at(0) != 0 && v.at(1) == 0)  // check if post-sample is 0, if so set high chi2
      {
        val.push_back(1000000);
      }
      else
      {
        val.push_back(std::numeric_limits<float>::quiet_NaN());
      }
      val.push_back(0);
      fit_values.push_back(val);
    }
    else
    {
      float maxheight = 0;
      int maxbin = 0;
      for (int i = 0; i < size1; i++)
      {
        if (v.at(i) > maxheight)
        {
          maxheight = v.at(i);
          maxbin = i;
        }
      }
      float pedestal = 1500;
      if (maxbin > 4)
      {
        pedestal = 0.5 * (v.at(maxbin - 4) + v.at(maxbin - 5));
      }
      else if (maxbin > 3)
      {
        pedestal = (v.at(maxbin - 4));
      }
      else
      {
        pedestal = 0.5 * (v.at(size1 - 3) + v.at(size1 - 2));
      }

      if ((_bdosoftwarezerosuppression && v.at(6) - v.at(0) < _nsoftwarezerosuppression) || (_maxsoftwarezerosuppression && maxheight - pedestal < _nsoftwarezerosuppression))
      {
        val.push_back(v.at(6) - v.at(0));
        val.push_back(std::numeric_limits<float>::quiet_NaN());
        val.push_back(v.at(0));
        if (v.at(0) != 0 && v.at(1) == 0)  // check if post-sample is 0, if so set high chi2
        {
          val.push_back(1000000);
        }
        else
        {
          val.push_back(std::numeric_limits<float>::quiet_NaN());
        }
        val.push_back(0);
        fit_values.push_back(val);
      }
      else
      {
        unsigned int nsamples = v.size();
        if (nsamples == 12)
        {
          // downstream onnx does not have a static input vector API,
          // so we need to make a copy
          std::vector<float> vtmp(v);
          val = onnxInference(onnxmodule, vtmp, 1, 12, 3);
          unsigned int nvals = val.size();
          for (unsigned int i = 0; i < nvals; i++)
          {
            val.at(i) = val.at(i) * m_Onnx_factor[i] + m_Onnx_offset[i];
          }
          val.push_back(2000);
          val.push_back(0);
          fit_values.push_back(val);
        }
        else
        {
          float v_diff = v[1] - v[0];
          std::vector<float> val1{v_diff, std::numeric_limits<float>::quiet_NaN(), v[1], std::numeric_limits<float>::quiet_NaN(), 0};
          fit_values.push_back(val1);
        }
      }
    }
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
