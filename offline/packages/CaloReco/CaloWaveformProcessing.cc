#include "CaloWaveformProcessing.h"
#include "CaloWaveformFitting.h"

#include <ffamodules/XploadInterface.h>

#include <phool/onnxlib.h>

#include <TSpline.h>

#include <cassert>
#include <iostream>
#include <string>


Ort::Session *onnxmodule;

void CaloWaveformProcessing::initialize_processing()
{
    char *calibrationsroot =  getenv("CALIBRATIONROOT");
    assert(calibrationsroot);
  if (m_processingtype == CaloWaveformProcessing::TEMPLATE)
  {
    std::string calibrations_repo_template = std::string(calibrationsroot) + "/WaveformProcessing/templates/" + m_template_input_file;
    url_template = XploadInterface::instance()->getUrl(m_template_input_file, calibrations_repo_template);
    m_TemplateFitter = new CaloWaveformFitting();
    m_TemplateFitter->initialize_processing(url_template);
  }
  if (m_processingtype == CaloWaveformProcessing::ONNX)
  {
    std::string calibrations_repo_model = std::string(calibrationsroot) + "/WaveformProcessing/models/" + m_model_name;
    url_onnx = XploadInterface::instance()->getUrl(m_model_name, calibrations_repo_model);
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
    fitresults = m_TemplateFitter->calo_processing_templatefit(waveformvector);
  }
  if (m_processingtype == CaloWaveformProcessing::ONNX)
  {
    fitresults = CaloWaveformProcessing::calo_processing_ONNX(waveformvector);
  }
  if (m_processingtype == CaloWaveformProcessing::FAST)
  {
    fitresults = CaloWaveformProcessing::calo_processing_fast(waveformvector);
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

void CaloWaveformProcessing::FastMax(float x0, float x1, float x2, float y0, float y1, float y2, float & xmax, float & ymax) {
  int n = 3;
  double xp[3] = {x0, x1, x2};
  double yp[3] = {y0, y1, y2};
  TSpline3 *sp = new TSpline3("", xp, yp, n, "b2e2", 0, 0);
  double X, Y, B, C, D;
  ymax = y1;
  xmax = x1;
  if (y0 > ymax) {
    ymax = y0;
    xmax = x0;
  }
  if (y2 > ymax) {
    ymax = y2;
    xmax = x2;
  }
  for (int i = 0; i <= 1; i++) {
    sp->GetCoeff(i, X, Y, B, C, D);
    if (D == 0) {

      if (C < 0) {
        //TSpline is a quadratic equation

        float root = -B / (2 * C) + X;
        if (root >= xp[i] && root <= xp[i + 1]) {
          float yvalue = sp->Eval(root);
          if (yvalue > ymax) {
            ymax = yvalue;
            xmax = root;
          }
        }
      }
    }
    else {
      //find x when derivative = 0
      float root = (-2 * C + sqrt(4 * C * C - 12 * B * D)) / (6 * D) + X;
      if (root >= xp[i] && root <= xp[i + 1]) {
        float yvalue = sp->Eval(root);
        if (yvalue > ymax) {
          ymax = yvalue;
          xmax = root;
        }
      }
      root = (-2 * C - sqrt(4 * C * C - 12 * B * D)) / (6 * D) + X;
      if (root >= xp[i] && root <= xp[i + 1]) {
        float yvalue = sp->Eval(root);
        if (yvalue > ymax) {
          ymax = yvalue;
          xmax = root;
        }
      }

    }
  }
  delete sp;
  return;
}
std::vector<std::vector<float>> CaloWaveformProcessing::calo_processing_fast(std::vector<std::vector<float>> chnlvector)
{
  std::vector<std::vector<float>> fit_values;
  int nchnls = chnlvector.size();
  for (int m = 0; m < nchnls; m++)
  {
    std::vector<float> v = chnlvector.at(m);
    int nsamples = v.size();

    double maxy = v.at(0);
    float amp = 0;
    float time = 0;
    float ped = 0;
    if (nsamples >= 3) {
      int maxx = 0;
      for (int i = 0; i < nsamples; i++) {
        if (i < 3) ped += v.at(i);
        if (v.at(i) > maxy) {
          maxy = v.at(i);
          maxx = i;
        }
      }
      ped /= 3;
      if (maxx == 0 || maxx == nsamples - 1) {
        amp = maxy;
        time = maxx;
      }
      else {
        FastMax(maxx - 1, maxx, maxx + 1, v.at(maxx - 1), v.at(maxx), v.at(maxx + 1), time, amp);

      }
    }
    amp -= ped;
    std::vector<float> val = {amp, time, ped};
    fit_values.push_back(val);
    val.clear();
  }
  return fit_values;
}
