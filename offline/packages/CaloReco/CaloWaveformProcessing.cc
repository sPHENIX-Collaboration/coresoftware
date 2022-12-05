#include "CaloWaveformProcessing.h"

#include <ffamodules/XploadInterface.h>

#include <phool/onnxlib.h>

#include <TF1.h>
#include <TFile.h>
#include <TProfile.h>

#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Fit/Fitter.h>
#include <Fit/UnBinData.h>
#include <HFitInterface.h>
#include <Math/WrappedMultiTF1.h>
#include <Math/WrappedTF1.h>
#include <ROOT/TThreadExecutor.hxx>
#include <ROOT/TThreadedObject.hxx>

#include <iostream>
#include <pthread.h>
#include <string>


//Define some items that must be defined globally in the .cc file
TProfile *CaloWaveformProcessing::h_template = nullptr;
Ort::Session *onnxmodule;

double CaloWaveformProcessing::template_function(double *x, double *par)
{
  Double_t v1 = par[0] * h_template->Interpolate(x[0] - par[1]) + par[2];
  return v1;
}

void CaloWaveformProcessing::initialize_processing()
{
    char *calibrationsroot =  getenv("CALIBRATIONROOT");
    assert(calibrationsroot);
  if (m_processingtype == CaloWaveformProcessing::TEMPLATE)
  {
    std::string calibrations_repo_template = std::string(calibrationsroot) + "/WaveformProcessing/templates/" + m_template_input_file;
    url_template = XploadInterface::instance()->getUrl(m_template_input_file, calibrations_repo_template);
    TFile *fin = TFile::Open(url_template.c_str());
    assert(fin);
    assert(fin->IsOpen());
    h_template = static_cast<TProfile *>(fin->Get("waveform_template"));
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
    fitresults = CaloWaveformProcessing::calo_processing_templatefit(waveformvector);
  }
  if (m_processingtype == CaloWaveformProcessing::ONNX)
  {
    fitresults = CaloWaveformProcessing::calo_processing_ONNX(waveformvector);
  }
  return fitresults;
}

std::vector<std::vector<float>> CaloWaveformProcessing::calo_processing_templatefit(std::vector<std::vector<float>> chnlvector)
{
  ROOT::TThreadExecutor t(_nthreads);
  auto func = [&](std::vector<float> &v) {
    int size1 = v.size() - 1;
    auto h = new TH1F(Form("h_%d", (int) round(v.at(size1))), "", size1, 0, size1);
    float maxheight = 0;
    int maxbin = 0;
    for (int i = 0; i < size1; i++)
    {
      h->SetBinContent(i + 1, v.at(i));
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
    auto f = new TF1(Form("f_%d", (int) round(v.at(size1))), template_function, 0, 31, 3);
    ROOT::Math::WrappedMultiTF1 *fitFunction = new ROOT::Math::WrappedMultiTF1(*f, 3);
    ROOT::Fit::BinData data(v.size() - 1, 1);
    ROOT::Fit::FillData(data, h);
    ROOT::Fit::Chi2Function *EPChi2 = new ROOT::Fit::Chi2Function(data, *fitFunction);
    ROOT::Fit::Fitter *fitter = new ROOT::Fit::Fitter();
    fitter->Config().MinimizerOptions().SetMinimizerType("GSLMultiFit");
    double params[] = {static_cast<double>(maxheight), static_cast<double>(maxbin - 5), static_cast<double>(pedestal)};
    fitter->Config().SetParamsSettings(3, params);
    fitter->FitFCN(*EPChi2, nullptr, data.Size(), true);
    for (int i = 0; i < 3; i++)
    {
      v.push_back(f->GetParameter(i));
    }
    h->Delete();
    f->Delete();
    delete fitFunction;
    delete fitter;
    delete EPChi2;
  };

  t.Foreach(func, chnlvector);

  int size3 = chnlvector.size();
  std::vector<std::vector<float>> fit_params;
  std::vector<float> fit_params_tmp;
  for (int i = 0; i < size3; i++)
  {
    std::vector<float> tv = chnlvector.at(i);
    int size2 = tv.size();
    for (int q = 3; q > 0; q--)
    {
      fit_params_tmp.push_back(tv.at(size2 - q));
    }
    fit_params.push_back(fit_params_tmp);
    fit_params_tmp.clear();
  }
  chnlvector.clear();
  return fit_params;
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
