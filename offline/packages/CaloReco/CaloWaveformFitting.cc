#include "CaloWaveformFitting.h"

#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TSpline.h>

#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Fit/Fitter.h>
#include <Fit/UnBinData.h>
#include <HFitInterface.h>
#include <Math/WrappedMultiTF1.h>
#include <Math/WrappedTF1.h>
#include <ROOT/TThreadExecutor.hxx>
#include <ROOT/TThreadedObject.hxx>

#include <pthread.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>

static ROOT::TThreadExecutor *t = new ROOT::TThreadExecutor(1);  // NOLINT(misc-use-anonymous-namespace)
double CaloWaveformFitting::template_function(double *x, double *par)
{
  Double_t v1 = (par[0] * h_template->Interpolate(x[0] - par[1])) + par[2];
  return v1;
}

CaloWaveformFitting::~CaloWaveformFitting()
{
  delete h_template;
}

void CaloWaveformFitting::initialize_processing(const std::string &templatefile)
{
  TFile *fin = TFile::Open(templatefile.c_str());
  assert(fin);
  assert(fin->IsOpen());
  h_template = dynamic_cast<TProfile *>(fin->Get("waveform_template"));
  h_template->SetDirectory(nullptr);
  fin->Close();
  delete fin;
  m_peakTimeTemp = h_template->GetBinCenter(h_template->GetMaximumBin());
  t = new ROOT::TThreadExecutor(_nthreads);
}

std::vector<std::vector<float>> CaloWaveformFitting::process_waveform(std::vector<std::vector<float>> waveformvector)
{
  int size1 = waveformvector.size();
  std::vector<std::vector<float>> fitresults;
  for (int i = 0; i < size1; i++)
  {
    waveformvector.at(i).push_back(i);
  }
  fitresults = calo_processing_templatefit(waveformvector);
  return fitresults;
}

std::vector<std::vector<float>> CaloWaveformFitting::calo_processing_templatefit(std::vector<std::vector<float>> chnlvector)
{
  auto func = [&](std::vector<float> &v)
  {
    int size1 = v.size() - 1;
    if (size1 == _nzerosuppresssamples)
    {
      v.push_back(v.at(1) - v.at(0));                        // returns peak sample - pedestal sample
      v.push_back(std::numeric_limits<float>::quiet_NaN());  // set time to qnan for ZS
      v.push_back(v.at(0));
      if (v.at(0) != 0 && v.at(1) == 0)  // check if post-sample is 0, if so set high chi2
      {
        v.push_back(1000000);
      }
      else
      {
        v.push_back(std::numeric_limits<float>::quiet_NaN());
      }
      v.push_back(0);
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
        v.push_back(v.at(6) - v.at(0));
        v.push_back(std::numeric_limits<float>::quiet_NaN());
        v.push_back(v.at(0));
        if (v.at(0) != 0 && v.at(1) == 0)  // check if post-sample is 0, if so set high chi2
        {
          v.push_back(1000000);
        }
        else
        {
          v.push_back(std::numeric_limits<float>::quiet_NaN());
        }
        v.push_back(0);
      }
      else
      {
        auto *h = new TH1F(std::string("h_" + std::to_string((int) round(v.at(size1)))).c_str(), "", size1, -0.5, size1 - 0.5);

        int ndata = 0;
        for (int i = 0; i < size1; ++i)
        {
          if ((v.at(i) == 16383) && _handleSaturation)
          {
            continue;
          }

          h->SetBinContent(i + 1, v.at(i));
          h->SetBinError(i + 1, 1);
          ndata++;
        }
        // if too many are saturated don't do the saturation recovery need enough ndf
        if (ndata < (size1 - 4))
        {
          ndata = size1;
          for (int i = 0; i < size1; ++i)
          {
            h->SetBinContent(i + 1, v.at(i));
            h->SetBinError(i + 1, 1);
          }
        }

        auto *f = new TF1(std::string("f_" + std::to_string((int) round(v.at(size1)))).c_str(), this, &CaloWaveformFitting::template_function, 0, 31, 3, "CaloWaveformFitting", "template_function");
        ROOT::Math::WrappedMultiTF1 *fitFunction = new ROOT::Math::WrappedMultiTF1(*f, 3);
        ROOT::Fit::BinData data(v.size() - 1, 1);
        ROOT::Fit::FillData(data, h);
        ROOT::Fit::Chi2Function *EPChi2 = new ROOT::Fit::Chi2Function(data, *fitFunction);
        ROOT::Fit::Fitter *fitter = new ROOT::Fit::Fitter();
        fitter->Config().MinimizerOptions().SetMinimizerType("GSLMultiFit");
        fitter->Config().MinimizerOptions().SetPrintLevel(-1);
        double params[] = {static_cast<double>(maxheight - pedestal), static_cast<double>(maxbin - m_peakTimeTemp), static_cast<double>(pedestal)};
        // double params[] = {static_cast<double>(maxheight - pedestal), 0, static_cast<double>(pedestal)};
        fitter->Config().SetParamsSettings(3, params);
        fitter->Config().ParSettings(1).SetLimits(-1 * m_peakTimeTemp, size1 - m_peakTimeTemp);  // set lim on time par
        if (m_setTimeLim)
        {
          fitter->Config().ParSettings(1).SetLimits(m_timeLim_low, m_timeLim_high);
        }
        fitter->FitFCN(*EPChi2, nullptr, data.Size(), true);
        ROOT::Fit::FitResult fitres = fitter->Result();
        // get the result status
        /*
        bool validfit = fitres.IsValid();
        if(!validfit)
        {
          std::cout<<"invalid fit"<<std::endl;
          for (int i = 0; i < size1; ++i)
        {
          std::cout<<v.at(i)<<std::endl;
        }
        }
        */
        double chi2min = fitres.MinFcnValue();
        // chi2min /= size1 - 3;  // divide by the number of dof
        chi2min /= ndata - 3;  // divide by the number of dof
        if (chi2min > _chi2threshold && (f->GetParameter(2) < _bfr_highpedestalthreshold || pedestal < _bfr_highpedestalthreshold) && (f->GetParameter(2) > _bfr_lowpedestalthreshold || pedestal > _bfr_lowpedestalthreshold) && _dobitfliprecovery)
        {
          std::vector<float> rv;  // temporary recovered waveform
          rv.reserve(size1);
          for (int i = 0; i < size1; i++)
          {
            rv.push_back(v.at(i));
          }
          unsigned int bits[3] = {8192, 4096, 2048};
          for (auto bit : bits)
          {
            for (int i = 0; i < size1; i++)
            {
              if (((unsigned int) rv.at(i) & bit) && ((unsigned int) rv.at(i) % bit > _bfr_lowpedestalthreshold))
              {
                rv.at(i) = rv.at(i) - bit;
              }
            }
          }
          for (int i = 0; i < size1; i++)
          {
            h->SetBinContent(i + 1, rv.at(i));
            h->SetBinError(i + 1, 1);
          }

          maxheight = 0;
          maxbin = 0;
          for (int i = 0; i < size1; i++)
          {
            if (rv.at(i) > maxheight)
            {
              maxheight = rv.at(i);
              maxbin = i;
            }
          }
          if (maxbin > 4)
          {
            pedestal = 0.5 * (rv.at(maxbin - 4) + rv.at(maxbin - 5));
          }
          else if (maxbin > 3)
          {
            pedestal = (rv.at(maxbin - 4));
          }
          else
          {
            pedestal = 0.5 * (rv.at(size1 - 3) + rv.at(size1 - 2));
          }

          auto *recover_f = new TF1(std::string("recover_f_" + std::to_string((int) round(v.at(size1)))).c_str(), this, &CaloWaveformFitting::template_function, 0, 31, 3, "CaloWaveformFitting", "template_function");
          ROOT::Math::WrappedMultiTF1 *recoverFitFunction = new ROOT::Math::WrappedMultiTF1(*recover_f, 3);
          ROOT::Fit::BinData recoverData(rv.size() - 1, 1);
          ROOT::Fit::FillData(recoverData, h);
          ROOT::Fit::Chi2Function *recoverEPChi2 = new ROOT::Fit::Chi2Function(recoverData, *recoverFitFunction);
          ROOT::Fit::Fitter *recoverFitter = new ROOT::Fit::Fitter();
          recoverFitter->Config().MinimizerOptions().SetMinimizerType("GSLMultiFit");
          double recover_params[] = {static_cast<double>(maxheight - pedestal), 0, static_cast<double>(pedestal)};
          recoverFitter->Config().SetParamsSettings(3, recover_params);
          recoverFitter->Config().ParSettings(1).SetLimits(-1 * m_peakTimeTemp, size1 - m_peakTimeTemp);  // set lim on time par
          recoverFitter->FitFCN(*recoverEPChi2, nullptr, recoverData.Size(), true);
          ROOT::Fit::FitResult recover_fitres = recoverFitter->Result();
          double recover_chi2min = recover_fitres.MinFcnValue();
          recover_chi2min /= size1 - 3;  // divide by the number of dof
          if (recover_chi2min < _chi2lowthreshold && recover_f->GetParameter(2) < _bfr_highpedestalthreshold && recover_f->GetParameter(2) > _bfr_lowpedestalthreshold)
          {
            for (int i = 0; i < size1; i++)
            {
              v.at(i) = rv.at(i);
            }
            for (int i = 0; i < 3; i++)
            {
              v.push_back(recover_f->GetParameter(i));
            }
            v.push_back(recover_chi2min);
            v.push_back(1);
          }
          else
          {
            for (int i = 0; i < 3; i++)
            {
              v.push_back(f->GetParameter(i));
            }
            v.push_back(chi2min);
            v.push_back(0);
          }
          recover_f->Delete();
          delete recoverFitFunction;
          delete recoverFitter;
          delete recoverEPChi2;
        }
        else
        {
          for (int i = 0; i < 3; i++)
          {
            v.push_back(f->GetParameter(i));
          }
          v.push_back(chi2min);
          v.push_back(0);
        }
        h->Delete();
        f->Delete();
        delete fitFunction;
        delete fitter;
        delete EPChi2;
      }
    }
  };

  t->Foreach(func, chnlvector);
  int size3 = chnlvector.size();
  std::vector<std::vector<float>> fit_params;
  std::vector<float> fit_params_tmp;
  for (int i = 0; i < size3; i++)
  {
    const std::vector<float> &tv = chnlvector.at(i);
    int size2 = tv.size();
    for (int q = 5; q > 0; q--)
    {
      fit_params_tmp.push_back(tv.at(size2 - q));
    }
    fit_params.push_back(fit_params_tmp);
    fit_params_tmp.clear();
  }
  chnlvector.clear();
  return fit_params;
}

void CaloWaveformFitting::FastMax(float x0, float x1, float x2, float y0, float y1, float y2, float &xmax, float &ymax)
{
  int n = 3;
  double xp[3] = {x0, x1, x2};
  double yp[3] = {y0, y1, y2};
  TSpline3 *sp = new TSpline3("", xp, yp, n, "b2e2", 0, 0);
  double X;
  double Y;
  double B;
  double C;
  double D;
  ymax = y1;
  xmax = x1;
  if (y0 > ymax)
  {
    ymax = y0;
    xmax = x0;
  }
  if (y2 > ymax)
  {
    ymax = y2;
    xmax = x2;
  }
  for (int i = 0; i <= 1; i++)
  {
    sp->GetCoeff(i, X, Y, B, C, D);
    if (D == 0)
    {
      if (C < 0)
      {
        // TSpline is a quadratic equation

        float root = (-B / (2 * C)) + X;
        if (root >= xp[i] && root <= xp[i + 1])
        {
          float yvalue = sp->Eval(root);
          if (yvalue > ymax)
          {
            ymax = yvalue;
            xmax = root;
          }
        }
      }
    }
    else
    {
      // find x when derivative = 0
      float root = ((-2 * C + sqrt((4 * C * C) - (12 * B * D))) / (6 * D)) + X;
      if (root >= xp[i] && root <= xp[i + 1])
      {
        float yvalue = sp->Eval(root);
        if (yvalue > ymax)
        {
          ymax = yvalue;
          xmax = root;
        }
      }
      root = (-2 * C - sqrt((4 * C * C) - (12 * B * D))) / (6 * D) + X;
      if (root >= xp[i] && root <= xp[i + 1])
      {
        float yvalue = sp->Eval(root);
        if (yvalue > ymax)
        {
          ymax = yvalue;
          xmax = root;
        }
      }
    }
  }
  delete sp;
  return;
}
std::vector<std::vector<float>> CaloWaveformFitting::calo_processing_fast(const std::vector<std::vector<float>> &chnlvector)
{
  std::vector<std::vector<float>> fit_values;
  int nchnls = chnlvector.size();
  for (int m = 0; m < nchnls; m++)
  {
    const std::vector<float> &v = chnlvector.at(m);
    int nsamples = v.size();

    double maxy = v.at(0);
    float amp = 0;
    float time = 0;
    float ped = 0;
    float chi2 = std::numeric_limits<float>::quiet_NaN();
    if (nsamples == 2)
    {
      amp = v.at(1);
      time = std::numeric_limits<float>::quiet_NaN();
      ped = v.at(0);
      if (v.at(0) != 0 && v.at(1) == 0)  // check if post-sample is 0, if so set high chi2
      {
        chi2 = 1000000;
      }
    }
    else if (nsamples >= 3)
    {
      int maxx = 0;
      for (int i = 0; i < nsamples; i++)
      {
        if (i < 3)
        {
          ped += v.at(i);
        }
        if (v.at(i) > maxy)
        {
          maxy = v.at(i);
          maxx = i;
        }
      }
      ped /= 3;
      // if maxx <=5 nsample >=10 use the last two sample for pedestal(for HCal TP)
      if (maxx <= 5 && nsamples >= 10)
      {
        ped = 0.5 * (v.at(nsamples - 2) + v.at(nsamples - 1));
      }
      if (maxx == 0 || maxx == nsamples - 1)
      {
        amp = maxy;
        time = maxx;
      }
      else
      {
        FastMax(maxx - 1, maxx, maxx + 1, v.at(maxx - 1), v.at(maxx), v.at(maxx + 1), time, amp);
      }
    }
    amp -= ped;
    std::vector<float> val = {amp, time, ped, chi2, 0};
    fit_values.push_back(val);
    val.clear();
  }
  return fit_values;
}

std::vector<std::vector<float>> CaloWaveformFitting::calo_processing_nyquist(const std::vector<std::vector<float>> &chnlvector)
{
  std::vector<std::vector<float>> fit_values;
  int nchnls = chnlvector.size();
  for (int m = 0; m < nchnls; m++)
  {
    std::vector<float> v = chnlvector.at(m);
    int nsamples = (int) v.size();

    if (nsamples == 2)
    {
      float chi2 = std::numeric_limits<float>::quiet_NaN();
      if (v.at(0) != 0 && v.at(1) == 0)  // check if post-sample is 0, if so set high chi2
      {
        chi2 = 1000000;
      }
      fit_values.push_back({v.at(1) - v.at(0), std::numeric_limits<float>::quiet_NaN(), v.at(0), chi2, 0});
      continue;
    }

    std::vector<float> result = NyquistInterpolation(v);
    fit_values.push_back(result);
  }
  return fit_values;
}
// mabye I can find a way to make it thread safe
std::vector<float> CaloWaveformFitting::NyquistInterpolation(std::vector<float> &vec_signal_samples)
{
  // int N = (int) vec_signal_samples.size();
  auto max_elem_iter = std::max_element(vec_signal_samples.begin(), vec_signal_samples.end());
  int maxx = std::distance(vec_signal_samples.begin(), max_elem_iter);
  float max = *max_elem_iter;

  float maxpos = maxx;
  float steplength = 0.5;

  while (steplength > 0.001)
  {
    // use 1.5 instead of 1 to avoid the floating point error...
    float starttime = maxpos - (1 * steplength);
    float endtime = maxpos + (1.5 * steplength);

    for (float i = starttime; i < endtime; i += steplength)
    {
      float yval = max;
      if (i != maxpos)
      {
        yval = psinc(i, vec_signal_samples);
      }
      if (yval > max)
      {
        max = yval;
        maxpos = i;
      }
    }
    steplength /= 2;
  }

  float pedestal = 0;

  if (maxpos > 5)
  {
    for (int i = 0; i < 3; i++)
    {
      pedestal += vec_signal_samples[i];
    }
    pedestal = pedestal / 3;
  }
  else if (maxpos > 4)
  {
    pedestal = (vec_signal_samples[0] + vec_signal_samples[1]) / 2;
  }
  // need more consideration for what is the most effieicnt
  else
  {
    pedestal = max;
    for (float i = maxpos - 5; i < maxpos; i += 0.1)
    {
      float yval = psinc(i, vec_signal_samples);
      pedestal = std::min(yval, pedestal);
    }
  }
  // calculate chi2 using the tempalte
  float chi2 = 0;
  double par[3] = {max - pedestal, maxpos - m_peakTimeTemp, pedestal};
  for (int i = 0; i < (int) vec_signal_samples.size(); i++)
  {
    double xval[1] = {(double) i};
    float diff = vec_signal_samples[i] - template_function(xval, par);
    chi2 += diff * diff;
  }
  std::vector<float> val = {max - pedestal, maxpos, pedestal, chi2, 0};
  return val;
}

// for odd N
double CaloWaveformFitting::Dkernelodd(double x, int N)
{
  double sum = 0;
  for (int k = 0; k < (N + 1) / 2; k++)
  {
    sum += 2 * std::cos(2 * M_PI * k * x / N);
  }
  sum -= 1;
  sum = sum / N;
  return sum;
}

// for even N
double CaloWaveformFitting::Dkernel(double x, int N)
{
  double sum = 0;
  for (int k = 0; k < N / 2; k++)
  {
    sum += 2 * std::cos(2 * M_PI * k * x / N);
  }
  sum -= 1;
  sum += std::cos(M_PI * x);
  sum = sum / N;
  return sum;
}

float CaloWaveformFitting::stablepsinc(float time, std::vector<float> &vec_signal_samples)
{
  int N = (int) vec_signal_samples.size();
  float sum = 0;
  if (N % 2 == 0)
  {
    for (int n = 0; n < N; n++)
    {
      sum += vec_signal_samples[n] * Dkernel(time - n, N);
    }
  }
  else
  {
    for (int n = 0; n < N; n++)
    {
      sum += vec_signal_samples[n] * Dkernelodd(time - n, N);
    }
  }
  return sum;
}

float CaloWaveformFitting::psinc(float time, std::vector<float> &vec_signal_samples)
{
  int N = (int) vec_signal_samples.size();

  if (std::abs(std::round(time) - time) < 1e-6)
  {
    if (time < 0 || time >= N)
    {
      return stablepsinc(time, vec_signal_samples);
    }

    return vec_signal_samples.at(std::round(time));
  }

  float sum = 0;
  if (N % 2 == 0)
  {
    for (int n = 0; n < N; n++)
    {
      double piu = M_PI * (time - n);
      double piuN = piu / N;
      sum += vec_signal_samples[n] * std::sin(piu) / (std::tan(piuN)) / N;
    }
  }
  else
  {
    for (int n = 0; n < N; n++)
    {
      double piu = M_PI * (time - n);
      double piuN = piu / N;
      sum += vec_signal_samples[n] * std::sin(piu) / (std::sin(piuN)) / N;
    }
  }

  return sum;
}

double CaloWaveformFitting::SignalShape_PowerLawExp(double *x, double *par)
{
  // par[0]: Amplitude
  // par[1]: Sample Start (t0)
  // par[2]: Power
  // par[3]: Decay
  // par[4]: Pedestal
  double pedestal = par[4];
  if (x[0] < par[1])
  {
    return pedestal;
  }
  double signal = par[0] * pow((x[0] - par[1]), par[2]) * exp(-(x[0] - par[1]) * par[3]);
  return pedestal + signal;
}

double CaloWaveformFitting::SignalShape_PowerLawDoubleExp(double *x, double *par)
{
  // par[0]: Amplitude
  // par[1]: Sample Start (t0)
  // par[2]: Power
  // par[3]: Peak Time 1
  // par[4]: Pedestal
  // par[5]: Amplitude ratio
  // par[6]: Peak Time 2
  double pedestal = par[4];
  if (x[0] < par[1])
  {
    return pedestal;
  }
  double signal = par[0] * pow((x[0] - par[1]), par[2]) *
                  (((1.0 - par[5]) / pow(par[3], par[2]) * exp(par[2])) *
                       exp(-(x[0] - par[1]) * (par[2] / par[3])) +
                   (par[5] / pow(par[6], par[2]) * exp(par[2])) *
                       exp(-(x[0] - par[1]) * (par[2] / par[6])));
  return pedestal + signal;
}


double CaloWaveformFitting::SignalShape_FermiExp(double *x, double *par)
{
  // par[0]: Amplitude
  // par[1]: Midpoint (t0)
  // par[2]: Rise width (w)
  // par[3]: Decay time (tau)
  // par[4]: Pedestal

  double tt  = x[0];
  double A  = par[0];
  double t0 = par[1];
  double w  = par[2];
  double tau = par[3];
  double ped = par[4];

  if (w <= 0 || tau <= 0)
  {
    return ped;
  }

  double fermi = 1.0 / (1.0 + exp(-(tt - t0) / w));

  double expo = exp(-(tt - t0) / tau);

  if (tt < t0)
  {
    expo = 1.0;  
  }

  double signal = A * fermi * expo;

  return ped + signal;
}


std::vector<std::vector<float>> CaloWaveformFitting::calo_processing_funcfit(const std::vector<std::vector<float>> &chnlvector)
{
  std::vector<std::vector<float>> fit_values;
  int nchnls = chnlvector.size();

  for (int m = 0; m < nchnls; m++)
  {
    const std::vector<float> &v = chnlvector.at(m);
    int nsamples = v.size();

    float amp = 0;
    float time = 0;
    float ped = 0;
    float chi2 = std::numeric_limits<float>::quiet_NaN();

    // Handle zero-suppressed samples (2-sample case)
    if (nsamples == _nzerosuppresssamples)
    {
      amp = v.at(1) - v.at(0);
      time = std::numeric_limits<float>::quiet_NaN();
      ped = v.at(0);
      if (v.at(0) != 0 && v.at(1) == 0)
      {
        chi2 = 1000000;
      }
      fit_values.push_back({amp, time, ped, chi2, 0});
      continue;
    }

    // Find peak position and estimate pedestal
    float maxheight = 0;
    int maxbin = 0;
    for (int i = 0; i < nsamples; i++)
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
      pedestal = v.at(maxbin - 4);
    }
    else
    {
      pedestal = 0.5 * (v.at(nsamples - 3) + v.at(nsamples - 2));
    }

    // Software zero suppression check
    if ((_bdosoftwarezerosuppression && v.at(6) - v.at(0) < _nsoftwarezerosuppression) ||
        (_maxsoftwarezerosuppression && maxheight - pedestal < _nsoftwarezerosuppression))
    {
      amp = v.at(6) - v.at(0);
      time = std::numeric_limits<float>::quiet_NaN();
      ped = v.at(0);
      if (v.at(0) != 0 && v.at(1) == 0)
      {
        chi2 = 1000000;
      }
      fit_values.push_back({amp, time, ped, chi2, 0});
      continue;
    }

    // Create histogram for fitting
    TH1F h("h_funcfit", "", nsamples, -0.5, nsamples - 0.5);
    int ndata = 0;
    for (int i = 0; i < nsamples; ++i)
    {
      if ((v.at(i) == 16383) && _handleSaturation)
      {
        continue;
      }
      h.SetBinContent(i + 1, v.at(i));
      h.SetBinError(i + 1, 1);
      ndata++;
    }

    // If too many saturated, use all data
    if (ndata < (nsamples - 4))
    {
      ndata = nsamples;
      for (int i = 0; i < nsamples; ++i)
      {
        h.SetBinContent(i + 1, v.at(i));
        h.SetBinError(i + 1, 1);
      }
    }

    double fit_amp = 0;
    double fit_time = 0;
    double fit_ped = 0;
    double chi2val = 0;
    int npar = 0;

    if (m_funcfit_type == POWERLAWEXP)
    {
      // Create fit function with 5 parameters
      TF1 f("f_powerlaw", SignalShape_PowerLawExp, 0, nsamples, 5);
      npar = 5;

      // Set initial parameters
      double risetime = m_powerlaw_power / m_powerlaw_decay;
      double par[5];
      par[0] = maxheight - pedestal;  // Amplitude
      par[1] = maxbin - risetime;     // t0
      par[1] = std::max<double>(par[1], 0);
      par[2] = m_powerlaw_power;  // Power
      par[3] = m_powerlaw_decay;  // Decay
      par[4] = pedestal;          // Pedestal

      f.SetParameters(par);
      f.SetParLimits(0, (maxheight - pedestal) * 0.5, (maxheight - pedestal) * 10);
      f.SetParLimits(1, 0, nsamples);
      f.SetParLimits(2, 0, 10.0);
      f.SetParLimits(3, 0, 10.0);
      f.SetParLimits(4, pedestal - std::abs(maxheight - pedestal), pedestal + std::abs(maxheight - pedestal));

      // Perform fit
      h.Fit(&f, "QRN0W", "", 0, nsamples);

      // Calculate peak amplitude and time from fit parameters
      // Peak height is (p0 * Power(p2/p3, p2)) / exp(p2)
      fit_amp = (f.GetParameter(0) * pow(f.GetParameter(2) / f.GetParameter(3), f.GetParameter(2))) / exp(f.GetParameter(2));
      // Peak time is t0 + power/decay
      fit_time = f.GetParameter(1) + f.GetParameter(2) / f.GetParameter(3);
      fit_ped = f.GetParameter(4);

      // Calculate chi2
      for (int i = 0; i < nsamples; i++)
      {
        if (h.GetBinContent(i + 1) > 0)
        {
          double diff = h.GetBinContent(i + 1) - f.Eval(i);
          chi2val += diff * diff;
        }
      }
    }
    else if(m_funcfit_type == POWERLAWDOUBLEEXP) // POWERLAWDOUBLEEXP
    {
      // Create fit function with 7 parameters
      TF1 f("f_doubleexp", SignalShape_PowerLawDoubleExp, 0, nsamples, 7);
      npar = 7;

      // Set initial parameters
      double risetime = 2.0;
      double par[7];
      par[0] = (maxheight - pedestal) * 0.7;  // Amplitude
      par[1] = maxbin - risetime;             // t0
      par[1] = std::max<double>(par[1], 0);
      par[2] = m_doubleexp_power;      // Power
      par[3] = m_doubleexp_peaktime1;  // Peak Time 1
      par[4] = pedestal;               // Pedestal
      par[5] = m_doubleexp_ratio;      // Amplitude ratio
      par[6] = m_doubleexp_peaktime2;  // Peak Time 2

      f.SetParameters(par);
      f.SetParLimits(0, (maxheight - pedestal) * -1.5, (maxheight - pedestal) * 1.5);
      f.SetParLimits(1, maxbin - 3 * risetime, maxbin + risetime);
      f.SetParLimits(2, 1, 5.0);
      f.SetParLimits(3, risetime * 0.5, risetime * 4);
      f.SetParLimits(4, pedestal - std::abs(maxheight - pedestal), pedestal + std::abs(maxheight - pedestal));
      f.SetParLimits(5, 0, 1);
      f.SetParLimits(6, risetime * 0.5, risetime * 4);

      // Perform fit
      h.Fit(&f, "QRN0W", "", 0, nsamples);

      // Find peak by evaluating the function
      double peakpos1 = f.GetParameter(3);
      double peakpos2 = f.GetParameter(6);
      double max_peakpos = f.GetParameter(1) + (peakpos1 > peakpos2 ? peakpos1 : peakpos2);
      max_peakpos = std::min<double>(max_peakpos, nsamples - 1);

      fit_time = f.GetMaximumX(f.GetParameter(1), max_peakpos);
      fit_amp = f.Eval(fit_time) - f.GetParameter(4);
      fit_ped = f.GetParameter(4);

      // Calculate chi2
      for (int i = 0; i < nsamples; i++)
      {
        if (h.GetBinContent(i + 1) > 0)
        {
          double diff = h.GetBinContent(i + 1) - f.Eval(i);
          chi2val += diff * diff;
        }
      }
    }
    else if(m_funcfit_type == FERMIEXP) // POWERLAWDOUBLEEXP
    {
      TF1 f("f_fermiexp", SignalShape_FermiExp, 0, nsamples, 5);
      npar = 5;

      // Set initial parameters
      double par[5];
      par[0] = maxheight - pedestal; // Amplitude
      par[1] = maxbin ;              // t0
      par[2] = 1.0;                  // width
      par[3] = 2.0;                  // Peak Time 1
      par[4] = pedestal;             // Pedestal

      f.SetParameters(par);
      f.SetParLimits(0, maxheight-pedestal, 3*(maxheight-pedestal));
      f.SetParLimits(1, maxbin-1, maxbin);
      f.SetParLimits(2, 0.025, 2.0);
      f.SetParLimits(3, 0.5, 4.0);
      f.SetParLimits(4, pedestal-500, pedestal+500);

      f.FixParameter(2, 0.2);  

      h.Fit(&f, "QRN0W", "", 0, nsamples);

      fit_time = f.GetParameter(1);
      fit_amp = f.GetParameter(0);
      fit_ped = f.GetParameter(4);

      // Calculate chi2
      for (int i = 0; i < nsamples; i++)
      {
        if (h.GetBinContent(i + 1) > 0)
        {
          double diff = h.GetBinContent(i + 1) - f.Eval(i);
          chi2val += diff * diff;
        }
      }
    }

    int ndf = ndata - npar;
    if (ndf > 0)
    {
      chi2val /= ndf;
    }
    else
    {
      chi2val = std::numeric_limits<double>::quiet_NaN();
    }

    fit_values.push_back({static_cast<float>(fit_amp), static_cast<float>(fit_time),
                          static_cast<float>(fit_ped), static_cast<float>(chi2val), 0});
  }

  return fit_values;
}
