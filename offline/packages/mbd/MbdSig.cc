#include "MbdSig.h"

#include <phool/phool.h>

#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TMath.h>
#include <TPad.h>
#include <TSpectrum.h>
#include <TSpline.h>
#include <TTree.h>

#include <fstream>
#include <iostream>
#include <limits>

using namespace std;

MbdSig::MbdSig(const int chnum, const int nsamp)
  : _ch{chnum}
  , _nsamples{nsamp}
  , f_ampl{0}
  , f_time{0}
  , f_time_offset{4.0}
  ,  // time shift from fit
  f_integral{0.}
  ,  // time shift from fit
  hRawPulse{nullptr}
  , hSubPulse{nullptr}
  , hpulse{nullptr}
  , gRawPulse{nullptr}
  , gSubPulse{nullptr}
  , gpulse{nullptr}
  , hPed0{nullptr}
  , ped0{0}
  ,  // ped average
  ped0rms{0}
  , use_ped0{0}
  , minped0samp{-9999}
  , maxped0samp{-9999}
  , minped0x{0.}
  , maxped0x{0.}
  ,
  // time_calib{0},
  h2Template{nullptr}
  , h2Residuals{nullptr}
  ,
  // range of good amplitudes for templates
  // units are usually in ADC counts
  hAmpl{nullptr}
  , hTime{nullptr}
  , template_npointsx{0}
  , template_npointsy{0}
  , template_begintime{0}
  , template_endtime{0}
  ,
  // template_min_good_amplitude{20.},
  // template_max_good_amplitude{4080},
  // template_min_xrange{0},
  // template_max_xrange{0},
  template_fcn{nullptr}
  , _verbose{0}
{
  // cout << "In MbdSig::MbdSig(" << _ch << "," << _nsamples << ")" << endl;
}

void MbdSig::Init()
{
  TString name;

  name = "hrawpulse";
  name += _ch;
  hRawPulse = new TH1F(name, name, _nsamples, -0.5, _nsamples - 0.5);
  name = "hsubpulse";
  name += _ch;
  hSubPulse = new TH1F(name, name, _nsamples, -0.5, _nsamples - 0.5);

  // gRawPulse = new TGraphErrors(_nsamples);
  gRawPulse = new TGraphErrors();
  name = "grawpulse";
  name += _ch;
  gRawPulse->SetName(name);
  // gSubPulse = new TGraphErrors(_nsamples);
  gSubPulse = new TGraphErrors();
  name = "gsubpulse";
  name += _ch;
  gSubPulse->SetName(name);

  hpulse = hRawPulse;  // hpulse,gpulse point to raw by default
  gpulse = gRawPulse;  // we switch to sub for default if ped is applied

  //ped0stats = std::make_unique<MbdRunningStats>(100);  // use the last 100 events for running pedestal
  ped0stats = new MbdRunningStats(100);  // use the last 100 events for running pedestal
  name = "hPed0_";
  name += _ch;
  hPed0 = new TH1F(name, name, 16384, -0.5, 16383.5);
  // hPed0 = new TH1F(name,name,10000,1,0); // automatically determine the range

  SetTemplateSize(900, 1000, -10., 20.);
  // SetTemplateSize(300,300,0.,15.);

  // Set pedestal function
  ped_fcn = new TF1("ped_fcn","[0]",0,2);
  ped_fcn->SetLineColor(3);
}

void MbdSig::SetTemplateSize(const Int_t nptsx, const Int_t nptsy, const Double_t begt, const Double_t endt)
{
  template_npointsx = nptsx;
  template_npointsy = nptsy;
  template_begintime = begt;
  template_endtime = endt;

  template_y.resize(template_npointsx);
  template_yrms.resize(template_npointsx);

  Double_t xbinwid = (template_endtime - template_begintime) / (template_npointsx - 1);
  Double_t ybinwid = (1.1 + 0.1) / template_npointsy;  // yscale... should we vary this?
  if (h2Template)
  {
    delete h2Template;
  }
  if (h2Residuals)
  {
    delete h2Residuals;
  }

  TString name = "h2Template";
  name += _ch;
  h2Template = new TH2F(name, name, template_npointsx, template_begintime - xbinwid / 2., template_endtime + xbinwid / 2,
                        template_npointsy, -0.1 + ybinwid / 2.0, 1.1 + ybinwid / 2.0);

  name = "h2Residuals";
  name += _ch;
  h2Residuals = new TH2F(name, name, template_npointsx, template_begintime - xbinwid / 2., template_endtime + xbinwid / 2,
                         80, -20, 20);

  /*
  int nbins[] = { template_npointsx, nbinsy };
  Double_t lowrange[] = { template_begintime-xbinwid/2.0, -0.1+ybinwid/2.0 };
  Double_t highrange[] = { template_endtime+xbinwid/2.0, 1.1+ybinwid/2.0 };
  h2Template = new THnSparseF(name,name,2,nbins,lowrange,highrange);
  */
  // h2Template->cd( gDirectory );
}

void MbdSig::SetMinMaxFitTime(const Double_t mintime, const Double_t maxtime)
{
  fit_min_time = mintime;
  fit_max_time = maxtime;
}

MbdSig::~MbdSig()
{
  delete hRawPulse;
  delete hSubPulse;
  delete gRawPulse;
  delete gSubPulse;
  delete hPed0;
  delete ped0stats;
  // h2Template->Write();
  delete h2Template;
  delete h2Residuals;
  delete hAmpl;
  delete hTime;
  delete template_fcn;
  delete ped_fcn;
}

// This sets y, and x to sample number (starts at 0)
void MbdSig::SetY(const Float_t* y, const int invert)
{
  if (hRawPulse == nullptr)
  {
    Init();
  }

  hpulse->Reset();
  f_ampl = -9999.;
  f_time = -9999.;

  for (int isamp = 0; isamp < _nsamples; isamp++)
  {
    hRawPulse->SetBinContent(isamp + 1, y[isamp]);
    gRawPulse->SetPoint(isamp, Double_t(isamp), y[isamp]);
  }

  // Apply pedestal
  if (use_ped0 != 0 || minped0samp >= 0 || minped0x != maxped0x || ped_presamp != 0)
  {
    // cout << "sub" << endl;

    if (minped0samp >= 0)
    {
      CalcEventPed0(minped0samp, maxped0samp);
    }
    else if (minped0x != maxped0x)
    {
      CalcEventPed0(minped0x, maxped0x);
    }
    else if (ped_presamp != 0)
    {
      CalcEventPed0_PreSamp(ped_presamp, ped_presamp_nsamps);
    }

    for (int isamp = 0; isamp < _nsamples; isamp++)
    {
      hSubPulse->SetBinContent(isamp + 1, invert * (y[isamp] - ped0));
      hSubPulse->SetBinError(isamp + 1, ped0rms);
      gSubPulse->SetPoint(isamp, (Double_t) isamp, invert * (y[isamp] - ped0));
      gSubPulse->SetPointError(isamp, 0., ped0rms);
    }
  }
}

void MbdSig::SetXY(const Float_t* x, const Float_t* y, const int invert)
{
  //_verbose = 100;
  if (hRawPulse == nullptr)
  {
    Init();
  }

  hRawPulse->Reset();
  hSubPulse->Reset();
  _status = 0;

  f_ampl = -9999.;
  f_time = -9999.;

  // cout << "_nsamples " << _nsamples << endl;
  // cout << "use_ped0 " << use_ped0 << "\t" << ped0 << endl;

  for (int isamp = 0; isamp < _nsamples; isamp++)
  {
    // cout << "aaa\t" << isamp << "\t" << x[isamp] << "\t" << y[isamp] << endl;
    hRawPulse->SetBinContent(isamp + 1, y[isamp]);
    gRawPulse->SetPoint(isamp, x[isamp], y[isamp]);
    gRawPulse->SetPointError(isamp, 0, 4.0);
  }
  if ( _verbose && _ch==9 )
  {
    gRawPulse->Draw("ap");
    gRawPulse->GetHistogram()->SetTitle(gRawPulse->GetName());
    gPad->SetGridy(1);
    PadUpdate();
  }

  if (use_ped0 != 0 || minped0samp >= 0 || minped0x != maxped0x || ped_presamp != 0)
  {
    if (minped0samp >= 0)
    {
      CalcEventPed0(minped0samp, maxped0samp);
    }
    else if (minped0x != maxped0x)
    {
      CalcEventPed0(minped0x, maxped0x);
    }
    else if (ped_presamp != 0)
    {
      CalcEventPed0_PreSamp(ped_presamp, ped_presamp_nsamps);
    }

    for (int isamp = 0; isamp < _nsamples; isamp++)
    {
      // How do we handle data which is not in samples, but is in time,
      // such as DRS4 data
      if ( _verbose && isamp==(_nsamples-1) )
      {
        cout << "bbb ch " << _ch << "\t" << isamp << "\t" << x[isamp] << "\t" << invert*(y[isamp]-ped0) << endl;
      }
      hSubPulse->SetBinContent(isamp + 1, invert * (y[isamp] - ped0));
      hSubPulse->SetBinError(isamp + 1, ped0rms);
      gSubPulse->SetPoint(isamp, x[isamp], invert * (y[isamp] - ped0));
      gSubPulse->SetPointError(isamp, 0., ped0rms);
    }
    if ( _verbose && _ch==9 )
    {
      cout << "SetXY: ch " << _ch << endl;
      gSubPulse->Print("ALL");
    }
  }

  _verbose = 0;
}

Double_t MbdSig::GetSplineAmpl()
{
  if (gSubPulse == nullptr)
  {
    cout << "gsub bad " << (uint64_t) gSubPulse << endl;
    return 0.;
  }

  TSpline3 s3("s3", gSubPulse);

  // First find maximum, to rescale
  f_ampl = -999999.;
  double step_size = 0.01;
  // cout << "step size " << step_size << endl;
  for (double ix = 0; ix < _nsamples; ix += step_size)
  {
    Double_t val = s3.Eval(ix);
    if (val > f_ampl)
    {
      f_ampl = val;
    }
  }

  return f_ampl;
}

void MbdSig::FillPed0(const Int_t sampmin, const Int_t sampmax)
{
  Double_t x, y;
  for (int isamp = sampmin; isamp <= sampmax; isamp++)
  {
    gRawPulse->GetPoint(isamp, x, y);
    // gRawPulse->Print("all");
    hPed0->Fill(y);

    /*
    // chiu taken out
    ped0stats->Push( y );
    ped0 = ped0stats->Mean();
    ped0rms = ped0stats->RMS();
    */

    // cout << "ped0 " << _ch << " " << n << "\t" << ped0 << endl;
    // cout << "ped0 " << _ch << "\t" << ped0 << endl;
  }
}

void MbdSig::FillPed0(const Double_t begin, const Double_t end)
{
  Double_t x, y;
  Int_t n = gRawPulse->GetN();
  for (int isamp = 0; isamp < n; isamp++)
  {
    gRawPulse->GetPoint(isamp, x, y);
    if (x >= begin && x <= end)
    {
      hPed0->Fill(y);

      /*
         ped0stats->Push( y );
         ped0 = ped0stats->Mean();
         ped0rms = ped0stats->RMS();
         */

      // cout << "ped0 " << _ch << " " << n << "\t" << x << "\t" << y << endl;
    }

    // quit if we are past the ped region
    if (x > end)
    {
      break;
    }
  }
}

void MbdSig::SetPed0(const Double_t mean, const Double_t rms)
{
  ped0 = mean;
  ped0rms = rms;
  use_ped0 = 1;
  hpulse = hSubPulse;
  gpulse = gSubPulse;
  // if ( _ch==8 ) cout << "_ch " << _ch << " Ped = " << ped0 << endl;
}

// Get Event by Event Ped0 if requested
void MbdSig::CalcEventPed0(const Int_t minpedsamp, const Int_t maxpedsamp)
{
  // if (_ch==8) cout << "In MbdSig::CalcEventPed0(int,int)" << endl;
  hPed0->Reset();
  // ped0stats->Clear();

  Double_t x, y;
  for (int isamp = minpedsamp; isamp <= maxpedsamp; isamp++)
  {
    gRawPulse->GetPoint(isamp, x, y);

    hPed0->Fill(y);
    // ped0stats->Push( y );
    // if ( _ch==8 ) cout << "ped0stats " << isamp << "\t" << y << endl;
  }

  // use straight mean for pedestal
  // Could consider using fit to hPed0 to remove outliers
  float mean = hPed0->GetMean();
  float rms = hPed0->GetRMS();

  // SetPed0( ped0stats->Mean(), ped0stats->RMS() );

  SetPed0(mean, rms);
  // if (_ch==8) cout << "ped0stats mean, rms " << mean << "\t" << rms << endl;
}

// Get Event by Event Ped0 if requested
void MbdSig::CalcEventPed0(const Double_t minpedx, const Double_t maxpedx)
{
  hPed0->Reset();
  // ped0stats->Clear();

  Double_t x, y;
  Int_t n = gRawPulse->GetN();

  for (int isamp = 0; isamp < n; isamp++)
  {
    gRawPulse->GetPoint(isamp, x, y);

    if (x >= minpedx && x <= maxpedx)
    {
      hPed0->Fill(y);
      // ped0stats->Push( y );
    }
  }

  // use straight mean for pedestal
  // Could consider using fit to hPed0 to remove outliers
  SetPed0(hPed0->GetMean(), hPed0->GetRMS());
  /*
     Double_t mean = ped0stats->Mean();
     Double_t rms = ped0stats->RMS();
     SetPed0( mean, rms );
     */
  // cout << "ped0stats " << mean << "\t" << rms << endl;
}

// Get Event by Event Ped0, num samples before peak
// presample is number of samples before peak, nsamps is how many samples
void MbdSig::CalcEventPed0_PreSamp(const int presample, const int nsamps)
{
  //_verbose = 100;
  hPed0->Reset();
  //ped0stats->Clear();

  Double_t x, y;
  // Int_t n = gRawPulse->GetN();
  // Int_t max = gRawPulse->GetHistogram()->GetMaximumBin();
  Long64_t max = ped_presamp_maxsamp;

  // if there is no maxsamp set, use the max found in this event
  if ( ped_presamp_maxsamp == -1 )
  {
    max = TMath::LocMax(gRawPulse->GetN(), gRawPulse->GetY());
  }
  Int_t minsamp = max - presample - nsamps + 1;
  Int_t maxsamp = max - presample;
  //std::cout << "CalcEventPed0_PreSamp: ch ctr max " << _ch << "\t" << _evt_counter << "\t" << max << std::endl;

  if (minsamp < 0)
  {
    minsamp = 0;
    _status = 1;  // bad pedestal
  }
  if (maxsamp < 0)
  {
    maxsamp = 0;
    _status = 1;
  }

  if ( _verbose )
  {
    std::cout << "hPed0 presamp ch " << _ch << std::endl;
  }

  double mean = ped0stats->Mean();
  //double rms = ped0stats->RMS();
  double rms = 4.0;

  ped_fcn->SetRange(minsamp-0.1,maxsamp+0.1);
  ped_fcn->SetParameter(0,1500.);
  if ( _verbose )
  {
    gRawPulse->Fit( ped_fcn, "R" );
    gRawPulse->Draw("ap");
    ped_fcn->Draw("same");
    PadUpdate();
  }
  else
  {
    //std::cout << PHWHERE << std::endl;
    gRawPulse->Fit( ped_fcn, "RNQ" );
  }

  double chi2 = ped_fcn->GetChisquare();
  double ndf = ped_fcn->GetNDF();

  //int nstats = ped0stats->NumDataValues();

  if ( chi2/ndf < 4.0 )
  {
    mean = ped_fcn->GetParameter(0);

    for (int isamp = minsamp; isamp <= maxsamp; isamp++)
    {
      gRawPulse->GetPoint(isamp, x, y);

      // exclude outliers
      if ( fabs(y-mean) < 5.0*rms )
      {
        hPed0->Fill(y);
        ped0stats->Push( y );
      }

      if ( _verbose )
      {
        std::cout << "CalcEventPed0_PreSamp: ped0stats " << _ch << "\t" << _evt_counter << "\t" 
          << "isamp " << isamp << "\t" << x << "\t" << y << std::endl;
      }
    }
  }
  else
  {
    if ( _verbose )
    {
      if ( _evt_counter<128*100 )
      {
        std::cout << "ccc ch " << _ch << "\t" << chi2/ndf << "\t" << mean << std::endl;
      }
      std::string junk;
      cin >> junk;
    }
    // use straight mean for pedestal
    // Could consider using fit to hPed0 to remove outliers
    mean = ped0stats->Mean();
    //rms = ped0stats->RMS();
    //Double_t mean = hPed0->GetMean();
    //Double_t rms = hPed0->GetRMS();
  }

  SetPed0(mean, rms);
  if (_verbose > 0 && _evt_counter < 128*10)
  {
    std::cout << "CalcEventPed0_PreSamp: ped0stats " << _ch << "\t" << _evt_counter << "\t" << mean << "\t" << rms << std::endl;
    cout << "CalcEventPed0_PreSamp: ped0stats " << ped0stats->RMS() << std::endl;
  }
  _evt_counter++;

  _verbose = 0;
}

Double_t MbdSig::LeadingEdge(const Double_t threshold)
{
  // Find first point above threshold
  // We also make sure the next point is above threshold
  // to get rid of a high fluctuation
  int n = gSubPulse->GetN();
  Double_t* x = gSubPulse->GetX();
  Double_t* y = gSubPulse->GetY();

  int sample = -1;
  for (int isamp = 0; isamp < n; isamp++)
  {
    if (y[isamp] > threshold)
    {
      if (isamp == (n - 1) || y[isamp + 1] > threshold)
      {
        sample = isamp;
        break;
      }
    }
  }
  if (sample < 1)
  {
    return -9999.;  // no signal above threshold
  }

  // Linear Interpolation of start time
  Double_t dx = x[sample] - x[sample - 1];
  Double_t dy = y[sample] - y[sample - 1];
  Double_t dt1 = y[sample] - threshold;

  Double_t t0 = x[sample] - dt1 * (dx / dy);

  return t0;
}

Double_t MbdSig::dCFD(const Double_t fraction_threshold)
{
  // Find first point above threshold
  // We also make sure the next point is above threshold
  // to get rid of a high fluctuation
  int n = gSubPulse->GetN();
  Double_t* x = gSubPulse->GetX();
  Double_t* y = gSubPulse->GetY();

  // Get max amplitude
  Double_t ymax = TMath::MaxElement(n, y);
  if (f_ampl == -9999.)
  {
    f_ampl = ymax;
  }

  Double_t threshold = fraction_threshold * ymax;  // get fraction of amplitude
  // cout << "threshold = " << threshold << "\tymax = " << ymax <<endl;

  int sample = -1;
  for (int isamp = 0; isamp < n; isamp++)
  {
    if (y[isamp] > threshold)
    {
      if (isamp == (n - 1) || y[isamp + 1] > threshold)
      {
        sample = isamp;
        break;
      }
    }
  }
  if (sample < 1)
  {
    return -9999.;  // no signal above threshold
  }

  // Linear Interpolation of start time
  Double_t dx = x[sample] - x[sample - 1];
  Double_t dy = y[sample] - y[sample - 1];
  Double_t dt1 = y[sample] - threshold;

  Double_t t0 = x[sample] - dt1 * (dx / dy);

  return t0;
}

Double_t MbdSig::MBDTDC(const Int_t max_samp)
{
  // Get the amplitude of a fixed sample (max_samp) to get time
  // Used in MBD Time Channels
  Double_t* y = gSubPulse->GetY();

  if (y == nullptr)
  {
    std::cout << "ERROR y == 0" << std::endl;
    return NAN;
  }

  f_time = y[max_samp];

  /*
  if ( _ch==128 )
  {
    std::cout << "msig\t" << _ch << "\t" << f_time << "\t" << f_ampl << std::endl;
    gSubPulse->Print("ALL");
  }
  */

  return f_time;
}

Double_t MbdSig::Integral(const Double_t xmin, const Double_t xmax)
{
  Int_t n = gSubPulse->GetN();
  Double_t* x = gSubPulse->GetX();
  Double_t* y = gSubPulse->GetY();

  f_integral = 0.;
  for (int ix = 0; ix < n; ix++)
  {
    if (x[ix] >= xmin && x[ix] <= xmax)
    {
      // Get dx
      Double_t dx = (x[ix + 1] - x[ix - 1]) / 2.0;
      f_integral += (y[ix] * dx);
    }
  }

  return f_integral;
}

void MbdSig::LocMax(Double_t& x_at_max, Double_t& ymax, Double_t xminrange, Double_t xmaxrange)
{
  // Find index of maximum peak
  Int_t n = gSubPulse->GetN();
  Double_t* x = gSubPulse->GetX();
  Double_t* y = gSubPulse->GetY();

  // if flipped or equal, we search the whole range
  if (xmaxrange <= xminrange)
  {
    xminrange = -DBL_MAX;
    xmaxrange = DBL_MAX;
  }

  x_at_max = -DBL_MAX;
  ymax = -DBL_MAX;

  for (int i = 0; i < n; i++)
  {
    // Skip if out of range
    if (x[i] < xminrange)
    {
      continue;
    }
    if (x[i] > xmaxrange)
    {
      break;
    }

    if (y[i] > ymax)
    {
      ymax = y[i];
      x_at_max = x[i];
    }
  }
}

void MbdSig::LocMin(Double_t& x_at_max, Double_t& ymin, Double_t xminrange, Double_t xmaxrange)
{
  // Find index of minimum peak (for neg signals)
  Int_t n = gSubPulse->GetN();
  Double_t* x = gSubPulse->GetX();
  Double_t* y = gSubPulse->GetY();

  // if flipped or equal, we search the whole range
  if (xmaxrange <= xminrange)
  {
    xminrange = -DBL_MAX;
    xmaxrange = DBL_MAX;
  }

  ymin = DBL_MAX;

  for (int i = 0; i < n; i++)
  {
    // Skip if out of range
    if (x[i] < xminrange)
    {
      continue;
    }
    if (x[i] > xmaxrange)
    {
      break;
    }

    if (y[i] < ymin)
    {
      ymin = y[i];
      x_at_max = x[i];
    }
  }

  // old way of getting locmax
  // int locmax = TMath::LocMin(n,y);
}

void MbdSig::Print()
{
  Double_t x, y;
  cout << "CH " << _ch << endl;
  for (int isamp = 0; isamp < _nsamples; isamp++)
  {
    gpulse->GetPoint(isamp, x, y);
    cout << isamp << "\t" << x << "\t" << y << endl;
  }
}

void MbdSig::PadUpdate()
{
  // Make sure TCanvas is created externally!
  if ( _verbose>5 )
  {
    gPad->Modified();
    gPad->Update();
    cout << _ch << " ? ";
    if ( _verbose>10 )
    {
      std::string junk;
      cin >> junk;

      if (junk[0] == 'w' || junk[0] == 's')
      {
        TString name = "ch";
        name += _ch;
        name += ".png";
        gPad->SaveAs(name);
      }
    }
  }
}

Double_t MbdSig::TemplateFcn(const Double_t* x, const Double_t* par)
{
  // par[0] is the amplitude (relative to the spline amplitude)
  // par[1] is the start time (in sample number)
  // x[0] units are in sample number
  Double_t xx = x[0] - par[1];
  Double_t f = 0.;

  //_verbose = 100;
  if ( _verbose )
  {
    std::cout << PHWHERE << " " << _ch << " x par0 par1 " << x[0] << "\t" << par[0] << "\t" << par[1] << std::endl;
  }


  // When fit is out of limits of good part of spline, ignore fit
  if (xx < template_begintime || xx > template_endtime || isnan(xx) )
  {
    TF1::RejectPoint();

    if ( isnan(xx) )
    {
      if (_verbose > 0)
      {
	std::cout << PHWHERE << " " << _evt_counter << ", " << _ch
		  << " ERROR x par0 par1 " << x[0] << "\t" << par[0] << "\t" << par[1] << std::endl;
	if ( x[0] == 0. )
	{
	  gSubPulse->Print("ALL");
	}
      }
      return 0.;
    }

    if (xx < template_begintime)
    {
      // Double_t x0,y0;
      Double_t y0 = template_y[0];
      return par[0] * y0;
    }
    else if (xx > template_endtime)
    {
      // Double_t x0,y0;
      Double_t y0 = template_y[template_npointsx - 1];
      return par[0] * y0;
    }
  }

  // Linear Interpolation of template
  Double_t x0 = 0.;
  Double_t y0 = 0.;
  Double_t x1 = 0.;
  Double_t y1 = 0.;

  // find the index in the vector which is closest to xx
  Double_t step = (template_endtime - template_begintime) / (template_npointsx - 1);
  Double_t index = (xx - template_begintime) / step;
  //std::cout << "xxx " << index << "\t" << xx << "\t" << template_begintime << "\t" << template_endtime << "\t" << step << std::endl;

  int ilow = TMath::FloorNint(index);
  int ihigh = TMath::CeilNint(index);
  if (ilow < 0 || ihigh >= template_npointsx)
  {
    if (_verbose > 0)
    {
      cout << "ERROR, ilow ihigh " << ilow << "\t" << ihigh << endl;
      cout << " " << xx << " " << x[0] << " " << par[1] << endl;
    }

    if (ilow < 0)
    {
      ilow = 0;
    }
    else if (ihigh >= template_npointsx)
    {
      ihigh = template_npointsx - 1;
    }
  }

  if (ilow == ihigh)
  {
    f = par[0] * template_y[ilow];
  }
  else
  {
    x0 = template_begintime + ilow * step;
    y0 = template_y[ilow];
    x1 = template_begintime + ihigh * step;
    y1 = template_y[ihigh];
    f = par[0] * (y0 + ((y1 - y0) / (x1 - x0)) * (xx - x0));  // linear interpolation
  }

  // reject points with very bad rms in shape
  if (template_yrms[ilow] >= 1.0 || template_yrms[ihigh] >= 1.0)
  {
    TF1::RejectPoint();
    // return f;
  }

  // Reject points where ADC saturates
  int samp_point = static_cast<int>(x[0]);
  Double_t temp_x, temp_y;
  gRawPulse->GetPoint(samp_point, temp_x, temp_y);
  if (temp_y > 16370)
  {
    // cout << "XXXX " << _ch << "\t" << samp_point << "\t" << temp_x << "\t" << temp_y << std::endl;
    TF1::RejectPoint();
  }

  _verbose = 0;
  return f;
}

// sampmax>0 means fit to the peak near sampmax
int MbdSig::FitTemplate( const Int_t sampmax )
{
  _verbose = 0;	// uncomment to see fits
  if (_verbose > 0)
  {
    cout << "Fitting ch " << _ch << endl;
  }

  // Check if channel is empty
  if (gSubPulse->GetN() == 0)
  {
    f_ampl = 0.;
    f_time = std::numeric_limits<Float_t>::quiet_NaN();
    cout << "ERROR, gSubPulse empty" << endl;
    return 1;
  }

  // Get x and y of maximum
  Double_t x_at_max{-1.};
  Double_t ymax{0.};
  if ( sampmax>=0 )
  {
    gSubPulse->GetPoint(sampmax, x_at_max, ymax);
    x_at_max -= 2.0;
  }
  else
  {
    ymax = TMath::MaxElement( gSubPulse->GetN(), gSubPulse->GetY() );
    x_at_max = TMath::LocMax( gSubPulse->GetN(), gSubPulse->GetY() );
  }

  // Threshold cut
  if ( ymax < 10. )
  {
    f_ampl = 0.;
    f_time = std::numeric_limits<Float_t>::quiet_NaN();
    if ( _verbose>10 )
    {
      std::cout << "skipping, ymax < 10" << std::endl;
      gSubPulse->Draw("ap");
      gSubPulse->GetHistogram()->SetTitle(gSubPulse->GetName());
      gPad->SetGridy(1);
      PadUpdate();
    }
    return 1;
  }

  template_fcn->SetParameters(ymax, x_at_max);
  // template_fcn->SetParLimits(1, fit_min_time, fit_max_time);
  // template_fcn->SetParLimits(1, 3, 15);
  // template_fcn->SetRange(template_min_xrange,template_max_xrange);
  template_fcn->SetRange(0, _nsamples);

  if (_verbose == 0)
  {
    //std::cout << PHWHERE << std::endl;
    gSubPulse->Fit(template_fcn, "RNQ");
  }
  else
  {
    std::cout << "doing fit" << std::endl;
    gSubPulse->Fit(template_fcn, "R");
    gSubPulse->Draw("ap");
    gSubPulse->GetHistogram()->SetTitle(gSubPulse->GetName());
    gPad->SetGridy(1);
    PadUpdate();
    //gSubPulse->Print("ALL");
  }

  // Get fit parameters
  f_ampl = template_fcn->GetParameter(0);
  f_time = template_fcn->GetParameter(1);
  if ( f_time<0. || f_time>_nsamples )
  {
    f_time = _nsamples*0.5;  // bad fit last time
  }

  // refit with new range to exclude after-pulses
  template_fcn->SetParameters( f_ampl, f_time );
  template_fcn->SetRange( 0., f_time+4.0 );

  if (_verbose == 0)
  {
    //std::cout << PHWHERE << std::endl;
    int fit_status = gSubPulse->Fit(template_fcn, "RNQ");
    if ( fit_status<0 && _verbose )
    {
      std::cout << PHWHERE << "\t" << fit_status << std::endl;
      gSubPulse->Print("ALL");
      gSubPulse->Draw("ap");
      gSubPulse->Fit(template_fcn, "R");
      std::cout << "ampl time before refit " << f_ampl << "\t" << f_time << std::endl;
      f_ampl = template_fcn->GetParameter(0);
      f_time = template_fcn->GetParameter(1);
      std::cout << "ampl time after  refit " << f_ampl << "\t" << f_time << std::endl;
      PadUpdate();
      std::string junk;
      std::cin >> junk;
    }
  }
  else
  {
    gSubPulse->Fit(template_fcn, "R");
    //gSubPulse->Print("ALL");
    std::cout << "ampl time before refit " << f_ampl << "\t" << f_time << std::endl;
    f_ampl = template_fcn->GetParameter(0);
    f_time = template_fcn->GetParameter(1);
    std::cout << "ampl time after  refit " << f_ampl << "\t" << f_time << std::endl;
  }

  f_ampl = template_fcn->GetParameter(0);
  f_time = template_fcn->GetParameter(1);

  //std::cout << "FitTemplate " << _ch << "\t" << f_ampl << "\t" << f_time << endl;
  if (_verbose > 0 && fabs(f_ampl) > 0.)
  //if ( f_time<0 || f_time>30 )
  {
    cout << "FitTemplate " << _ch << "\t" << f_ampl << "\t" << f_time << endl;
    gSubPulse->Draw("ap");
    gSubPulse->GetHistogram()->SetTitle(gSubPulse->GetName());
    gPad->SetGridy(1);
    template_fcn->SetLineColor(4);
    template_fcn->Draw("same");
    PadUpdate();
  }

  _verbose = 0;
  return 1;
}

int MbdSig::SetTemplate(const std::vector<float>& shape, const std::vector<float>& sherr)
{
  template_y = shape;
  template_yrms = sherr;

  if (_verbose)
  {
    std::cout << "SHAPE " << _ch << "\t" << template_y.size() << std::endl;
    for (size_t i = 0; i < template_y.size(); i++)
    {
      if (i % 10 == 0)
      {
        std::cout << i << ":\t" << std::endl;
      }
      std::cout << " " << template_y[i];
    }
    std::cout << std::endl;
  }

  if (template_fcn == nullptr)
  {
    TString name = "template_fcn";
    name += _ch;
    // template_fcn = new TF1(name,this,&MbdSig::TemplateFcn,template_min_xrange,template_max_xrange,2,"MbdSig","TemplateFcn");
    // template_fcn = new TF1(name,this,&MbdSig::TemplateFcn,-10,20,2,"MbdSig","TemplateFcn");
    template_fcn = new TF1(name, this, &MbdSig::TemplateFcn, 0, _nsamples, 2, "MbdSig", "TemplateFcn");
    template_fcn->SetParameters(1, 10);
    template_fcn->SetParName(0, "ampl");
    template_fcn->SetParName(1, "time");
    SetTemplateSize(900, 1000, -10., 20.);

    if (_verbose)
    {
      std::cout << "SHAPE " << _ch << std::endl;
      template_fcn->Draw("acp");
      gPad->Modified();
      gPad->Update();
    }
  }

  return 1;
}
