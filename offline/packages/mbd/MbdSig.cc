#include "MbdSig.h"

#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH2.h>
#include <TSpline.h>
#include <TPad.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TSpectrum.h>

#include <iostream>
#include <fstream>
#include <limits>

using namespace std;


MbdSig::MbdSig(const int chnum, const int nsamp) :
  ch{chnum},
  nsamples{nsamp},
  f_ampl{0},
  f_time{0},
  f_time_offset{4.0}, // time shift from fit
  f_integral{0.},     // time shift from fit
  hRawPulse{nullptr},
  hSubPulse{nullptr},
  hpulse{nullptr},
  gRawPulse{nullptr},
  gSubPulse{nullptr},
  gpulse{nullptr},
  hPed0{nullptr},
  ped0{0},            // ped average
  ped0rms{0},
  use_ped0{0},
  minped0samp{-9999},
  maxped0samp{-9999},
  minped0x{0.},
  maxped0x{0.},
  //time_calib{0},
  h2Template{nullptr},
  h2Residuals{nullptr},
  // range of good amplitudes for templates
  // units are usually in ADC counts
  hAmpl{nullptr},
  hTime{nullptr},
  template_npointsx{0},
  template_npointsy{0},
  template_begintime{0},
  template_endtime{0},
  //template_min_good_amplitude{20.},
  //template_max_good_amplitude{4080},
  //template_min_xrange{0},
  //template_max_xrange{0},
  template_fcn{nullptr},
  verbose{0}
{
  //cout << "In MbdSig::MbdSig(" << ch << "," << nsamples << ")" << endl;

}


void MbdSig::Init()
{
  TString name;

  name = "hrawpulse"; name += ch;
  hRawPulse = new TH1F(name,name,nsamples,-0.5,nsamples-0.5);
  name = "hsubpulse"; name += ch;
  hSubPulse = new TH1F(name,name,nsamples,-0.5,nsamples-0.5);

  //gRawPulse = new TGraphErrors(nsamples);
  gRawPulse = new TGraphErrors();
  name = "grawpulse"; name += ch;
  gRawPulse->SetName(name);
  //gSubPulse = new TGraphErrors(nsamples);
  gSubPulse = new TGraphErrors();
  name = "gsubpulse"; name += ch;
  gSubPulse->SetName(name);

  hpulse = hRawPulse;   // hpulse,gpulse point to raw by default
  gpulse = gRawPulse;   // we switch to sub for default if ped is applied

  //ped0stats = new RunningStats();
  name = "hPed0_"; name += ch;
  hPed0 = new TH1F(name,name,16384,-0.5,16383.5);
  //hPed0 = new TH1F(name,name,10000,1,0); // automatically determine the range

  SetTemplateSize(900,1000,-10.,20.);
  //SetTemplateSize(300,300,0.,15.);
}

void  MbdSig::SetTemplateSize(const Int_t nptsx, const Int_t nptsy, const Double_t begt, const Double_t endt)
{
  template_npointsx = nptsx;
  template_npointsy = nptsy;
  template_begintime = begt;
  template_endtime = endt;

  template_y.resize(template_npointsx);
  template_yrms.resize(template_npointsx);

  Double_t xbinwid = (template_endtime - template_begintime)/(template_npointsx-1);
  Double_t ybinwid = (1.1+0.1)/template_npointsy;  // yscale... should we vary this?
  if ( h2Template ) delete h2Template;
  if ( h2Residuals ) delete h2Residuals;

  TString name = "h2Template"; name += ch;
  h2Template = new TH2F(name,name,template_npointsx,template_begintime-xbinwid/2.,template_endtime+xbinwid/2,
      template_npointsy,-0.1+ybinwid/2.0,1.1+ybinwid/2.0);
 
  name = "h2Residuals"; name += ch;
  h2Residuals = new TH2F(name,name,template_npointsx,template_begintime-xbinwid/2.,template_endtime+xbinwid/2,
      80,-20,20);
 
  /*
  int nbins[] = { template_npointsx, nbinsy };
  Double_t lowrange[] = { template_begintime-xbinwid/2.0, -0.1+ybinwid/2.0 };
  Double_t highrange[] = { template_endtime+xbinwid/2.0, 1.1+ybinwid/2.0 };
  h2Template = new THnSparseF(name,name,2,nbins,lowrange,highrange);
  */
  //h2Template->cd( gDirectory );

}

void  MbdSig::SetMinMaxFitTime(const Double_t mintime, const Double_t maxtime)
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
  //delete ped0stats;
  delete hPed0;
  //h2Template->Write();
  delete h2Template;
  delete h2Residuals;
  delete hAmpl;
  delete hTime;
  delete template_fcn;

}

// This sets y, and x to sample number (starts at 0)
void MbdSig::SetY(const Float_t *y, const int invert)
{
  if ( hRawPulse == nullptr )
  {
    Init();
  }

  hpulse->Reset();
  f_ampl = -9999.;
  f_time = -9999.;

  for (int isamp=0; isamp<nsamples; isamp++)
  {
    hRawPulse->SetBinContent( isamp+1, y[isamp] );
    gRawPulse->SetPoint( isamp, Double_t(isamp), y[isamp] );
  }

  // Apply pedestal
  if ( use_ped0 != 0 || minped0samp >= 0 || minped0x != maxped0x || ped_presamp != 0 )
  {
    //cout << "sub" << endl;

    if ( minped0samp >= 0 )
    {
      CalcEventPed0(minped0samp,maxped0samp);
    }
    else if ( minped0x != maxped0x )
    {
      CalcEventPed0(minped0x,maxped0x);
    }
    else if ( ped_presamp != 0 )
    {
      CalcEventPed0_PreSamp(ped_presamp,ped_presamp_nsamps);
    }

    for (int isamp=0; isamp<nsamples; isamp++)
    {
      hSubPulse->SetBinContent( isamp+1, invert*(y[isamp]-ped0) );
      hSubPulse->SetBinError( isamp+1, ped0rms );
      gSubPulse->SetPoint( isamp, (Double_t)isamp, invert*(y[isamp]-ped0) );
      gSubPulse->SetPointError( isamp, 0., ped0rms );
    }
  }
}

void MbdSig::SetXY(const Float_t *x, const Float_t *y, const int invert)
{
  if ( hRawPulse == nullptr )
  {
    Init();
  }

  hRawPulse->Reset();
  hSubPulse->Reset();
  _status = 0;

  f_ampl = -9999.;
  f_time = -9999.;

  //cout << "nsamples " << nsamples << endl;
  //cout << "use_ped0 " << use_ped0 << "\t" << ped0 << endl;

  for( int isamp=0; isamp<nsamples; isamp++ )
  {
    //cout << "aaa\t" << isamp << "\t" << x[isamp] << "\t" << y[isamp] << endl;
    hRawPulse->SetBinContent( isamp+1, y[isamp] );
    gRawPulse->SetPoint( isamp, x[isamp], y[isamp] );
  }

  if ( use_ped0 != 0 || minped0samp >= 0 || minped0x != maxped0x || ped_presamp!=0 )
  {
    if ( minped0samp >= 0 )
    {
      CalcEventPed0(minped0samp,maxped0samp);
    }
    else if ( minped0x != maxped0x )
    {
      CalcEventPed0(minped0x,maxped0x);
    }
    else if ( ped_presamp != 0 )
    {
      CalcEventPed0_PreSamp(ped_presamp,ped_presamp_nsamps);
    }

    for (int isamp=0; isamp<nsamples; isamp++)
    {
      // How do we handle data which is not in samples, but is in time,
      // such as DRS4 data
      //if ( isamp==(nsamples-1) ) cout << "bbb\t" << isamp << "\t" << x[isamp] << "\t" << invert*(y[isamp]-ped0) << endl;
      hSubPulse->SetBinContent( isamp+1, invert*(y[isamp]-ped0) );
      hSubPulse->SetBinError( isamp+1, ped0rms );
      gSubPulse->SetPoint( isamp, x[isamp], invert*(y[isamp]-ped0) );
      gSubPulse->SetPointError( isamp, 0., ped0rms );
    }
  }
}

Double_t MbdSig::GetSplineAmpl()
{
  if ( gSubPulse==nullptr)
  {
    cout << "gsub bad " << (uint64_t)gSubPulse << endl;
    return 0.;
  }

  TSpline3 s3("s3",gSubPulse);

  // First find maximum, to rescale
  f_ampl = -999999.;
  double step_size = 0.01;
  //cout << "step size " << step_size << endl;
  for (double ix=0; ix<nsamples; ix += step_size)
  {
    Double_t val = s3.Eval(ix);
    if ( val > f_ampl )
    {
      f_ampl = val;
    }
  }

  return f_ampl;
}

void MbdSig::FillPed0(const Int_t sampmin, const Int_t sampmax)
{
  Double_t x, y;
  for (int isamp=sampmin; isamp<=sampmax; isamp++)
  {
    gRawPulse->GetPoint(isamp,x,y);
    //gRawPulse->Print("all");
    hPed0->Fill( y );

    /*
    // chiu taken out
    ped0stats->Push( y );
    ped0 = ped0stats->Mean();
    ped0rms = ped0stats->RMS();
    */

    //cout << "ped0 " << ch << " " << n << "\t" << ped0 << endl;
    //cout << "ped0 " << ch << "\t" << ped0 << endl;
  }

}


void MbdSig::FillPed0(const Double_t begin, const Double_t end)
{
  Double_t x, y;
  Int_t n = gRawPulse->GetN();
  for (int isamp=0; isamp<n; isamp++)
  {
    gRawPulse->GetPoint(isamp,x,y);
    if ( x>=begin && x<=end )
    {
      hPed0->Fill( y );

      /*
         ped0stats->Push( y );
         ped0 = ped0stats->Mean();
         ped0rms = ped0stats->RMS();
         */

      //cout << "ped0 " << ch << " " << n << "\t" << x << "\t" << y << endl;
    }

    // quit if we are past the ped region
    if ( x>end ) break;
  }

}


void MbdSig::SetPed0(const Double_t mean, const Double_t rms)
{
  ped0 = mean;
  ped0rms = rms;
  use_ped0 = 1;
  hpulse = hSubPulse;
  gpulse = gSubPulse;
  //if ( ch==8 ) cout << "ch " << ch << " Ped = " << ped0 << endl;
}

// Get Event by Event Ped0 if requested
void MbdSig::CalcEventPed0(const Int_t minpedsamp, const Int_t maxpedsamp)
{
  //if (ch==8) cout << "In MbdSig::CalcEventPed0(int,int)" << endl;
  hPed0->Reset();
  //ped0stats->Clear();

  Double_t x, y;
  for (int isamp=minpedsamp; isamp<=maxpedsamp; isamp++)
  {
    gRawPulse->GetPoint(isamp,x,y);

    hPed0->Fill(y);
    //ped0stats->Push( y );
    //if ( ch==8 ) cout << "ped0stats " << isamp << "\t" << y << endl;
  }


  // use straight mean for pedestal
  // Could consider using fit to hPed0 to remove outliers
  float mean = hPed0->GetMean();
  float rms = hPed0->GetRMS();

  //SetPed0( ped0stats->Mean(), ped0stats->RMS() );

  SetPed0( mean, rms );
  //if (ch==8) cout << "ped0stats mean, rms " << mean << "\t" << rms << endl;
}

// Get Event by Event Ped0 if requested
void MbdSig::CalcEventPed0(const Double_t minpedx, const Double_t maxpedx)
{
  hPed0->Reset();
  //ped0stats->Clear();

  Double_t x, y;
  Int_t n = gRawPulse->GetN();

  for (int isamp=0; isamp<n; isamp++)
  {
    gRawPulse->GetPoint(isamp,x,y);

    if ( x>= minpedx && x<= maxpedx)
    {
      hPed0->Fill(y);
      //ped0stats->Push( y );
    }
  }

  // use straight mean for pedestal
  // Could consider using fit to hPed0 to remove outliers
  SetPed0( hPed0->GetMean(), hPed0->GetRMS() );
  /*
     Double_t mean = ped0stats->Mean();
     Double_t rms = ped0stats->RMS();
     SetPed0( mean, rms );
     */
  //cout << "ped0stats " << mean << "\t" << rms << endl;
}

// Get Event by Event Ped0, num samples before peak
// presample is number of samples before peak, nsamps is how many samples
void MbdSig::CalcEventPed0_PreSamp(const int presample, const int nsamps)
{
  hPed0->Reset();
  //ped0stats->Clear();

  Double_t x, y;
  //Int_t n = gRawPulse->GetN();
  //Int_t max = gRawPulse->GetHistogram()->GetMaximumBin();
  Long64_t max = TMath::LocMax(gRawPulse->GetN(),gRawPulse->GetY());
  Int_t minsamp = max - presample - nsamps + 1;
  Int_t maxsamp = max - presample;
  //cout << "CalcEventPed0_PreSamp: " << max << endl;

  if ( minsamp<0 )
  {
    minsamp = 0;
    _status = 1;  // bad pedestal
  }
  if ( maxsamp<0 )
  {
    maxsamp = 0;
    _status = 1;
  }

  for (int isamp=minsamp; isamp<=maxsamp; isamp++)
  {
    gRawPulse->GetPoint(isamp,x,y);

    hPed0->Fill(y);
    //ped0stats->Push( y );
  }

  // use straight mean for pedestal
  // Could consider using fit to hPed0 to remove outliers
  //Double_t mean = ped0stats->Mean();
  //Double_t rms = ped0stats->RMS();
  Double_t mean = hPed0->GetMean();
  Double_t rms = hPed0->GetRMS();
  SetPed0( mean, rms );
  static int counter = 0;
  if (verbose>0 && counter<10)
  {
    cout << "CalcEventPed0_PreSamp: ped0stats " << mean << "\t" << rms << endl;
    counter++;
  }
}

Double_t MbdSig::LeadingEdge(const Double_t threshold)
{
  // Find first point above threshold
  // We also make sure the next point is above threshold
  // to get rid of a high fluctuation
  int n = gSubPulse->GetN();
  Double_t *x = gSubPulse->GetX();
  Double_t *y = gSubPulse->GetY();

  int sample = -1;
  for (int isamp=0; isamp<n; isamp++)
  {
    if ( y[isamp] > threshold )
    {
      if ( isamp==(n-1) || y[isamp+1] > threshold )
      {
        sample = isamp;
        break;
      }
    }
  }
  if ( sample < 1 ) return -9999.;  // no signal above threshold

  // Linear Interpolation of start time
  Double_t dx = x[sample] - x[sample-1];
  Double_t dy = y[sample] - y[sample-1];
  Double_t dt1 = y[sample] - threshold;

  Double_t t0 = x[sample] - dt1*(dx/dy);

  return t0;
}

Double_t MbdSig::dCFD(const Double_t fraction_threshold)
{
  // Find first point above threshold
  // We also make sure the next point is above threshold
  // to get rid of a high fluctuation
  int n = gSubPulse->GetN();
  Double_t *x = gSubPulse->GetX();
  Double_t *y = gSubPulse->GetY();

  // Get max amplitude
  Double_t ymax = TMath::MaxElement(n,y);
  if ( f_ampl == -9999. ) f_ampl = ymax;

  Double_t threshold = fraction_threshold * ymax; // get fraction of amplitude
  //cout << "threshold = " << threshold << "\tymax = " << ymax <<endl;

  int sample = -1;
  for (int isamp=0; isamp<n; isamp++)
  {
    if ( y[isamp] > threshold )
    {
      if ( isamp==(n-1) || y[isamp+1] > threshold )
      {
        sample = isamp;
        break;
      }
    }
  }
  if ( sample < 1 ) return -9999.;  // no signal above threshold

  // Linear Interpolation of start time
  Double_t dx = x[sample] - x[sample-1];
  Double_t dy = y[sample] - y[sample-1];
  Double_t dt1 = y[sample] - threshold;

  Double_t t0 = x[sample] - dt1*(dx/dy);

  return t0;
}

Double_t MbdSig::MBD(const Int_t max_samp)
{
  // Get the amplitude of the sample number to get time
  Double_t *y = gSubPulse->GetY();

  if ( y==0 ) { 
    std::cout << "ERROR y == 0" << std::endl; 
    return NAN;
  }

  // SHOULD INCLUDE TIME CALIBRATION HERE
  f_time = y[max_samp];

  // Get max amplitude, and set it if it hasn't already been set
  int n = gSubPulse->GetN();
  Double_t ymax = TMath::MaxElement(n,y);
  if ( f_ampl == -9999. ) f_ampl = ymax;

  return f_time;
}

Double_t MbdSig::Integral(const Double_t xmin, const Double_t xmax)
{
  Int_t n = gSubPulse->GetN();
  Double_t* x = gSubPulse->GetX();
  Double_t* y = gSubPulse->GetY();

  f_integral = 0.;
  for (int ix=0; ix<n; ix++)
  {
    if (x[ix]>=xmin && x[ix]<=xmax)
    {
      // Get dx
      Double_t dx = (x[ix+1]-x[ix-1])/2.0;
      f_integral += (y[ix]*dx);
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
  if ( xmaxrange <= xminrange )
  {
    xminrange = -DBL_MAX;
    xmaxrange = DBL_MAX;
  }

  x_at_max = -DBL_MAX;
  ymax = -DBL_MAX;

  for (int i=0; i<n; i++)
  {
    // Skip if out of range
    if ( x[i] < xminrange ) continue;
    if ( x[i] > xmaxrange ) break;

    if ( y[i] > ymax )
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
  if ( xmaxrange <= xminrange )
  {
    xminrange = -DBL_MAX;
    xmaxrange = DBL_MAX;
  }

  ymin = DBL_MAX;

  for (int i=0; i<n; i++)
  {
    // Skip if out of range
    if ( x[i] < xminrange ) continue;
    if ( x[i] > xmaxrange ) break;

    if ( y[i] < ymin )
    {
      ymin = y[i];
      x_at_max = x[i];
    }
  }

  // old way of getting locmax
  //int locmax = TMath::LocMin(n,y);
}

void MbdSig::Print()
{
  Double_t x, y;
  cout << "CH " << ch << endl;
  for (int isamp=0; isamp<nsamples; isamp++)
  {
    gpulse->GetPoint(isamp,x,y);
    cout << isamp << "\t" << x << "\t" << y << endl;
  }
}

void MbdSig::PadUpdate()
{
  // Make sure TCanvas is created externally!
  gPad->Modified();
  gPad->Update();
  cout << ch << " ? ";
  TString junk;
  cin >> junk;

  if (junk[0] == 'w' || junk[0] == 's')
  {
    TString name = "ch"; name += ch; name += ".png";
    gPad->SaveAs( name );
  }
}

Double_t MbdSig::TemplateFcn(const Double_t *x, const Double_t *par)
{
  // par[0] is the amplitude (relative to the spline amplitude)
  // par[1] is the start time (in sample number)
  // x[0] units are in sample number
  Double_t xx = x[0]-par[1];
  Double_t f = 0.;

  //verbose = 100;

  // When fit is out of limits of good part of spline, ignore fit
  if ( xx<template_begintime || xx>template_endtime )
  {
    TF1::RejectPoint();
    if ( xx < template_begintime )
    {
      //Double_t x0,y0;
      Double_t y0 = template_y[0];
      return par[0]*y0;
    }
    else if ( xx > template_endtime )
    {
      //Double_t x0,y0;
      Double_t y0 = template_y[template_npointsx-1];
      return par[0]*y0;
    }
  }

  // Linear Interpolation of template
  Double_t x0 = 0.;
  Double_t y0 = 0.;
  Double_t x1 = 0.;
  Double_t y1 = 0.;

  // find the index in the vector which is closest to xx
  Double_t step = (template_endtime - template_begintime) / (template_npointsx-1);
  Double_t index = (xx - template_begintime)/step;

  int ilow = TMath::FloorNint( index );
  int ihigh = TMath::CeilNint( index );
  if ( ilow < 0 || ihigh >= template_npointsx )
  {
    if ( verbose>0 )
    {
      cout << "ERROR, ilow ihigh " << ilow << "\t" << ihigh << endl;
      cout << " " << xx << " " << x[0] << " " << par[1] << endl;
    }

    if ( ilow<0 )
    {
      ilow = 0;
    }
    else if ( ihigh >= template_npointsx )
    {
      ihigh = template_npointsx - 1;
    }
  }

  if ( ilow==ihigh )
  {
    f = par[0]*template_y[ilow];
  }
  else
  {
    x0 = template_begintime + ilow*step;
    y0 = template_y[ilow];
    x1 = template_begintime + ihigh*step;
    y1 = template_y[ihigh];
    f = par[0]*(y0+((y1-y0)/(x1-x0))*(xx-x0));  // linear interpolation
  }

  // reject points with very bad rms in shape
  if ( template_yrms[ilow]>=1.0 || template_yrms[ihigh]>=1.0 )
  {
    TF1::RejectPoint();
    //return f;
  }

  // Reject points where ADC saturates
  int samp_point = static_cast<int>( x[0] );
  Double_t temp_x, temp_y;
  gRawPulse->GetPoint( samp_point, temp_x, temp_y );
  if ( temp_y > 16370 )
  {
    //cout << "XXXX " << ch << "\t" << samp_point << "\t" << temp_x << "\t" << temp_y << std::endl;
    TF1::RejectPoint();
  }

  //verbose = 0;
  return f;
}

int MbdSig::FitTemplate()
{
  //verbose = 100;	// uncomment to see fits
  if ( verbose>0 ) cout << "Fitting ch " << ch << endl;

  // Check if channel is empty
  if ( gSubPulse->GetN() == 0 )
  {
    f_ampl = -9999.;
    f_time = -9999.;
    //cout << "gSubPulse empty" << endl;
    return 1;
  }

  // Get x-position of maximum
  Double_t x_at_max, ymax;
  //LocMax(x_at_max, ymax);

  Int_t xsamp = (fit_min_time + fit_max_time)/2 + 2;  // use max samp
  gSubPulse->GetPoint( xsamp, x_at_max, ymax );

  template_fcn->SetParameters(ymax, x_at_max);
  //template_fcn->SetParLimits(1, fit_min_time, fit_max_time);
  //template_fcn->SetParLimits(1, 3, 15);
  //template_fcn->SetRange(template_min_xrange,template_max_xrange);
  template_fcn->SetRange(0,nsamples);

  if ( verbose==0 ) gSubPulse->Fit(template_fcn,"RNQ");
  else              gSubPulse->Fit(template_fcn,"R");

  // Get fit parameters
  f_ampl = template_fcn->GetParameter(0);
  f_time = template_fcn->GetParameter(1);

  if ( verbose>0 && fabs(f_ampl) > 0. )
  {
    cout << "FitTemplate " << ch << "\t" << f_ampl << "\t" << f_time << endl;
    gSubPulse->Draw("ap");
    template_fcn->SetLineColor(4);
    template_fcn->Draw("same");
    PadUpdate();
  }

  //verbose = 0;
  return 1;
}

int MbdSig::SetTemplate(const std::vector<float>& shape, const std::vector<float>& sherr)
{
  template_y = shape;
  template_yrms = sherr;

  if ( verbose )
  {
    std::cout << "SHAPE " << ch << "\t" << template_y.size() << std::endl;
    for ( size_t i=0; i<template_y.size(); i++)
    {
      if ( i%10 == 0 ) std::cout << i << ":\t" << std::endl;
      std::cout << " " << template_y[i];
    }
    std::cout << std::endl;
  }

  if ( template_fcn == nullptr )
  {
    TString name = "template_fcn"; name += ch;
    //template_fcn = new TF1(name,this,&MbdSig::TemplateFcn,template_min_xrange,template_max_xrange,2,"MbdSig","TemplateFcn");
    //template_fcn = new TF1(name,this,&MbdSig::TemplateFcn,-10,20,2,"MbdSig","TemplateFcn");
    template_fcn = new TF1(name,this,&MbdSig::TemplateFcn,0,nsamples,2,"MbdSig","TemplateFcn");
    template_fcn->SetParameters(1,10);
    SetTemplateSize(900,1000,-10.,20.);

    if ( verbose )
    {
      std::cout << "SHAPE " << ch << std::endl;
      template_fcn->Draw("acp");
      gPad->Modified();
      gPad->Update();
    }
  }

  return 1;
}

