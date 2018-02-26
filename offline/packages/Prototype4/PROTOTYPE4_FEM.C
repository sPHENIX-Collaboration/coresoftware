// $Id: $

/*!
 * \file PROTOTYPE4_FEM.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PROTOTYPE4_FEM.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TVirtualFitter.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

int PROTOTYPE4_FEM::GetChannelNumber(std::string caloname, int i_column, int i_row)
{
  if (caloname == "HCALIN")
  {
    //    > Anthony Hodges
    //    > PhD. Student, Georgia State University
    //    > Nuclear and High Energy Physics
    //    > ahodges21@student.gsu.edu

    static const int hbdchanIHC[4][4] = {{4, 8, 12, 16},
                                         {3, 7, 11, 15},
                                         {2, 6, 10, 14},
                                         {1, 5, 9, 13}};

    assert(i_row >= 0);
    assert(i_row < NCH_IHCAL_ROWS);
    assert(i_column >= 0);
    assert(i_column < NCH_IHCAL_COLUMNS);

    return hbdchanIHC[i_row][i_column] + 64;
  }
  else if (caloname == "HCALOUT")
  {
    // Place holder

    return 112 + 8 * i_column + 2 * i_row;
  }
  else if (caloname == "EMCAL")
  {
    assert(i_row >= 0);
    assert(i_row < NCH_EMCAL_ROWS);
    assert(i_column >= 0);
    assert(i_column < NCH_EMCAL_COLUMNS);

    //    > Anthony Hodges
    //    > PhD. Student, Georgia State University
    //    > Nuclear and High Energy Physics
    //    > ahodges21@student.gsu.edu

    // mapping taken from John Haggerty's emcalall.C found here:
    // /gpfs/mnt/gpfs02/sphenix/data/data01/caladc/wd/wd409/macros
    // This map here takes in a channel number and gives you the corresponding canvas position
    // Presumably we want the opposite, put in canvas position, output channel number
    // So we'll work backwards
    //    Float_t canmap[64] = {6, 7, 14, 15, 4, 5, 12, 13, 2, 3, 10, 11, 0, 1, 8, 9, 22, 23, 30, 31, 20, 21, 28, 29,
    //                          18, 19, 26, 27, 16, 17, 24, 25, 38, 39, 46, 47, 36, 37, 44, 45, 34, 35, 42, 43,
    //                          32, 33, 40, 41, 54, 55, 62, 63, 52, 53, 60, 61, 50, 51, 58, 59, 48, 49, 56, 57};
    //    static int hbdchanEMC[8][8];  //I'm gonna fill this boy with the above channels
    //    for (int chan = 0; chan < 64; chan++)
    //    {
    //      hbdchanEMC[(int) floor(canmap[chan] / 8)][((int) canmap[chan]) % 8] = chan;
    //    }

    // Revision from Martin with static reverse:
    //    This is now lining up towers from 0....63, and tells you
    //
    //    for tower i, what is the actual ADC index I have to go to?
    //
    //    tower[0] -> chvector[0] = 12  is then the adc channel nr.
    //
    //
    //
    //
    //
    //
    //     static const int  chvector[]=  {12 , 13 , 8 , 9 , 4 ,5 ,0 ,1 ,14 ,15
    //    ,10 ,11 ,6 ,7 ,2 ,3 ,28 ,29 ,24 ,25 ,20 ,21 ,
    //                                      16 ,17 ,30 ,31 ,26 ,27 ,22 ,23 ,18 ,19
    //    ,44 ,45 ,40 ,41 ,36 ,37 ,32 ,33 ,46 ,47 ,42 ,43 ,38 ,
    //                                      39 ,34 ,35 ,60 ,61 ,56 ,57 ,52 ,53 ,48
    //    ,49 ,62 ,63 ,58 ,59 ,54 ,55 ,50 ,51};
    const static int canmap[64] = {
        12, 13, 8, 9, 4, 5, 0, 1, 14, 15, 10, 11, 6, 7, 2, 3, 28, 29, 24, 25, 20, 21,
        16, 17, 30, 31, 26, 27, 22, 23, 18, 19, 44, 45, 40, 41, 36, 37, 32, 33, 46, 47, 42, 43, 38,
        39, 34, 35, 60, 61, 56, 57, 52, 53, 48, 49, 62, 63, 58, 59, 54, 55, 50, 51
    };

    const int linear_channel_ID = (NCH_EMCAL_ROWS - i_row - 1) * NCH_EMCAL_COLUMNS + i_column;

    assert(linear_channel_ID >= 0);
    assert(linear_channel_ID < 64);

    return canmap[linear_channel_ID];
  }

  std::cout << "PROTOTYPE4_FEM::GetHBDCh - invalid input caloname " << caloname
            << " i_column " << i_column << " i_row " << i_row << std::endl;
  exit(1);
  return -9999;
}

bool PROTOTYPE4_FEM::SampleFit_PowerLawExp(  //
    const std::vector<double> &samples,      //
    double &peak,                            //
    double &peak_sample,                     //
    double &pedstal,                         //
    const int verbosity)
{
  int peakPos = 0.;

  assert(samples.size() == NSAMPLES);

  TGraph gpulse(NSAMPLES);
  for (int i = 0; i < NSAMPLES; i++)
  {
    (gpulse.GetX())[i] = i;

    (gpulse.GetY())[i] = samples[i];
  }

  double pedestal = gpulse.GetY()[0];  //(double) PEDESTAL;
  double peakval = pedestal;
  const double risetime = 4;

  for (int iSample = 0; iSample < NSAMPLES; iSample++)
  {
    if (abs(gpulse.GetY()[iSample] - pedestal) > abs(peakval - pedestal))
    {
      peakval = gpulse.GetY()[iSample];
      peakPos = iSample;
    }
  }
  peakval -= pedestal;

  // fit function
  TF1 fits("f_SignalShape_PowerLawExp", SignalShape_PowerLawExp, 0., NSAMPLES, 5);

  double par[6] =
      {0};
  par[0] = peakval;  // /3.;
  par[1] = peakPos - risetime;
  if (par[1] < 0.)
    par[1] = 0.;
  par[2] = 4.;
  par[3] = 1.5;
  par[4] = pedestal;
  par[5] = 0;
  fits.SetParameters(par);
  fits.SetParNames("Amplitude", "Peak Sample", "Power", "Decay", "Pedstal", "Baseline shift");
  fits.SetParLimits(0, peakval * 0.5, peakval * 10);
  fits.SetParLimits(1, 0, NSAMPLES);
  fits.SetParLimits(2, 0, 10.);
  fits.SetParLimits(3, 0, 10);
  fits.SetParLimits(4, pedestal - abs(peakval), pedestal + abs(peakval));
  //  fits.SetParLimits(5, - abs(peakval),  + abs(peakval));
  //  fits.FixParameter(5, 0);

  //Saturation correction - Abhisek
  for (int ipoint = 0; ipoint < gpulse.GetN(); ipoint++)
    if ((gpulse.GetY())[ipoint] <= 10 or (gpulse.GetY())[ipoint] >= ((1 << 14) - 10))  // drop point if touching max or low limit on ADCs
    {
      gpulse.RemovePoint(ipoint);
      ipoint--;
    }

  if (verbosity <= 1)
    gpulse.Fit(&fits, "MQRN0W", "goff", 0., (double) NSAMPLES);
  else
    gpulse.Fit(&fits, "MRN0VW", "goff", 0., (double) NSAMPLES);

  if (verbosity)
  {
    static int id = 0;
    ++id;

    string c_name(string("PROTOTYPE4_FEM_SampleFit_PowerLawExp_") + to_string(id));

    TCanvas *canvas = new TCanvas(
        c_name.c_str(), c_name.c_str());
    canvas->Update();

    TGraph *g_plot = static_cast<TGraph *>(gpulse.DrawClone("ap*l"));
    g_plot->SetTitle("ADC data and fit;Sample number;ADC value");
    fits.DrawClone("same");
    fits.Print();
    canvas->Update();
    //    sleep(1);
  }

  //  peak = fits.GetParameter(0); // not exactly peak height
  peak = (fits.GetParameter(0) * pow(fits.GetParameter(2) / fits.GetParameter(3), fits.GetParameter(2))) / exp(fits.GetParameter(2));  // exact peak height is (p0*Power(p2/p3,p2))/Power(E,p2)

  //  peak_sample = fits.GetParameter(1); // signal start time
  peak_sample = fits.GetParameter(1) + fits.GetParameter(2) / fits.GetParameter(3);  // signal peak time

  // peak integral = p0*Power(p3,-1 - p2)*Gamma(1 + p2). Note yet used in output

  pedstal = fits.GetParameter(4);

  return true;
}

double
PROTOTYPE4_FEM::SignalShape_PowerLawExp(double *x, double *par)
{
  double pedestal = par[4];
  //                        + ((x[0] - 1.5 * par[1]) > 0) * par[5];  // quick fix on exting tails on the signal function
  if (x[0] < par[1])
    return pedestal;
  //double  signal = (-1)*par[0]*pow((x[0]-par[1]),par[2])*exp(-(x[0]-par[1])*par[3]);
  double signal = par[0] * pow((x[0] - par[1]), par[2]) * exp(-(x[0] - par[1]) * par[3]);
  return pedestal + signal;
}

bool PROTOTYPE4_FEM::SampleFit_PowerLawDoubleExp(  //
    const std::vector<double> &samples,            //
    double &peak,                                  //
    double &peak_sample,                           //
    double &pedstal,                               //
    const int verbosity)
{
  int peakPos = 0.;

  assert(samples.size() == NSAMPLES);

  TGraph gpulse(NSAMPLES);
  for (int i = 0; i < NSAMPLES; i++)
  {
    (gpulse.GetX())[i] = i;

    (gpulse.GetY())[i] = samples[i];
  }

  double pedestal = gpulse.GetY()[0];  //(double) PEDESTAL;
  double peakval = pedestal;
  const double risetime = 4;

  for (int iSample = 0; iSample < NSAMPLES; iSample++)
  {
    if (abs(gpulse.GetY()[iSample] - pedestal) > abs(peakval - pedestal))
    {
      peakval = gpulse.GetY()[iSample];
      peakPos = iSample;
    }
  }
  peakval -= pedestal;

  // fit function
  TF1 fits("f_SignalShape_PowerLawDoubleExp", SignalShape_PowerLawDoubleExp, 0., NSAMPLES, 7);

  double par[10] =
      {0};
  par[0] = peakval;  // /3.;
  par[1] = peakPos - risetime;
  if (par[1] < 0.)
    par[1] = 0.;
  par[2] = 2.;
  par[3] = 1.5;
  par[4] = pedestal;
  par[5] = peakval / 20;
  par[6] = 1;
  fits.SetParameters(par);
  fits.SetParNames("Amplitude1", "Peak Sample", "Power", "Decay1", "Pedestal", "Amplitude2", "Decay2");
  fits.SetParLimits(0, peakval * 0.5, peakval * 5);
  fits.SetParLimits(1, 0, NSAMPLES);
  fits.SetParLimits(2, 1, 10.);
  fits.SetParLimits(3, 1., 2.5);
  fits.SetParLimits(4, pedestal - abs(peakval), pedestal + abs(peakval));
  fits.SetParLimits(5, 0, peakval * 1.5);
  fits.SetParLimits(6, 0, 1);

  //Saturation correction - Abhisek
  for (int ipoint = 0; ipoint < gpulse.GetN(); ipoint++)
    if ((gpulse.GetY())[ipoint] <= 10 or (gpulse.GetY())[ipoint] >= ((1 << 14) - 10))  // drop point if touching max or low limit on ADCs
    {
      gpulse.RemovePoint(ipoint);
      ipoint--;
    }

  if (verbosity <= 1)
    gpulse.Fit(&fits, "MQRN0W", "goff", 0., (double) NSAMPLES);
  else
    gpulse.Fit(&fits, "MRN0VW", "goff", 0., (double) NSAMPLES);

  if (verbosity)
  {
    static int id = 0;
    ++id;

    string c_name(string("PROTOTYPE4_FEM_SampleFit_PowerLawDoubleExp_") + to_string(id));

    TCanvas *canvas = new TCanvas(
        c_name.c_str(), c_name.c_str());
    canvas->Update();

    TGraph *g_plot = static_cast<TGraph *>(gpulse.DrawClone("ap*l"));
    g_plot->SetTitle("ADC data and fit;Sample number;ADC value");

    fits.SetLineColor(kMagenta);
    fits.DrawClone("same");
    fits.Print();

    TF1 f1("f_SignalShape_PowerLawExp1", SignalShape_PowerLawExp, 0., NSAMPLES, 5);
    f1.SetParameters(
        fits.GetParameter(0),
        fits.GetParameter(1),
        fits.GetParameter(2),
        fits.GetParameter(3),
        fits.GetParameter(4));
    f1.SetLineColor(kBlue);
    f1.DrawClone("same");

    TF1 f2("f_SignalShape_PowerLawExp2", SignalShape_PowerLawExp, 0., NSAMPLES, 5);
    f2.SetParameters(
        fits.GetParameter(5),
        fits.GetParameter(1),
        fits.GetParameter(2),
        fits.GetParameter(6),
        fits.GetParameter(4));
    f2.SetLineColor(kRed);
    f2.DrawClone("same");

    canvas->Update();

    //    sleep(1);
  }

  //  peak = fits.GetParameter(0); // not exactly peak height
  //  peak = (fits.GetParameter(0) * pow(fits.GetParameter(2) / fits.GetParameter(3), fits.GetParameter(2)))    //
  //             / exp(fits.GetParameter(2))                                                                    //
  //         + (fits.GetParameter(5) * pow(fits.GetParameter(2) / fits.GetParameter(6), fits.GetParameter(2)))  //
  //               / exp(fits.GetParameter(2));
  // exact peak height is (p0*Power(p2/p3,p2))/Power(E,p2)

  //  peak_sample = fits.GetParameter(1); // signal start time
  //  peak_sample = fits.GetParameter(1) + fits.GetParameter(2) / fits.GetParameter(3);  // signal peak time

  // peak integral = p0*Power(p3,-1 - p2)*Gamma(1 + p2). Note yet used in output

  pedstal = fits.GetParameter(4);
  peak_sample = fits.GetMaximumX(fits.GetParameter(1), fits.GetParameter(1) + fits.GetParameter(2) / fits.GetParameter(6));
  peak = fits.Eval(peak_sample) - pedstal;

  if (verbosity)
  {
    cout << "PROTOTYPE4_FEM::SampleFit_PowerLawDoubleExp - "
         << "peak_sample = " << peak_sample << ", "
         << "peak = " << peak << ", "
         << "pedstal = " << pedstal << endl;
  }

  return true;
}

double
PROTOTYPE4_FEM::SignalShape_PowerLawDoubleExp(double *x, double *par)
{
  double pedestal = par[4];
  //                        + ((x[0] - 1.5 * par[1]) > 0) * par[5];  // quick fix on exting tails on the signal function
  if (x[0] < par[1])
    return pedestal;
  //double  signal = (-1)*par[0]*pow((x[0]-par[1]),par[2])*exp(-(x[0]-par[1])*par[3]);
  double signal =                                   //
      pow((x[0] - par[1]), par[2])                  //
      * (par[0] * exp(-(x[0] - par[1]) * par[3])    //
         + par[5] * exp(-(x[0] - par[1]) * par[6])  //
         );
  return pedestal + signal;
}
