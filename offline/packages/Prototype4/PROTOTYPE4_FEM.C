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
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

int PROTOTYPE4_FEM::GetHBDCh(std::string caloname, int i_column, int i_row)
{
  if (caloname == "HCALIN")
  {
    return 64 + 8 * i_column + 2 * i_row;
  }
  else if (caloname == "HCALOUT")
  {
    return 112 + 8 * i_column + 2 * i_row;
  }
  else if (caloname == "EMCAL")
  {
    // EMcal cable mapping from John haggerty
    assert(i_column >= 0);
    assert(i_column < NCH_EMCAL_COLUMNS);
    assert(i_row >= 0);
    assert(i_row < NCH_EMCAL_ROWS);

    static int canmap[] =
        {
            //            > https://docdb.sphenix.bnl.gov/0000/000034/001/T1044-2017a-2.xlsx
            //            Sean and John spot checked a number of these towers at BNL.  There could be mistakes, but cosmics look reasonable.
            //            Front view (looking downstream, same as above) but in HBD channel numbers
            //            3  , 2  , 19 , 18 , 35 , 34 , 51 , 50,
            //            1  , 0  , 17 , 16 , 33 , 32 , 49 , 48,
            //            7  , 6  , 23 , 22 , 39 , 38 , 55 , 54,
            //            5  , 4  , 21 , 20 , 37 , 36 , 53 , 52,
            //            11 , 10 , 27 , 26 , 43 , 42 , 59 , 58,
            //            9  , 8  , 25 , 24 , 41 , 40 , 57 , 56,
            //            15 , 14 , 31 , 30 , 47 , 46 , 63 , 62,
            //            13 , 12 , 29 , 28 , 45 , 44 , 61 , 60,
            3, 2, 19, 18, 35, 34, 51, 50,
            1, 0, 17, 16, 33, 32, 49, 48,
            7, 6, 23, 22, 39, 38, 55, 54,
            5, 4, 21, 20, 37, 36, 53, 52,
            11, 10, 27, 26, 43, 42, 59, 58,
            9, 8, 25, 24, 41, 40, 57, 56,
            15, 14, 31, 30, 47, 46, 63, 62,
            13, 12, 29, 28, 45, 44, 61, 60,
            0};

    const int tower_index = i_column + NCH_EMCAL_COLUMNS * (NCH_EMCAL_ROWS - 1 - i_row);

    assert(tower_index >= 0);
    assert(tower_index < NCH_EMCAL_ROWS * NCH_EMCAL_COLUMNS);

    return canmap[tower_index];
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
  fits.SetParNames("Amplitude", "Peak Sample", "Power", "Decay","Pedstal","Baseline shift");
  fits.SetParLimits(0, peakval * 0.9, peakval * 1.1);
  fits.SetParLimits(1, 0, NSAMPLES);
  fits.SetParLimits(2, 2, 4.);
  fits.SetParLimits(3, 1., 2.);
  fits.SetParLimits(4, pedestal - abs(peakval), pedestal + abs(peakval));
  //  fits.SetParLimits(5, - abs(peakval),  + abs(peakval));
//  fits.FixParameter(5, 0);

  //Saturation correction - Abhisek
  for (int ipoint = 0; ipoint < gpulse.GetN(); ipoint++)
    if ((gpulse.GetY())[ipoint] <=10  or (gpulse.GetY())[ipoint] >= ((1<<14) - 10))  // drop point if touching max or low limit on ADCs
    {
      gpulse.RemovePoint(ipoint);
      ipoint--;
    }

  gpulse.Fit(&fits, "MQRN0", "goff", 0., (double) NSAMPLES);

  if (verbosity)
  {
    static int id = 0;
    ++id;

    TCanvas *canvas = new TCanvas(
        (string("PROTOTYPE4_FEM_SampleFit_PowerLawExp") + to_string(id)).c_str(), "PROTOTYPE4_FEM::SampleFit_PowerLawExp");
    TGraph * g_plot = static_cast<TGraph *>(gpulse.DrawClone("ap*l"));
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
