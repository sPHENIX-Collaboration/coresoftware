#include "TPCDaqDefs.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVirtualFitter.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

using namespace std;

namespace TPCDaqDefs
{
//! TPC v1 FEE test stand decoder
namespace FEEv1
{
SampleFit_PowerLawDoubleExp_PDFMaker::SampleFit_PowerLawDoubleExp_PDFMaker()
{
  gStyle->SetOptFit(1111);

  m_canvas = new TCanvas("SampleFit_PowerLawDoubleExp_PDFMaker", "SampleFit_PowerLawDoubleExp_PDFMaker");
  m_pavedtext = new TPaveText(.05, .1, .95, .8);

  m_pavedtext->AddText("SampleFit_PowerLawDoubleExp Fit output");
  m_pavedtext->AddText("A double-component power-law exponential fit of time-dependent ADC pulses.");
  m_pavedtext->AddText("Magenta curve is the sum of the two component, the red and blue curves.");
  m_pavedtext->AddText("Red dot denote the max points");
  m_pavedtext->Draw();

  m_canvas->Print("SampleFit_PowerLawDoubleExp.pdf(");  //open multiplage PDF
}
SampleFit_PowerLawDoubleExp_PDFMaker::~SampleFit_PowerLawDoubleExp_PDFMaker()
{
  if (m_pavedtext) delete m_pavedtext;
  if (m_canvas) delete m_canvas;

  m_canvas = new TCanvas("SampleFit_PowerLawDoubleExp_PDFMaker", "SampleFit_PowerLawDoubleExp_PDFMaker");
  m_pavedtext = new TPaveText(.05, .1, .95, .8);

  m_pavedtext->AddText("SampleFit_PowerLawDoubleExp Fit output");
  m_pavedtext->AddText("End of pages");
  m_pavedtext->Draw();

  m_canvas->Print("SampleFit_PowerLawDoubleExp.pdf)");  //close multiplage PDF
}

void SampleFit_PowerLawDoubleExp_PDFMaker::MakeSectionPage(const string &title)
{
  if (m_pavedtext) delete m_pavedtext;
  if (m_canvas) delete m_canvas;

  m_canvas = new TCanvas("SampleFit_PowerLawDoubleExp_PDFMaker", "SampleFit_PowerLawDoubleExp_PDFMaker");

  m_pavedtext = new TPaveText(.05, .1, .95, .8);

  m_pavedtext->AddText(title.c_str());
  m_pavedtext->Draw();

  m_canvas->Print("SampleFit_PowerLawDoubleExp.pdf");
}

bool SampleFit_PowerLawDoubleExp(        //
    const std::vector<double> &samples,  //
    double &peak,                        //
    double &peak_sample,                 //
    double &pedestal,                    //
    std::map<int, double> &parameters_io,
    const int verbosity)
{
  static const int n_parameter = 7;

  // inital guesses
  int peakPos = 0.;

  //  assert(samples.size() == n_samples);
  const int n_samples = samples.size();

  TGraph gpulse(n_samples);
  for (int i = 0; i < n_samples; i++)
  {
    (gpulse.GetX())[i] = i;

    (gpulse.GetY())[i] = samples[i];
  }

  //Saturation correction - Abhisek
  //  for (int ipoint = 0; ipoint < gpulse.GetN(); ipoint++)
  //    if ((gpulse.GetY())[ipoint] >= ((1 << 10) - 10)  // drop point if touching max or low limit on ADCs
  //        or (not isnormal((gpulse.GetY())[ipoint])))
  //    {
  //      gpulse.RemovePoint(ipoint);
  //      ipoint--;
  //    }

  pedestal = gpulse.GetY()[0];  //(double) PEDESTAL;
  double peakval = pedestal;
  const double risetime = 1.5;

  for (int iSample = 0; iSample < n_samples - risetime * 3; iSample++)
  {
    if (abs(gpulse.GetY()[iSample] - pedestal) > abs(peakval - pedestal))
    {
      peakval = gpulse.GetY()[iSample];
      peakPos = iSample;
    }
  }
  peakval -= pedestal;

  if (verbosity)
  {
    cout << "SampleFit_PowerLawDoubleExp - "
         << "pedestal = " << pedestal << ", "
         << "peakval = " << peakval << ", "
         << "peakPos = " << peakPos << endl;
  }

  // build default value
  struct default_values_t
  {
    default_values_t(double default_value, double min_value, double max_value)
      : def(default_value)
      , min(min_value)
      , max(max_value)
    {
    }
    double def;
    double min;
    double max;
  };

  vector<default_values_t> default_values(n_parameter, default_values_t(numeric_limits<double>::signaling_NaN(), numeric_limits<double>::signaling_NaN(), numeric_limits<double>::signaling_NaN()));

  default_values[0] = default_values_t(peakval * .7, peakval * -1.5, peakval * 1.5);
  default_values[1] = default_values_t(peakPos - risetime, peakPos - 3 * risetime, peakPos + risetime);
  default_values[2] = default_values_t(5., 1, 10.);
  default_values[3] = default_values_t(risetime, risetime * .2, risetime * 10);
  default_values[4] = default_values_t(pedestal, pedestal - abs(peakval), pedestal + abs(peakval));
  //  default_values[5] = default_values_t(0.3, 0, 1);
  //  default_values[6] = default_values_t(5, risetime * .2, risetime * 10);
  default_values[5] = default_values_t(0, 0, 0);  // disable 2nd component
  default_values[6] = default_values_t(risetime, risetime, risetime);

  // fit function
  TF1 fits("f_SignalShape_PowerLawDoubleExp", SignalShape_PowerLawDoubleExp, 0., n_samples, n_parameter);
  fits.SetParNames("Amplitude", "Sample Start", "Power", "Peak Time 1", "Pedestal", "Amplitude ratio", "Peak Time 2");

  for (int i = 0; i < n_parameter; ++i)
  {
    if (parameters_io.find(i) == parameters_io.end())
    {
      fits.SetParameter(i, default_values[i].def);

      if (default_values[i].min < default_values[i].max)
      {
        fits.SetParLimits(i, default_values[i].min, default_values[i].max);
      }
      else
      {
        fits.FixParameter(i, default_values[i].def);
      }

      if (verbosity)
      {
        cout << "SampleFit_PowerLawDoubleExp - parameter [" << i << "]: "
             << "default value = " << default_values[i].def
             << ", min value = " << default_values[i].min
             << ", max value = " << default_values[i].max << endl;
      }
    }
    else
    {
//      fits.SetParLimits(i, parameters_io[i], parameters_io[i]);
      fits.SetParameter(i, parameters_io[i]);
      fits.FixParameter(i, parameters_io[i]);

      if (verbosity)
      {
        cout << "SampleFit_PowerLawDoubleExp - parameter [" << i << "]: fixed to " << parameters_io[i] << endl;
      }
    }
  }

  if (verbosity <= 1)
    gpulse.Fit(&fits, "QRN0W", "goff", 0., (double) n_samples);
  else
    gpulse.Fit(&fits, "RN0VW+", "goff", 0., (double) n_samples);

  // store results
  pedestal = fits.GetParameter(4);

  const double peakpos1 = fits.GetParameter(3);
  const double peakpos2 = fits.GetParameter(6);
  double max_peakpos = fits.GetParameter(1) + (peakpos1 > peakpos2 ? peakpos1 : peakpos2);
  if (max_peakpos > n_samples - 1) max_peakpos = n_samples - 1;

  if (fits.GetParameter(0) > 0)
    peak_sample = fits.GetMaximumX(fits.GetParameter(1), max_peakpos);
  else
    peak_sample = fits.GetMinimumX(fits.GetParameter(1), max_peakpos);

  peak = fits.Eval(peak_sample) - pedestal;

  if (verbosity)
  {
    static int id = 0;
    ++id;

    string c_name(string("SampleFit_PowerLawDoubleExp_") + to_string(id));

    TCanvas *canvas = new TCanvas(
        c_name.c_str(), c_name.c_str());
    canvas->Update();

    TGraph *g_plot = static_cast<TGraph *>(gpulse.DrawClone("ap*l"));
    g_plot->SetTitle((string("ADC data and fit #") + to_string(id) + string(";Sample number;ADC value")).c_str());

    fits.SetLineColor(kMagenta);
    fits.DrawClone("same");
    fits.Print();

    TF1 f1("f_SignalShape_PowerLawExp1", SignalShape_PowerLawExp, 0., n_samples, 5);
    f1.SetParameters(
        fits.GetParameter(0) * (1 - fits.GetParameter(5)) / pow(fits.GetParameter(3), fits.GetParameter(2)) * exp(fits.GetParameter(2)),
        fits.GetParameter(1),
        fits.GetParameter(2),
        fits.GetParameter(2) / fits.GetParameter(3),
        fits.GetParameter(4));
    f1.SetLineColor(kBlue);
    f1.DrawClone("same");

    TF1 f2("f_SignalShape_PowerLawExp2", SignalShape_PowerLawExp, 0., n_samples, 5);
    f2.SetParameters(
        fits.GetParameter(0) * fits.GetParameter(5) / pow(fits.GetParameter(6), fits.GetParameter(2)) * exp(fits.GetParameter(2)),
        fits.GetParameter(1),
        fits.GetParameter(2),
        fits.GetParameter(2) / fits.GetParameter(6),
        fits.GetParameter(4));
    f2.SetLineColor(kRed);
    f2.DrawClone("same");

    TGraph g_max(1);

    g_max.GetX()[0] = peak_sample;
    g_max.GetY()[0] = peak + pedestal;

    g_max.SetMarkerStyle(kFullCircle);
    g_max.SetMarkerSize(2);
    g_max.SetMarkerColor(kRed);

    static_cast<TGraph *>(g_max.DrawClone("p"));

    canvas->Update();

    //    if (id == 1)
    //    {
    //      canvas->Print("SampleFit_PowerLawDoubleExp.pdf(");
    //    }
    canvas->Print("SampleFit_PowerLawDoubleExp.pdf");
  }

  for (int i = 0; i < n_parameter; ++i)
  {
    parameters_io[i] = fits.GetParameter(i);
  }

  if (verbosity)
  {
    cout << "SampleFit_PowerLawDoubleExp - "
         << "peak_sample = " << peak_sample << ", "
         << "max_peakpos = " << max_peakpos << ", "
         << "fits.GetParameter(1) = " << fits.GetParameter(1) << ", "
         << "peak = " << peak << ", "
         << "pedestal = " << pedestal << endl;
  }

  return true;
}

double
SignalShape_PowerLawExp(double *x, double *par)
{
  double pedestal = par[4];
  //                        + ((x[0] - 1.5 * par[1]) > 0) * par[5];  // quick fix on exting tails on the signal function
  if (x[0] < par[1])
    return pedestal;
  //double  signal = (-1)*par[0]*pow((x[0]-par[1]),par[2])*exp(-(x[0]-par[1])*par[3]);
  double signal = par[0] * pow((x[0] - par[1]), par[2]) * exp(-(x[0] - par[1]) * par[3]);
  return pedestal + signal;
}

double
SignalShape_PowerLawDoubleExp(double *x, double *par)
{
  double pedestal = par[4];
  //                        + ((x[0] - 1.5 * par[1]) > 0) * par[5];  // quick fix on exting tails on the signal function
  if (x[0] < par[1])
    return pedestal;
  //double  signal = (-1)*par[0]*pow((x[0]-par[1]),par[2])*exp(-(x[0]-par[1])*par[3]);
  //  peak / pow(fits.GetParameter(2) / fits.GetParameter(3), fits.GetParameter(2)) * exp(fits.GetParameter(2)) = fits.GetParameter(0);  // exact peak height is (p0*Power(p2/p3,p2))/Power(E,p2)
  //  fits.GetParameter(2) / peak_shift =  fits.GetParameter(3);  // signal peak time

  double signal =                                                                                         //
      par[0]                                                                                              //
      * pow((x[0] - par[1]), par[2])                                                                      //
      * (((1. - par[5]) / pow(par[3], par[2]) * exp(par[2])) * exp(-(x[0] - par[1]) * (par[2] / par[3]))  //
         + (par[5] / pow(par[6], par[2]) * exp(par[2])) * exp(-(x[0] - par[1]) * (par[2] / par[6]))       //
        );
  return pedestal + signal;
}

pair<int, int> SAMPAChan2PadXY(uint32_t fee_channel)
{
  static const int pad_map_to_xy[512][2] = {
      {-1, 0},   //0
      {-1, 1},   //1
      {-1, 2},   //2
      {1, 2},    //3
      {1, 1},    //4
      {1, 4},    //5
      {1, 3},    //6
      {2, 1},    //7
      {1, 5},    //8
      {2, 3},    //9
      {2, 2},    //10
      {2, 5},    //11
      {2, 4},    //12
      {3, 2},    //13
      {3, 1},    //14
      {3, 4},    //15
      {3, 3},    //16
      {4, 1},    //17
      {3, 5},    //18
      {4, 3},    //19
      {4, 2},    //20
      {4, 5},    //21
      {4, 4},    //22
      {5, 2},    //23
      {5, 1},    //24
      {5, 4},    //25
      {5, 3},    //26
      {6, 1},    //27
      {5, 5},    //28
      {6, 3},    //29
      {6, 2},    //30
      {6, 5},    //31
      {6, 4},    //32
      {7, 2},    //33
      {7, 1},    //34
      {7, 4},    //35
      {7, 3},    //36
      {8, 1},    //37
      {7, 5},    //38
      {8, 3},    //39
      {8, 2},    //40
      {8, 5},    //41
      {8, 4},    //42
      {9, 2},    //43
      {9, 1},    //44
      {9, 4},    //45
      {9, 3},    //46
      {10, 1},   //47
      {9, 5},    //48
      {10, 3},   //49
      {10, 2},   //50
      {10, 5},   //51
      {10, 4},   //52
      {11, 2},   //53
      {11, 1},   //54
      {11, 4},   //55
      {11, 3},   //56
      {12, 1},   //57
      {11, 5},   //58
      {12, 3},   //59
      {12, 2},   //60
      {12, 5},   //61
      {12, 4},   //62
      {13, 2},   //63
      {13, 1},   //64
      {13, 4},   //65
      {13, 3},   //66
      {14, 1},   //67
      {13, 5},   //68
      {14, 3},   //69
      {14, 2},   //70
      {14, 5},   //71
      {14, 4},   //72
      {15, 2},   //73
      {15, 1},   //74
      {15, 4},   //75
      {15, 3},   //76
      {16, 1},   //77
      {15, 5},   //78
      {16, 3},   //79
      {16, 2},   //80
      {16, 5},   //81
      {16, 4},   //82
      {17, 2},   //83
      {17, 1},   //84
      {17, 4},   //85
      {17, 3},   //86
      {18, 1},   //87
      {17, 5},   //88
      {18, 3},   //89
      {18, 2},   //90
      {18, 5},   //91
      {18, 4},   //92
      {19, 2},   //93
      {19, 1},   //94
      {19, 4},   //95
      {19, 3},   //96
      {20, 1},   //97
      {19, 5},   //98
      {20, 3},   //99
      {20, 2},   //100
      {20, 5},   //101
      {20, 4},   //102
      {21, 2},   //103
      {21, 1},   //104
      {21, 4},   //105
      {21, 3},   //106
      {22, 1},   //107
      {21, 5},   //108
      {22, 3},   //109
      {22, 2},   //110
      {22, 5},   //111
      {22, 4},   //112
      {23, 2},   //113
      {23, 1},   //114
      {23, 4},   //115
      {23, 3},   //116
      {24, 1},   //117
      {23, 5},   //118
      {24, 3},   //119
      {24, 2},   //120
      {24, 5},   //121
      {24, 4},   //122
      {-1, 3},   //123
      {-1, 4},   //124
      {-1, 5},   //125
      {-1, 6},   //126
      {-1, 7},   //127
      {-9, 2},   //128
      {-3, 1},   //129
      {-3, 2},   //130
      {25, 2},   //131
      {25, 1},   //132
      {25, 4},   //133
      {25, 3},   //134
      {26, 1},   //135
      {25, 5},   //136
      {26, 3},   //137
      {26, 2},   //138
      {26, 5},   //139
      {26, 4},   //140
      {27, 2},   //141
      {27, 1},   //142
      {27, 4},   //143
      {27, 3},   //144
      {28, 1},   //145
      {27, 5},   //146
      {28, 3},   //147
      {28, 2},   //148
      {28, 5},   //149
      {28, 4},   //150
      {29, 2},   //151
      {29, 1},   //152
      {29, 4},   //153
      {29, 3},   //154
      {30, 1},   //155
      {29, 5},   //156
      {30, 3},   //157
      {30, 2},   //158
      {30, 5},   //159
      {30, 4},   //160
      {31, 2},   //161
      {31, 1},   //162
      {31, 4},   //163
      {31, 3},   //164
      {32, 1},   //165
      {31, 5},   //166
      {32, 3},   //167
      {32, 2},   //168
      {32, 5},   //169
      {32, 4},   //170
      {33, 2},   //171
      {33, 1},   //172
      {33, 4},   //173
      {33, 3},   //174
      {34, 1},   //175
      {33, 5},   //176
      {34, 3},   //177
      {34, 2},   //178
      {34, 5},   //179
      {34, 4},   //180
      {35, 2},   //181
      {35, 1},   //182
      {35, 4},   //183
      {35, 3},   //184
      {36, 1},   //185
      {35, 5},   //186
      {36, 3},   //187
      {36, 2},   //188
      {36, 5},   //189
      {36, 4},   //190
      {37, 2},   //191
      {37, 1},   //192
      {37, 4},   //193
      {37, 3},   //194
      {38, 1},   //195
      {37, 5},   //196
      {38, 3},   //197
      {38, 2},   //198
      {38, 5},   //199
      {38, 4},   //200
      {39, 2},   //201
      {39, 1},   //202
      {39, 4},   //203
      {39, 3},   //204
      {40, 1},   //205
      {39, 5},   //206
      {40, 3},   //207
      {40, 2},   //208
      {40, 5},   //209
      {40, 4},   //210
      {41, 2},   //211
      {41, 1},   //212
      {41, 4},   //213
      {41, 3},   //214
      {42, 1},   //215
      {41, 5},   //216
      {42, 3},   //217
      {42, 2},   //218
      {42, 5},   //219
      {42, 4},   //220
      {43, 2},   //221
      {43, 1},   //222
      {43, 4},   //223
      {43, 3},   //224
      {44, 1},   //225
      {43, 5},   //226
      {44, 3},   //227
      {44, 2},   //228
      {44, 5},   //229
      {44, 4},   //230
      {45, 2},   //231
      {45, 1},   //232
      {45, 4},   //233
      {45, 3},   //234
      {46, 1},   //235
      {45, 5},   //236
      {46, 3},   //237
      {46, 2},   //238
      {46, 5},   //239
      {46, 4},   //240
      {47, 2},   //241
      {47, 1},   //242
      {47, 4},   //243
      {47, 3},   //244
      {48, 1},   //245
      {47, 5},   //246
      {48, 3},   //247
      {48, 2},   //248
      {48, 5},   //249
      {48, 4},   //250
      {-3, 3},   //251
      {-3, 4},   //252
      {-3, 5},   //253
      {-3, 6},   //254
      {-3, 7},   //255
      {-5, 0},   //256
      {-5, 1},   //257
      {-5, 2},   //258
      {48, 6},   //259
      {-5, 3},   //260
      {48, 8},   //261
      {48, 7},   //262
      {48, 10},  //263
      {48, 9},   //264
      {47, 7},   //265
      {47, 6},   //266
      {47, 9},   //267
      {47, 8},   //268
      {46, 6},   //269
      {47, 10},  //270
      {46, 8},   //271
      {46, 7},   //272
      {46, 10},  //273
      {46, 9},   //274
      {45, 7},   //275
      {45, 6},   //276
      {45, 9},   //277
      {45, 8},   //278
      {44, 6},   //279
      {45, 10},  //280
      {44, 8},   //281
      {44, 7},   //282
      {44, 10},  //283
      {44, 9},   //284
      {43, 7},   //285
      {43, 6},   //286
      {43, 9},   //287
      {43, 8},   //288
      {42, 6},   //289
      {43, 10},  //290
      {42, 8},   //291
      {42, 7},   //292
      {42, 10},  //293
      {42, 9},   //294
      {41, 7},   //295
      {41, 6},   //296
      {41, 9},   //297
      {41, 8},   //298
      {40, 6},   //299
      {41, 10},  //300
      {40, 8},   //301
      {40, 7},   //302
      {40, 10},  //303
      {40, 9},   //304
      {39, 7},   //305
      {39, 6},   //306
      {39, 9},   //307
      {39, 8},   //308
      {38, 6},   //309
      {39, 10},  //310
      {38, 8},   //311
      {38, 7},   //312
      {38, 10},  //313
      {38, 9},   //314
      {37, 7},   //315
      {37, 6},   //316
      {37, 9},   //317
      {37, 8},   //318
      {36, 6},   //319
      {37, 10},  //320
      {36, 8},   //321
      {36, 7},   //322
      {36, 10},  //323
      {36, 9},   //324
      {35, 7},   //325
      {35, 6},   //326
      {35, 9},   //327
      {35, 8},   //328
      {34, 6},   //329
      {35, 10},  //330
      {34, 8},   //331
      {34, 7},   //332
      {34, 10},  //333
      {34, 9},   //334
      {33, 7},   //335
      {33, 6},   //336
      {33, 9},   //337
      {33, 8},   //338
      {32, 6},   //339
      {33, 10},  //340
      {32, 8},   //341
      {32, 7},   //342
      {32, 10},  //343
      {32, 9},   //344
      {31, 7},   //345
      {31, 6},   //346
      {31, 9},   //347
      {31, 8},   //348
      {30, 6},   //349
      {31, 10},  //350
      {30, 8},   //351
      {30, 7},   //352
      {30, 10},  //353
      {30, 9},   //354
      {29, 7},   //355
      {29, 6},   //356
      {29, 9},   //357
      {29, 8},   //358
      {28, 6},   //359
      {29, 10},  //360
      {28, 8},   //361
      {28, 7},   //362
      {28, 10},  //363
      {28, 9},   //364
      {27, 7},   //365
      {27, 6},   //366
      {27, 9},   //367
      {27, 8},   //368
      {26, 6},   //369
      {27, 10},  //370
      {26, 8},   //371
      {26, 7},   //372
      {26, 10},  //373
      {26, 9},   //374
      {25, 7},   //375
      {25, 6},   //376
      {25, 9},   //377
      {25, 8},   //378
      {-9, 0},   //379
      {25, 10},  //380
      {-5, 5},   //381
      {-5, 6},   //382
      {-5, 7},   //383
      {-7, 0},   //384
      {-7, 1},   //385
      {-7, 2},   //386
      {24, 6},   //387
      {-7, 3},   //388
      {24, 8},   //389
      {24, 7},   //390
      {24, 10},  //391
      {24, 9},   //392
      {23, 7},   //393
      {23, 6},   //394
      {23, 9},   //395
      {23, 8},   //396
      {22, 6},   //397
      {23, 10},  //398
      {22, 8},   //399
      {22, 7},   //400
      {22, 10},  //401
      {22, 9},   //402
      {21, 7},   //403
      {21, 6},   //404
      {21, 9},   //405
      {21, 8},   //406
      {20, 6},   //407
      {21, 10},  //408
      {20, 8},   //409
      {20, 7},   //410
      {20, 10},  //411
      {20, 9},   //412
      {19, 7},   //413
      {19, 6},   //414
      {19, 9},   //415
      {19, 8},   //416
      {18, 6},   //417
      {19, 10},  //418
      {18, 8},   //419
      {18, 7},   //420
      {18, 10},  //421
      {18, 9},   //422
      {17, 7},   //423
      {17, 6},   //424
      {17, 9},   //425
      {17, 8},   //426
      {16, 6},   //427
      {17, 10},  //428
      {16, 8},   //429
      {16, 7},   //430
      {16, 10},  //431
      {16, 9},   //432
      {15, 7},   //433
      {15, 6},   //434
      {15, 9},   //435
      {15, 8},   //436
      {14, 6},   //437
      {15, 10},  //438
      {14, 8},   //439
      {14, 7},   //440
      {14, 10},  //441
      {14, 9},   //442
      {13, 7},   //443
      {13, 6},   //444
      {13, 9},   //445
      {13, 8},   //446
      {12, 6},   //447
      {13, 10},  //448
      {12, 8},   //449
      {12, 7},   //450
      {12, 10},  //451
      {12, 9},   //452
      {11, 7},   //453
      {11, 6},   //454
      {11, 9},   //455
      {11, 8},   //456
      {10, 6},   //457
      {11, 10},  //458
      {10, 8},   //459
      {10, 7},   //460
      {10, 10},  //461
      {10, 9},   //462
      {9, 7},    //463
      {9, 6},    //464
      {9, 9},    //465
      {9, 8},    //466
      {8, 6},    //467
      {9, 10},   //468
      {8, 8},    //469
      {8, 7},    //470
      {8, 10},   //471
      {8, 9},    //472
      {7, 7},    //473
      {7, 6},    //474
      {7, 9},    //475
      {7, 8},    //476
      {6, 6},    //477
      {7, 10},   //478
      {6, 8},    //479
      {6, 7},    //480
      {6, 10},   //481
      {6, 9},    //482
      {5, 7},    //483
      {5, 6},    //484
      {5, 9},    //485
      {5, 8},    //486
      {4, 6},    //487
      {5, 10},   //488
      {4, 8},    //489
      {4, 7},    //490
      {4, 10},   //491
      {4, 9},    //492
      {3, 7},    //493
      {3, 6},    //494
      {3, 9},    //495
      {3, 8},    //496
      {2, 6},    //497
      {3, 10},   //498
      {2, 8},    //499
      {2, 7},    //500
      {2, 10},   //501
      {2, 9},    //502
      {1, 7},    //503
      {1, 6},    //504
      {1, 9},    //505
      {1, 8},    //506
      {0, 6},    //507
      {1, 10},   //508
      {-7, 5},   //509
      {-7, 6},   //510
      {-7, 7}};

  static const int sampa_chan_to_pad[] = {
      287,  //0
      286,  //1
      285,  //2
      284,  //3
      283,  //4
      282,  //5
      281,  //6
      280,  //7
      279,  //8
      278,  //9
      277,  //10
      276,  //11
      275,  //12
      274,  //13
      273,  //14
      272,  //15
      271,  //16
      270,  //17
      269,  //18
      268,  //19
      267,  //20
      266,  //21
      265,  //22
      264,  //23
      263,  //24
      262,  //25
      261,  //26
      260,  //27
      259,  //28
      258,  //29
      257,  //30
      256,  //31
      319,  //32
      318,  //33
      317,  //34
      316,  //35
      315,  //36
      314,  //37
      313,  //38
      312,  //39
      311,  //40
      310,  //41
      309,  //42
      308,  //43
      307,  //44
      306,  //45
      305,  //46
      304,  //47
      303,  //48
      302,  //49
      301,  //50
      300,  //51
      299,  //52
      298,  //53
      297,  //54
      296,  //55
      295,  //56
      294,  //57
      293,  //58
      292,  //59
      291,  //60
      290,  //61
      289,  //62
      288,  //63
      351,  //64
      350,  //65
      349,  //66
      348,  //67
      347,  //68
      346,  //69
      345,  //70
      344,  //71
      343,  //72
      342,  //73
      341,  //74
      340,  //75
      339,  //76
      338,  //77
      337,  //78
      336,  //79
      335,  //80
      334,  //81
      333,  //82
      332,  //83
      331,  //84
      330,  //85
      329,  //86
      328,  //87
      327,  //88
      326,  //89
      325,  //90
      324,  //91
      323,  //92
      322,  //93
      321,  //94
      320,  //95
      383,  //96
      382,  //97
      381,  //98
      380,  //99
      379,  //100
      378,  //101
      377,  //102
      376,  //103
      375,  //104
      374,  //105
      373,  //106
      372,  //107
      371,  //108
      370,  //109
      369,  //110
      368,  //111
      367,  //112
      366,  //113
      365,  //114
      364,  //115
      363,  //116
      362,  //117
      361,  //118
      360,  //119
      359,  //120
      358,  //121
      357,  //122
      356,  //123
      355,  //124
      354,  //125
      353,  //126
      352,  //127
      415,  //128
      414,  //129
      413,  //130
      412,  //131
      411,  //132
      410,  //133
      409,  //134
      408,  //135
      407,  //136
      406,  //137
      405,  //138
      404,  //139
      403,  //140
      402,  //141
      401,  //142
      400,  //143
      399,  //144
      398,  //145
      397,  //146
      396,  //147
      395,  //148
      394,  //149
      393,  //150
      392,  //151
      391,  //152
      390,  //153
      389,  //154
      388,  //155
      387,  //156
      386,  //157
      385,  //158
      384,  //159
      447,  //160
      446,  //161
      445,  //162
      444,  //163
      443,  //164
      442,  //165
      441,  //166
      440,  //167
      439,  //168
      438,  //169
      437,  //170
      436,  //171
      435,  //172
      434,  //173
      433,  //174
      432,  //175
      431,  //176
      430,  //177
      429,  //178
      428,  //179
      427,  //180
      426,  //181
      425,  //182
      424,  //183
      423,  //184
      422,  //185
      421,  //186
      420,  //187
      419,  //188
      418,  //189
      417,  //190
      416,  //191
      479,  //192
      478,  //193
      477,  //194
      476,  //195
      475,  //196
      474,  //197
      473,  //198
      472,  //199
      471,  //200
      470,  //201
      469,  //202
      468,  //203
      467,  //204
      466,  //205
      465,  //206
      464,  //207
      463,  //208
      462,  //209
      461,  //210
      460,  //211
      459,  //212
      458,  //213
      457,  //214
      456,  //215
      455,  //216
      454,  //217
      453,  //218
      452,  //219
      451,  //220
      450,  //221
      449,  //222
      448,  //223
      511,  //224
      510,  //225
      509,  //226
      508,  //227
      507,  //228
      506,  //229
      505,  //230
      504,  //231
      503,  //232
      502,  //233
      501,  //234
      500,  //235
      499,  //236
      498,  //237
      497,  //238
      496,  //239
      495,  //240
      494,  //241
      493,  //242
      492,  //243
      491,  //244
      490,  //245
      489,  //246
      488,  //247
      487,  //248
      486,  //249
      485,  //250
      484,  //251
      483,  //252
      482,  //253
      481,  //254
      480,  //255
      0     //256 , the invalid ID
  };

  if (fee_channel >= 256)
    fee_channel = 256;

  const int pad_number = sampa_chan_to_pad[fee_channel];
  const int pad_x = pad_map_to_xy[pad_number][0];
  const int pad_y = pad_map_to_xy[pad_number][1];

  return make_pair(pad_x, pad_y);
}

}  // namespace FEEv1

}  // namespace TPCDaqDefs
