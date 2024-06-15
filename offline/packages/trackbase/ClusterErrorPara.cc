#include "ClusterErrorPara.h"
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv4.h>

#include <TF1.h>
#include <cmath>
#include <iostream>
#include <string>

namespace
{

  //! convenience square method
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }
}  // namespace

ClusterErrorPara::ClusterErrorPara()
{
  f0 = new TF1("f0", "pol1", 0, 10);
  f0->SetParameter(0, 0.0163943);
  f0->SetParameter(1, 0.0192931);

  f1 = new TF1("f1", "pol2", 0, 10);
  f1->SetParameter(0, 0.0119384);
  f1->SetParameter(1, 0.0253197);
  f1->SetParameter(2, 0.0404213);

  f2 = new TF1("f2", "pol2", 0, 10);
  f2->SetParameter(0, 0.0107316);
  f2->SetParameter(1, 0.0294968);
  f2->SetParameter(2, 0.0414098);
  // f2->SetParameter(3,9.75877);

  fz0 = new TF1("fz0", "pol2", -2, 2);
  fz0->SetParameter(0, 0.0520278);
  fz0->SetParameter(1, -0.00578699);
  fz0->SetParameter(2, 0.0156972);

  fz1 = new TF1("fz1", "pol4", -2, 2);
  fz1->SetParameter(0, 0.0383233);
  fz1->SetParameter(1, -0.00577128);
  fz1->SetParameter(2, 0.0770914);
  fz1->SetParameter(3, -0.0818139);
  fz1->SetParameter(4, 0.050305);

  fz2 = new TF1("fz2", "pol2", -2, 2);
  fz2->SetParameter(0, 0.0371611);
  fz2->SetParameter(1, -0.000694558);
  fz2->SetParameter(2, 0.0437917);

  fmm_55_2 = new TF1("fmm_55_2", "pol2", -2, 2);
  fmm_55_2->SetParameter(0, 0.0430592);
  fmm_55_2->SetParameter(1, -0.000177174);
  fmm_55_2->SetParameter(2, 0.0914288);

  fmm_56_2 = new TF1("fmm_56_2", "pol2", -2, 2);
  fmm_56_2->SetParameter(0, 0.00363897);
  fmm_56_2->SetParameter(1, 0.0109713);
  fmm_56_2->SetParameter(2, 0.032354);

  fmm_3 = new TF1("fmm_3", "pol2", -2, 2);
  fmm_3->SetParameter(0, 0.00305396);
  fmm_3->SetParameter(1, 0.00505814);
  fmm_3->SetParameter(2, 0.0395137);

  fadcz0 = new TF1("fadcz0", "pol5", 0, 20000);
  fadcz0->SetParameter(0, 2.08854);
  fadcz0->SetParameter(1, -0.0536847);
  fadcz0->SetParameter(2, 0.000989393);
  fadcz0->SetParameter(3, -9.54492e-06);
  fadcz0->SetParameter(4, 4.42178e-08);
  fadcz0->SetParameter(5, -7.79669e-11);

  fadcz1 = new TF1("fadcz1", "pol5", 0, 20000);
  fadcz1->SetParameter(0, 2.35278);
  fadcz1->SetParameter(1, -0.0535903);
  fadcz1->SetParameter(2, 0.00088052);
  fadcz1->SetParameter(3, -7.75203e-06);
  fadcz1->SetParameter(4, 3.35361e-08);
  fadcz1->SetParameter(5, -5.61371e-11);

  fadcz2 = new TF1("fadcz2", "pol5", 0, 20000);
  fadcz2->SetParameter(0, 2.53191);
  fadcz2->SetParameter(1, -0.062285);
  fadcz2->SetParameter(2, 0.00103893);
  fadcz2->SetParameter(3, -9.18354e-06);
  fadcz2->SetParameter(4, 3.9802e-08);
  fadcz2->SetParameter(5, -6.67137e-11);

  fadcz0fine = new TF1("fadcz0fine", "[0]+([1]/pow(x-[2],2))", 0, 20000);
  fadcz0fine->SetParameter(0, 9.63983e-01);
  fadcz0fine->SetParameter(1, 2.68585e+01);
  fadcz0fine->SetParameter(2, -4.78664e+00);

  fadcz1fine = new TF1("fadcz1fine", "[0]+([1]/pow(x-[2],2))", 0, 20000);
  fadcz1fine->SetParameter(0, 9.85546e-01);
  fadcz1fine->SetParameter(1, 1.12622e+02);
  fadcz1fine->SetParameter(2, -1.26552e+01);

  fadcz2fine = new TF1("fadcz2fine", "[0]+([1]/pow(x-[2],2))", 0, 20000);
  fadcz2fine->SetParameter(0, 9.71125e-01);
  fadcz2fine->SetParameter(1, 6.67244e+01);
  fadcz2fine->SetParameter(2, -3.55034e+00);

  fadcphi0 = new TF1("fadcphi0", "pol4", 0, 20000);
  fadcphi0->SetParameter(0, 1.79273);
  fadcphi0->SetParameter(1, -0.0306044);
  fadcphi0->SetParameter(2, 0.000355984);
  fadcphi0->SetParameter(3, -2.09e-06);
  fadcphi0->SetParameter(4, 4.26161e-09);
  //  fadcphi0->SetParameter(5,-4.22758e-11);

  fadcphi0fine = new TF1("fadcphi0fine", "pol2", 0, 20000);
  fadcphi0fine->SetParameter(0, 1.02625);
  fadcphi0fine->SetParameter(1, -0.00167294);
  fadcphi0fine->SetParameter(2, 2.2912e-5);

  fadcphi1 = new TF1("fadcphi1", "pol4", 0, 20000);
  fadcphi1->SetParameter(0, 2.12873);
  fadcphi1->SetParameter(1, -0.0369604);
  fadcphi1->SetParameter(2, 0.00042828);
  fadcphi1->SetParameter(3, -2.3665e-06);
  fadcphi1->SetParameter(4, 4.87683e-09);

  fadcphi1fine = new TF1("fadcphi1fine", "pol4", 0, 20000);
  fadcphi1fine->SetParameter(0, 1.11749);
  fadcphi1fine->SetParameter(1, -0.00354277);
  fadcphi1fine->SetParameter(2, 5.60236e-05);
  fadcphi1fine->SetParameter(3, -4.46412e-07);
  fadcphi1fine->SetParameter(4, 1.22689e-09);

  fadcphi2 = new TF1("fadcphi2", "pol5", 0, 20000);
  fadcphi2->SetParameter(0, 2.29);
  fadcphi2->SetParameter(1, -0.0474362);
  fadcphi2->SetParameter(2, 0.000717789);
  fadcphi2->SetParameter(3, -6.00737e-06);
  fadcphi2->SetParameter(4, 2.52007e-08);
  fadcphi2->SetParameter(5, -4.14747e-11);

  fadcphi2fine1 = new TF1("fadcphi2fine1", "pol4", 0, 20000);
  fadcphi2fine1->SetParameter(0, 1.39404);
  fadcphi2fine1->SetParameter(1, -0.0202245);
  fadcphi2fine1->SetParameter(2, 0.000394666);
  fadcphi2fine1->SetParameter(3, -3.37831e-06);
  fadcphi2fine1->SetParameter(4, 1.05017e-08);

  fadcphi2fine2 = new TF1("fadcphi2fine2", "pol1", 0, 20000);
  fadcphi2fine2->SetParameter(0, 0.997);
  fadcphi2fine2->SetParameter(1, 0.00047);

  f0fine = new TF1("f0fine", "pol2", 0, 20000);
  f0fine->SetParameter(0, 0.98611);
  f0fine->SetParameter(1, -0.169505);
  f0fine->SetParameter(2, 1.12907);

  f1fine = new TF1("f1fine", "pol3", 0, 20000);
  f1fine->SetParameter(0, 0.968625);
  f1fine->SetParameter(1, -0.38894);
  f1fine->SetParameter(2, 3.36493);
  f1fine->SetParameter(3, -6.72275);
  /*
  f2fine = new TF1("f2fine","pol4",0,20000);
  f2fine->SetParameter(0,1.23748);
  f2fine->SetParameter(1,-2.56956);
  f2fine->SetParameter(2,15.8147);
  f2fine->SetParameter(3,-42.4668);
  f2fine->SetParameter(4,43.6083);
  */
  f2fine = new TF1("f2fine", "pol5", 0, 20000);
  f2fine->SetLineColor(kBlue);
  f2fine->SetParameter(0, 1.14119);
  f2fine->SetParameter(1, -2.81483);
  f2fine->SetParameter(2, 19.1877);
  f2fine->SetParameter(3, -57.214);
  f2fine->SetParameter(4, 72.2359);
  f2fine->SetParameter(5, -20.3802);

  fz0fine = new TF1("fz0fine", "pol2", 0, 20000);
  fz0fine->SetParameter(0, 0.96933);
  fz0fine->SetParameter(1, -0.0458534);
  fz0fine->SetParameter(2, 0.231419);

  fz1fine = new TF1("fz1fine", "pol3", 0, 20000);
  fz1fine->SetParameter(0, 0.886262);
  fz1fine->SetParameter(1, -0.0818167);
  fz1fine->SetParameter(2, 0.805824);
  fz1fine->SetParameter(3, -0.425423);

  fz2fine = new TF1("fz2fine", "pol5", 0, 20000);
  fz2fine->SetLineColor(kBlue);
  fz2fine->SetParameter(0, 0.880153);
  fz2fine->SetParameter(1, 0.552461);
  fz2fine->SetParameter(2, -2.57007);
  fz2fine->SetParameter(3, 7.509);
  fz2fine->SetParameter(4, -9.23698);
  fz2fine->SetParameter(5, 4.23039);

  static const double invsqrt12 = 1. / std::sqrt(12);

  pitcherr_phi_mvtx = 0.002688 * invsqrt12;
  pitcherr_phi_intt = 0.0078 * invsqrt12;
  pitcherr_phi_mm1 = 0.1 * invsqrt12;
  pitcherr_phi_mm2 = 31.6 * invsqrt12;
  pitcherr_z_mvtx = 0.002924 * invsqrt12;
  pitcherr_z_intt = 1.6 * invsqrt12;
  pitcherr_z_mm1 = 54.2 * invsqrt12;
  pitcherr_z_mm2 = 0.2 * invsqrt12;
  pull_fine_phi[0] = 1.055549;
  pull_fine_phi[1] = 1.049109;
  pull_fine_phi[2] = 1.043427;
  pull_fine_phi[3] = 1.200730;
  pull_fine_phi[4] = 1.157397;
  pull_fine_phi[5] = 1.357292;
  pull_fine_phi[6] = 1.367113;
  pull_fine_phi[7] = 0.999609;
  pull_fine_phi[8] = 1.020689;
  pull_fine_phi[9] = 1.032753;
  pull_fine_phi[10] = 1.018574;
  pull_fine_phi[11] = 1.035928;
  pull_fine_phi[12] = 1.024522;
  pull_fine_phi[13] = 1.035439;
  pull_fine_phi[14] = 1.020409;
  pull_fine_phi[15] = 1.035957;
  pull_fine_phi[16] = 1.033798;
  pull_fine_phi[17] = 1.032753;
  pull_fine_phi[18] = 1.030573;
  pull_fine_phi[19] = 1.035131;
  pull_fine_phi[20] = 1.031460;
  pull_fine_phi[21] = 1.039708;
  pull_fine_phi[22] = 1.015784;
  pull_fine_phi[23] = 1.013265;
  pull_fine_phi[24] = 0.997548;
  pull_fine_phi[25] = 0.990299;
  pull_fine_phi[26] = 0.999642;
  pull_fine_phi[27] = 1.005518;
  pull_fine_phi[28] = 1.012272;
  pull_fine_phi[29] = 1.016701;
  pull_fine_phi[30] = 1.013685;
  pull_fine_phi[31] = 1.025059;
  pull_fine_phi[32] = 1.030672;
  pull_fine_phi[33] = 1.031959;
  pull_fine_phi[34] = 1.040878;
  pull_fine_phi[35] = 1.045917;
  pull_fine_phi[36] = 1.052541;
  pull_fine_phi[37] = 1.056798;
  pull_fine_phi[38] = 1.000142;
  pull_fine_phi[39] = 1.028008;
  pull_fine_phi[40] = 0.980509;
  pull_fine_phi[41] = 0.981643;
  pull_fine_phi[42] = 0.988486;
  pull_fine_phi[43] = 0.987750;
  pull_fine_phi[44] = 0.989536;
  pull_fine_phi[45] = 0.992817;
  pull_fine_phi[46] = 0.993260;
  pull_fine_phi[47] = 0.994378;
  pull_fine_phi[48] = 0.988871;
  pull_fine_phi[49] = 0.985053;
  pull_fine_phi[50] = 0.985631;
  pull_fine_phi[51] = 0.987169;
  pull_fine_phi[52] = 0.992417;
  pull_fine_phi[53] = 0.996130;
  pull_fine_phi[54] = 0.996966;
  pull_fine_phi[55] = 0.943431;
  pull_fine_phi[56] = 0.000000;
  pull_fine_phi[57] = 0.000000;
  pull_fine_phi[58] = 0.000000;
  pull_fine_phi[59] = 0.000000;

  pull_fine_phi[3] *= 1.007551;
  pull_fine_phi[4] *= 1.006760;
  pull_fine_phi[5] *= 1.026019;
  pull_fine_phi[6] *= 1.030869;
  pull_fine_phi[55] *= 0.988674;

  pull_fine_z[0] = 1.117666;
  pull_fine_z[1] = 1.119458;
  pull_fine_z[2] = 1.123506;
  pull_fine_z[3] = 1.971984;
  pull_fine_z[4] = 0.643933;
  pull_fine_z[5] = 0.638478;
  pull_fine_z[6] = 2.039448;
  pull_fine_z[7] = 1.047709;
  pull_fine_z[8] = 1.021305;
  pull_fine_z[9] = 1.014697;
  pull_fine_z[10] = 1.011649;
  pull_fine_z[11] = 1.014099;
  pull_fine_z[12] = 1.018134;
  pull_fine_z[13] = 1.004908;
  pull_fine_z[14] = 1.016988;
  pull_fine_z[15] = 1.007906;
  pull_fine_z[16] = 1.006559;
  pull_fine_z[17] = 1.006695;
  pull_fine_z[18] = 1.010511;
  pull_fine_z[19] = 1.008526;
  pull_fine_z[20] = 1.012481;
  pull_fine_z[21] = 1.005427;
  pull_fine_z[22] = 1.008053;
  pull_fine_z[23] = 0.984082;
  pull_fine_z[24] = 1.016605;
  pull_fine_z[25] = 1.020224;
  pull_fine_z[26] = 1.014571;
  pull_fine_z[27] = 1.017648;
  pull_fine_z[28] = 1.021230;
  pull_fine_z[29] = 1.022464;
  pull_fine_z[30] = 1.026075;
  pull_fine_z[31] = 1.028448;
  pull_fine_z[32] = 1.025604;
  pull_fine_z[33] = 1.025361;
  pull_fine_z[34] = 1.031499;
  pull_fine_z[35] = 1.036946;
  pull_fine_z[36] = 1.037130;
  pull_fine_z[37] = 1.036079;
  pull_fine_z[38] = 0.966569;
  pull_fine_z[39] = 3.431897;
  pull_fine_z[40] = 1.000976;
  pull_fine_z[41] = 1.000437;
  pull_fine_z[42] = 0.999851;
  pull_fine_z[43] = 1.003534;
  pull_fine_z[44] = 1.002914;
  pull_fine_z[45] = 1.003099;
  pull_fine_z[46] = 1.008661;
  pull_fine_z[47] = 1.015236;
  pull_fine_z[48] = 1.009287;
  pull_fine_z[49] = 1.011966;
  pull_fine_z[50] = 1.005105;
  pull_fine_z[51] = 1.014257;
  pull_fine_z[52] = 1.020858;
  pull_fine_z[53] = 1.018221;
  pull_fine_z[54] = 1.012372;
  pull_fine_z[55] = 1.007561;
  pull_fine_z[56] = 0.000000;
  pull_fine_z[57] = 0.000000;
  pull_fine_z[58] = 0.000000;
  pull_fine_z[59] = 0.000000;

  pull_fine_z[3] *= 0.691672;
  pull_fine_z[4] *= 1.203607;
  pull_fine_z[5] *= 1.360947;
  pull_fine_z[6] *= 0.583844;
  pull_fine_z[39] *= 0.928246;
  pull_fine_z[3] *= 0.498090;
  pull_fine_z[4] *= 1.090848;
  pull_fine_z[5] *= 0.967529;
  pull_fine_z[6] *= 0.820550;
  pull_fine_z[55] *= 1.548348;
  pull_fine_phi[0] *= 0.976810;
  pull_fine_phi[1] *= 0.979091;
  pull_fine_phi[2] *= 0.985423;
  pull_fine_phi[3] *= 0.848807;
  pull_fine_phi[4] *= 0.879689;
  pull_fine_phi[5] *= 0.724431;
  pull_fine_phi[6] *= 0.717564;
  pull_fine_phi[7] *= 1.087427;
  pull_fine_phi[8] *= 1.025145;
  pull_fine_phi[9] *= 1.020610;
  pull_fine_phi[10] *= 1.026286;
  pull_fine_phi[11] *= 1.019533;
  pull_fine_phi[12] *= 1.022289;
  pull_fine_phi[13] *= 1.020049;
  pull_fine_phi[14] *= 1.020093;
  pull_fine_phi[15] *= 1.017107;
  pull_fine_phi[16] *= 1.017730;
  pull_fine_phi[17] *= 1.014773;
  pull_fine_phi[18] *= 1.016468;
  pull_fine_phi[19] *= 1.014107;
  pull_fine_phi[20] *= 1.015574;
  pull_fine_phi[21] *= 1.015261;
  pull_fine_phi[22] *= 1.055006;
  pull_fine_phi[23] *= 1.081166;
  pull_fine_phi[24] *= 1.037136;
  pull_fine_phi[25] *= 1.033333;
  pull_fine_phi[26] *= 1.032861;
  pull_fine_phi[27] *= 1.032109;
  pull_fine_phi[28] *= 1.035788;
  pull_fine_phi[29] *= 1.034689;
  pull_fine_phi[30] *= 1.034438;
  pull_fine_phi[31] *= 1.036618;
  pull_fine_phi[32] *= 1.035785;
  pull_fine_phi[33] *= 1.039290;
  pull_fine_phi[34] *= 1.036625;
  pull_fine_phi[35] *= 1.038621;
  pull_fine_phi[36] *= 1.038185;
  pull_fine_phi[37] *= 1.036850;
  pull_fine_phi[38] *= 1.059049;
  pull_fine_phi[39] *= 1.122510;
  pull_fine_phi[40] *= 1.045074;
  pull_fine_phi[41] *= 1.041234;
  pull_fine_phi[42] *= 1.041916;
  pull_fine_phi[43] *= 1.038968;
  pull_fine_phi[44] *= 1.037942;
  pull_fine_phi[45] *= 1.037845;
  pull_fine_phi[46] *= 1.035315;
  pull_fine_phi[47] *= 1.036548;
  pull_fine_phi[48] *= 1.032359;
  pull_fine_phi[49] *= 1.030394;
  pull_fine_phi[50] *= 1.031460;
  pull_fine_phi[51] *= 1.032054;
  pull_fine_phi[52] *= 1.031856;
  pull_fine_phi[53] *= 1.031389;
  pull_fine_phi[54] *= 1.030922;
  pull_fine_phi[55] *= 1.010915;
  pull_fine_phi[56] *= 1.271132;
  pull_fine_phi[57] *= 0.000000;
  pull_fine_z[0] *= 0.903223;
  pull_fine_z[1] *= 0.906998;
  pull_fine_z[2] *= 0.909652;
  pull_fine_z[3] *= 1.379617;
  pull_fine_z[4] *= 0.877200;
  pull_fine_z[5] *= 0.944639;
  pull_fine_z[6] *= 1.130620;
  pull_fine_z[7] *= 1.074286;
  pull_fine_z[8] *= 1.031207;
  pull_fine_z[9] *= 1.026923;
  pull_fine_z[10] *= 1.029683;
  pull_fine_z[11] *= 1.023258;
  pull_fine_z[12] *= 1.024277;
  pull_fine_z[13] *= 1.022564;
  pull_fine_z[14] *= 1.023237;
  pull_fine_z[15] *= 1.020460;
  pull_fine_z[16] *= 1.022844;
  pull_fine_z[17] *= 1.022628;
  pull_fine_z[18] *= 1.020984;
  pull_fine_z[19] *= 1.018908;
  pull_fine_z[20] *= 1.019514;
  pull_fine_z[21] *= 1.020316;
  pull_fine_z[22] *= 1.044147;
  pull_fine_z[23] *= 1.054320;
  pull_fine_z[24] *= 1.031116;
  pull_fine_z[25] *= 1.029622;
  pull_fine_z[26] *= 1.026101;
  pull_fine_z[27] *= 1.025649;
  pull_fine_z[28] *= 1.027153;
  pull_fine_z[29] *= 1.025389;
  pull_fine_z[30] *= 1.024853;
  pull_fine_z[31] *= 1.026242;
  pull_fine_z[32] *= 1.027198;
  pull_fine_z[33] *= 1.027609;
  pull_fine_z[34] *= 1.024212;
  pull_fine_z[35] *= 1.026990;
  pull_fine_z[36] *= 1.026987;
  pull_fine_z[37] *= 1.025806;
  pull_fine_z[38] *= 1.036276;
  pull_fine_z[39] *= 1.029688;
  pull_fine_z[40] *= 1.028911;
  pull_fine_z[41] *= 1.026570;
  pull_fine_z[42] *= 1.023386;
  pull_fine_z[43] *= 1.025188;
  pull_fine_z[44] *= 1.024090;
  pull_fine_z[45] *= 1.022645;
  pull_fine_z[46] *= 1.023954;
  pull_fine_z[47] *= 1.022967;
  pull_fine_z[48] *= 1.021187;
  pull_fine_z[49] *= 1.022650;
  pull_fine_z[50] *= 1.021918;
  pull_fine_z[51] *= 1.021044;
  pull_fine_z[52] *= 1.021689;
  pull_fine_z[53] *= 1.020743;
  pull_fine_z[54] *= 1.022534;
  pull_fine_z[55] *= 0.784192;
  pull_fine_z[56] *= 0.767927;
  pull_fine_z[57] *= 0.000000;
  pull_fine_z[55] *= 0.876383;
  pull_fine_z[56] *= 0.766922;
  pull_fine_z[3] *= 2.000261;
  pull_fine_z[4] *= 1.127752;
  pull_fine_z[5] *= 0.804010;
  pull_fine_z[6] *= 0.567351;
}

//_________________________________________________________________________________
ClusterErrorPara::error_t ClusterErrorPara::get_clusterv5_modified_error(TrkrCluster* cluster, double /*unused*/, TrkrDefs::cluskey key)
{
  int layer = TrkrDefs::getLayer(key);

  double phierror = cluster->getRPhiError();
  double zerror = cluster->getZError();
  if (TrkrDefs::getTrkrId(key) == TrkrDefs::tpcId)
  {
    if (layer == 7 || layer == 22 || layer == 23 || layer == 38 || layer == 39)
    {
      phierror *= 4;
      zerror *= 4;
    }
    if (cluster->getEdge() >= 3)
    {
      phierror *= 4;
    }
    if (cluster->getOverlap() >= 2)
    {
      phierror *= 2;
    }
    if (cluster->getPhiSize() == 1)
    {
      phierror *= 10;
    }
    if (cluster->getPhiSize() >= 5)
    {
      phierror *= 10;
    }

    if (phierror > 0.1)
    {
      phierror = 0.1;
    }
    if (phierror < 0.0005)
    {
      phierror = 0.1;
    }
  }
  return std::make_pair(square(phierror), square(zerror));
}

//_________________________________________________________________________________
ClusterErrorPara::error_t ClusterErrorPara::get_cluster_error(TrkrCluster* cluster, double cluster_r, TrkrDefs::cluskey key, float qOverR, float slope)
{
  float r = cluster_r;
  float R = TMath::Abs(1.0 / qOverR);
  double alpha = (r * r) / (2 * r * R);
  double beta = TMath::Abs(atan(slope));
  return get_cluster_error(cluster, key, alpha, beta);
}

double ClusterErrorPara::tpc_phi_error(int layer, double alpha, TrkrCluster* cluster)
{
  double phierror = 0;

  int sector = -1;
  if (layer >= 7 && layer < 23)
  {
    sector = 0;
  }
  else if (layer >= 23 && layer < 39)
  {
    sector = 1;
  }
  else if (layer >= 39 && layer < 55)
  {
    sector = 2;
  }
  phierror = 0.0005;

  if (sector != -1)
  {
    phierror = 0.0005;
  }

  if (sector == 0)
  {
    // phierror = 0.019886;
    phierror = f0->Eval(alpha);
    if (cluster->getMaxAdc() != 0)
    {
      if (cluster->getMaxAdc() > 150)
      {
        phierror *= 0.54 * 0.9;
      }
      else
      {
        phierror *= fadcphi0->Eval(cluster->getMaxAdc());
        phierror *= fadcphi0fine->Eval(cluster->getMaxAdc());
      }
    }
    if (cluster->getEdge() >= 5)
    {
      phierror *= 2;
    }

    if (cluster->getPhiSize() == 1)
    {
      phierror *= 1.5;
    }
    if (cluster->getPhiSize() == 4)
    {
      phierror *= 1.5;
    }
    if (cluster->getPhiSize() >= 5)
    {
      phierror *= 2.5;
    }

    phierror *= f0fine->Eval(alpha);
  }

  if (sector == 1)
  {
    // phierror = 0.018604;
    phierror = f1->Eval(alpha);
    if (cluster->getMaxAdc() != 0)
    {
      if (cluster->getMaxAdc() > 160)
      {
        phierror *= 0.6;
      }
      else
      {
        phierror *= fadcphi1->Eval(cluster->getMaxAdc());
      }
      if (cluster->getEdge() >= 5)
      {
        phierror *= 2;
      }
      if (cluster->getMaxAdc() > 140)
      {
        phierror *= 0.95;
      }
      else
      {
        phierror *= fadcphi1fine->Eval(cluster->getMaxAdc());
      }
    }
    phierror *= 0.975;
    if (cluster->getPhiSize() == 1)
    {
      phierror *= 4;
    }
    if (cluster->getPhiSize() == 4)
    {
      phierror *= 1.2;
    }
    if (cluster->getPhiSize() >= 5)
    {
      phierror *= 2;
    }

    phierror *= f1fine->Eval(alpha);
  }

  if (sector == 2)
  {
    // phierror = 0.02043;

    phierror = f2->Eval(alpha);
    if (cluster->getMaxAdc())
    {
      if (cluster->getMaxAdc() > 170)
      {
        phierror *= 0.6 * 0.95;
      }
      else
      {
        phierror *= fadcphi2->Eval(cluster->getMaxAdc());
        if (cluster->getMaxAdc() < 100)
        {
          phierror *= fadcphi2fine1->Eval(cluster->getMaxAdc());
        }
      }
    }
    if (cluster->getEdge() >= 5)
    {
      phierror *= 2;
    }

    if (cluster->getPhiSize() == 1)
    {
      phierror *= 10;
    }
    if (cluster->getPhiSize() >= 6)
    {
      phierror *= 10;
    }

    phierror *= f2fine->Eval(alpha);
  }
  if (layer == 7)
  {
    phierror *= (3 * 0.72);
  }
  if (layer == 22 || layer == 23)
  {
    phierror *= (3 * 0.8);
  }
  if (layer == 38)
  {
    phierror *= (3 * 0.9);
  }
  if (layer == 39)
  {
    phierror *= (3 * 1.05);
  }

  if (phierror > 0.1)
  {
    phierror = 0.1;
  }
  if (phierror < 0.0005)
  {
    phierror = 0.1;
  }

  return phierror;
}
double ClusterErrorPara::tpc_z_error(int layer, double beta, TrkrCluster* cluster)
{
  double zerror = 0.05;

  int sector = -1;
  if (layer >= 7 && layer < 23)
  {
    sector = 0;
  }
  else if (layer >= 23 && layer < 39)
  {
    sector = 1;
  }
  else if (layer >= 39 && layer < 55)
  {
    sector = 2;
  }
  if (sector != -1)
  {
    zerror = 0.05;
  }

  if (sector == 0)
  {
    zerror = fz0->Eval(beta);
    if (cluster->getMaxAdc() > 180)
    {
      zerror *= 0.5;
    }
    else
    {
      zerror *= fadcz0->Eval(cluster->getMaxAdc());
    }
    zerror *= fz0fine->Eval(beta);
    zerror *= fadcz0fine->Eval(cluster->getMaxAdc());
  }

  if (sector == 1)
  {
    zerror = fz1->Eval(beta);
    if (cluster->getMaxAdc() > 180)
    {
      zerror *= 0.6;
    }
    else
    {
      zerror *= fadcz1->Eval(cluster->getMaxAdc());
    }
    zerror *= fz1fine->Eval(beta);
    zerror *= fadcz1fine->Eval(cluster->getMaxAdc());
    zerror *= 0.98;
    //    zerror *= 1.05913
  }
  if (sector == 2)
  {
    zerror = fz2->Eval(beta);
    if (cluster->getMaxAdc() > 170)
    {
      zerror *= 0.6;
    }
    else
    {
      zerror *= fadcz2->Eval(cluster->getMaxAdc());
    }
    zerror *= fz2fine->Eval(beta);
    zerror *= fadcz2fine->Eval(cluster->getMaxAdc());
    // zerrror *= 1.15575;
  }
  if (layer == 7)
  {
    zerror *= (3.5 * 1.13);
  }
  if (layer == 22)
  {
    zerror *= (3.5 * 0.9);
  }
  if (layer == 23)
  {
    zerror *= (3.5 * 0.95);
  }
  if (layer == 38)
  {
    zerror *= (3.5 * 0.85);
  }

  if (zerror > 0.2)
  {
    zerror = 0.2;
  }
  if (zerror < 0.0002)
  {
    zerror = 0.2;
  }

  return zerror;
}

double ClusterErrorPara::mm_phi_error(int layer, double alpha, TrkrCluster* cluster)
{
  double phierror = 0;
  if (layer == 55)
  {
    if (cluster->getPhiSize() == 1)
    {
      phierror = pitcherr_phi_mm1;
    }
    else if (cluster->getPhiSize() == 2)
    {
      phierror = fmm_55_2->Eval(alpha);
    }
    else if (cluster->getPhiSize() >= 3)
    {
      phierror = fmm_3->Eval(alpha);
    }
    phierror *= scale_mm_0;
  }
  else if (layer == 56)
  {
    phierror = pitcherr_phi_mm2;
  }
  return phierror;
}

double ClusterErrorPara::mm_z_error(int layer, double beta, TrkrCluster* cluster)
{
  double zerror = 0;
  if (layer == 55)
  {
    zerror = pitcherr_z_mm1;
  }
  else if (layer == 56)
  {
    if (cluster->getZSize() == 1)
    {
      zerror = pitcherr_z_mm2;
    }
    else if (cluster->getZSize() == 2)
    {
      zerror = fmm_56_2->Eval(beta);
    }
    else if (cluster->getZSize() >= 3)
    {
      zerror = fmm_3->Eval(beta);
    }
    zerror *= scale_mm_1;
  }
  return zerror;
}

double ClusterErrorPara::mvtx_phi_error(TrkrCluster* cluster)
{
  double phierror = 0;
  static constexpr std::array<double, 7> scalefactors_mvtx_phi = {{0.36, 0.6, 0.37, 0.49, 0.4, 0.37, 0.33}};
  phierror = pitcherr_phi_mvtx;  // 0.00077595876
  //  std::cout << " phistart: " << phierror  << " phisize: " << cluster->getPhiSize() << "zsize: " << cluster->getZSize() << std::endl;
  if (cluster->getPhiSize() == 1 && cluster->getZSize() == 1)
  {
    phierror = pitcherr_phi_mvtx * scalefactors_mvtx_phi[0];
  }
  else if (cluster->getPhiSize() == 2 && cluster->getZSize() == 1)
  {
    phierror = pitcherr_phi_mvtx * scalefactors_mvtx_phi[1];
  }
  else if (cluster->getPhiSize() == 1 && cluster->getZSize() == 2)
  {
    phierror = pitcherr_phi_mvtx * scalefactors_mvtx_phi[2];
  }
  else if (cluster->getPhiSize() == 2 && cluster->getZSize() == 2)
  {
    phierror = pitcherr_phi_mvtx * scalefactors_mvtx_phi[0];
  }
  else if (cluster->getPhiSize() == 2 && cluster->getZSize() == 3)
  {
    phierror = pitcherr_phi_mvtx * scalefactors_mvtx_phi[1];
  }
  else if (cluster->getPhiSize() == 3 && cluster->getZSize() == 2)
  {
    phierror = pitcherr_phi_mvtx * scalefactors_mvtx_phi[2];
  }
  else if (cluster->getPhiSize() == 3 && cluster->getZSize() == 3)
  {
    phierror = pitcherr_phi_mvtx * scalefactors_mvtx_phi[3];
  }
  else
  {
    phierror = pitcherr_phi_mvtx * cluster->getPhiSize();
  }
  //  std::cout << " phi after size: " << phierror  << std::endl;
  phierror *= scale_mvtx;
  // else zerror = pitcherr_z_mvtx*cluster->getZSize();
  //  std::cout << " phi: " << phierror  << std::endl;
  return phierror;
}

double ClusterErrorPara::mvtx_phi_error(const TrkrCluster* cluster)
{
  double phierror = 0;
  static constexpr std::array<double, 7> scalefactors_mvtx_phi = {{0.36, 0.6, 0.37, 0.49, 0.4, 0.37, 0.33}};
  phierror = pitcherr_phi_mvtx;  // 0.00077595876
  //  std::cout << " phistart: " << phierror  << " phisize: " << cluster->getPhiSize() << "zsize: " << cluster->getZSize() << std::endl;
  if (cluster->getPhiSize() == 1 && cluster->getZSize() == 1)
  {
    phierror = pitcherr_phi_mvtx * scalefactors_mvtx_phi[0];
  }
  else if (cluster->getPhiSize() == 2 && cluster->getZSize() == 1)
  {
    phierror = pitcherr_phi_mvtx * scalefactors_mvtx_phi[1];
  }
  else if (cluster->getPhiSize() == 1 && cluster->getZSize() == 2)
  {
    phierror = pitcherr_phi_mvtx * scalefactors_mvtx_phi[2];
  }
  else if (cluster->getPhiSize() == 2 && cluster->getZSize() == 2)
  {
    phierror = pitcherr_phi_mvtx * scalefactors_mvtx_phi[0];
  }
  else if (cluster->getPhiSize() == 2 && cluster->getZSize() == 3)
  {
    phierror = pitcherr_phi_mvtx * scalefactors_mvtx_phi[1];
  }
  else if (cluster->getPhiSize() == 3 && cluster->getZSize() == 2)
  {
    phierror = pitcherr_phi_mvtx * scalefactors_mvtx_phi[2];
  }
  else if (cluster->getPhiSize() == 3 && cluster->getZSize() == 3)
  {
    phierror = pitcherr_phi_mvtx * scalefactors_mvtx_phi[3];
  }
  else
  {
    phierror = pitcherr_phi_mvtx * cluster->getPhiSize();
  }
  //  std::cout << " phi after size: " << phierror  << std::endl;
  phierror *= scale_mvtx;
  // else zerror = pitcherr_z_mvtx*cluster->getZSize();
  //  std::cout << " phi: " << phierror  << std::endl;
  return phierror;
}

double ClusterErrorPara::mvtx_z_error(TrkrCluster* cluster)
{
  double zerror = 0;

  zerror = pitcherr_z_mvtx;
  static constexpr std::array<double, 4> scalefactors_z = {{0.47, 0.48, 0.71, 0.55}};
  if (cluster->getZSize() == 2 && cluster->getPhiSize() == 2)
  {
    zerror = pitcherr_z_mvtx * scalefactors_z[0];
  }
  else if (cluster->getZSize() == 2 && cluster->getPhiSize() == 3)
  {
    zerror = pitcherr_z_mvtx * scalefactors_z[1];
  }
  else if (cluster->getZSize() == 3 && cluster->getPhiSize() == 2)
  {
    zerror = pitcherr_z_mvtx * scalefactors_z[2];
  }
  else if (cluster->getZSize() == 3 && cluster->getPhiSize() == 3)
  {
    zerror = pitcherr_z_mvtx * scalefactors_z[3];
  }
  zerror *= scale_mvtx_z;
  return zerror;
}
double ClusterErrorPara::mvtx_z_error(const TrkrCluster* cluster)
{
  double zerror = 0;

  zerror = pitcherr_z_mvtx;
  static constexpr std::array<double, 4> scalefactors_z = {{0.47, 0.48, 0.71, 0.55}};
  if (cluster->getZSize() == 2 && cluster->getPhiSize() == 2)
  {
    zerror = pitcherr_z_mvtx * scalefactors_z[0];
  }
  else if (cluster->getZSize() == 2 && cluster->getPhiSize() == 3)
  {
    zerror = pitcherr_z_mvtx * scalefactors_z[1];
  }
  else if (cluster->getZSize() == 3 && cluster->getPhiSize() == 2)
  {
    zerror = pitcherr_z_mvtx * scalefactors_z[2];
  }
  else if (cluster->getZSize() == 3 && cluster->getPhiSize() == 3)
  {
    zerror = pitcherr_z_mvtx * scalefactors_z[3];
  }
  zerror *= scale_mvtx_z;
  return zerror;
}

double ClusterErrorPara::intt_phi_error(int layer, TrkrCluster* cluster)
{
  double phierror = 0;
  static constexpr std::array<double, 3> scalefactors_intt_phi = {{0.85, 0.4, 0.33}};

  if (cluster->getPhiSize() == 1 && layer < 5)
  {
    phierror = pitcherr_phi_intt * scalefactors_intt_phi[0];
  }
  else if (cluster->getPhiSize() == 2 && layer < 5)
  {
    phierror = pitcherr_phi_intt * scalefactors_intt_phi[1];
  }
  else if (cluster->getPhiSize() == 2 && layer > 4)
  {
    phierror = pitcherr_phi_intt * scalefactors_intt_phi[2];
  }
  else
  {
    phierror = pitcherr_phi_intt * cluster->getPhiSize();
  }
  if (layer == 3)
  {
    phierror *= scale_intt_3;
  }
  if (layer == 4)
  {
    phierror *= scale_intt_4;
  }
  if (layer == 5)
  {
    phierror *= scale_intt_5;
  }
  if (layer == 6)
  {
    phierror *= scale_intt_6;
  }

  return phierror;
}

double ClusterErrorPara::intt_z_error(TrkrCluster* cluster)
{
  double zerror = 0;
  zerror = pitcherr_z_intt * cluster->getZSize();
  return zerror;
}
double ClusterErrorPara::intt_phi_error(int layer, const TrkrCluster* cluster)
{
  double phierror = 0;
  static constexpr std::array<double, 3> scalefactors_intt_phi = {{0.85, 0.4, 0.33}};

  if (cluster->getPhiSize() == 1 && layer < 5)
  {
    phierror = pitcherr_phi_intt * scalefactors_intt_phi[0];
  }
  else if (cluster->getPhiSize() == 2 && layer < 5)
  {
    phierror = pitcherr_phi_intt * scalefactors_intt_phi[1];
  }
  else if (cluster->getPhiSize() == 2 && layer > 4)
  {
    phierror = pitcherr_phi_intt * scalefactors_intt_phi[2];
  }
  else
  {
    phierror = pitcherr_phi_intt * cluster->getPhiSize();
  }
  if (layer == 3)
  {
    phierror *= scale_intt_3;
  }
  if (layer == 4)
  {
    phierror *= scale_intt_4;
  }
  if (layer == 5)
  {
    phierror *= scale_intt_5;
  }
  if (layer == 6)
  {
    phierror *= scale_intt_6;
  }

  return phierror;
}

double ClusterErrorPara::intt_z_error(const TrkrCluster* cluster)
{
  double zerror = 0;
  zerror = pitcherr_z_intt * cluster->getZSize();
  return zerror;
}

//_________________________________________________________________________________
ClusterErrorPara::error_t ClusterErrorPara::get_cluster_error(TrkrCluster* cluster, TrkrDefs::cluskey key, double alpha, double beta)
{
  int layer = TrkrDefs::getLayer(key);
  double phierror = 0;
  double zerror = 0;
  switch (TrkrDefs::getTrkrId(key))
  {
  default:
    break;

  case TrkrDefs::micromegasId:
    phierror = mm_phi_error(layer, alpha, cluster);
    zerror = mm_z_error(layer, beta, cluster);

    break;
  case TrkrDefs::mvtxId:

    phierror = mvtx_phi_error(cluster);
    zerror = mvtx_z_error(cluster);
    //      std::cout << "   z: " << zerror  << " l: " << layer  << std::endl;
    break;

  case TrkrDefs::inttId:

    phierror = intt_phi_error(layer, cluster);
    zerror = intt_z_error(cluster);
    break;

  case TrkrDefs::tpcId:
    phierror = tpc_phi_error(layer, alpha, cluster);
    zerror = tpc_z_error(layer, beta, cluster);
    /*
    std::cout << " phi: " << phierror << " | " << phierror2 << " l: " << layer << " a " << alpha << " maxadc: " << cluster->getMaxAdc() << std::endl;
    std::cout << "   z: " << zerror  << " | " << zerror2 << " l: " << layer << " b " << beta << " maxadc: " << cluster->getMaxAdc() << " adc: " << cluster->getAdc() << std::endl;
    */
    break;
  }
  if (pull_fine_phi[layer] != 0)
  {
    phierror *= pull_fine_phi[layer];
  }
  if (pull_fine_z[layer] != 0)
  {
    zerror *= pull_fine_z[layer];
  }

  if (phierror == 0)
  {
    phierror = 100;
  }
  if (zerror == 0)
  {
    zerror = 100;
  }

  return std::make_pair(square(phierror), square(zerror));
  //  return std::make_pair(phierror,zerror);
}
//_________________________________________________________________________________

ClusterErrorPara::error_t ClusterErrorPara::get_simple_cluster_error(TrkrCluster* cluster, double cluster_r, TrkrDefs::cluskey key)
{
  double alpha = 0.17;
  double beta = 0.4;
  if (cluster_r > 100000)
  {
    alpha = 0.17001;
  }
  int layer = TrkrDefs::getLayer(key);
  double phierror = 0;
  double zerror = 0;
  switch (TrkrDefs::getTrkrId(key))
  {
  default:
    break;

  case TrkrDefs::micromegasId:
    phierror = mm_phi_error(layer, 0.27, cluster);
    zerror = mm_z_error(layer, 0.4, cluster);

    break;
  case TrkrDefs::mvtxId:

    phierror = mvtx_phi_error(cluster);
    zerror = mvtx_z_error(cluster);
    //      std::cout << "   z: " << zerror  << " l: " << layer  << std::endl;
    break;

  case TrkrDefs::inttId:

    phierror = intt_phi_error(layer, cluster);
    zerror = intt_z_error(cluster);
    break;

  case TrkrDefs::tpcId:
    phierror = tpc_phi_error(layer, alpha, cluster);
    zerror = tpc_z_error(layer, beta, cluster);

    break;
  }
  if (pull_fine_phi[layer] != 0)
  {
    phierror *= pull_fine_phi[layer];
  }
  if (pull_fine_z[layer] != 0)
  {
    zerror *= pull_fine_z[layer];
  }

  if (phierror == 0)
  {
    phierror = 100;
  }
  if (zerror == 0)
  {
    zerror = 100;
  }

  return std::make_pair(square(phierror), square(zerror));
}

//_________________________________________________________________________________
ClusterErrorPara::error_t ClusterErrorPara::get_fix_tpc_cluster_error(TrkrCluster* cluster, TrkrDefs::cluskey key)
{
  int layer = TrkrDefs::getLayer(key);
  double phierror = 0;
  double zerror = 0;

  phierror = mvtx_phi_error(cluster);
  zerror = mvtx_z_error(cluster);

  if (pull_fine_phi[layer] != 0)
  {
    phierror *= pull_fine_phi[layer];
  }
  if (pull_fine_z[layer] != 0)
  {
    zerror *= pull_fine_z[layer];
  }

  if (phierror == 0)
  {
    phierror = 100;
  }
  if (zerror == 0)
  {
    zerror = 100;
  }

  return std::make_pair(square(phierror), square(zerror));
}

//_________________________________________________________________________________
ClusterErrorPara::error_t ClusterErrorPara::get_si_cluster_error(const TrkrCluster* cluster, TrkrDefs::cluskey key)
{
  int layer = TrkrDefs::getLayer(key);
  double phierror = 0;
  double zerror = 0;
  switch (TrkrDefs::getTrkrId(key))
  {
  default:
    break;

  case TrkrDefs::mvtxId:

    phierror = mvtx_phi_error(cluster);
    zerror = mvtx_z_error(cluster);
    //      std::cout << "   z: " << zerror  << " l: " << layer  << std::endl;
    break;

  case TrkrDefs::inttId:

    phierror = intt_phi_error(layer, cluster);
    zerror = intt_z_error(cluster);
    break;
  }
  if (pull_fine_phi[layer] != 0)
  {
    phierror *= pull_fine_phi[layer];
  }
  if (pull_fine_z[layer] != 0)
  {
    zerror *= pull_fine_z[layer];
  }

  if (phierror == 0)
  {
    phierror = 100;
  }
  if (zerror == 0)
  {
    zerror = 100;
  }

  return std::make_pair(square(phierror), square(zerror));
}
