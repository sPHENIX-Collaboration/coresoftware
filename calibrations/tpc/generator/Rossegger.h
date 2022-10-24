#ifndef __ROSSEGGER_H__
#define __ROSSEGGER_H__

//
//  Hello Space Charge Fans:  (although you should hate space charge)
//
//    This is a code that implements the calculations of space charge contributions
//  In a cylindrical TPC.  It uses the solutions discovered by Stefan Rossegger for
//  the ALICE experiment.  These solutions use various Bessel functions as a series
//  solution to the "point charge in a conducting cylinder" problem imposed by
//  The typical configuration of a TPC found at a collider.
//
//    These calculations have only a single dimension [length] and we shall choose
//  the units of cm throughout the calculations.
//
//                                                TKH
//                                                12-3-2015
//

#include <cmath>
#include <cstdio>
#include <map>
#include <string>

class TH2;
class TH3;

#define NumberOfOrders 15  // Convergence problems after 15; Rossegger used 30

class Rossegger
{
 public:
  explicit Rossegger(std::string filename);
  Rossegger(double a = 30, double b = 80, double L = 80, double epsilon = 1E-4);
  virtual ~Rossegger() {}

  void Verbosity(int v)
  {
    printf("verbosity set to %d.  was %d\n", v, verbosity);
    verbosity = v;
    return;
  };
  double Rmn(int m, int n, double r);               //Rmn function from Rossegger
  double Rmn_for_zeroes(int m, double x);           //Rmn function from Rossegger, as used to find Betamn zeroes.
  double Rmn1(int m, int n, double r);              //Rmn1 function from Rossegger
  double Rmn2(int m, int n, double r);              //Rmn2 function from Rossegger
  double RPrime(int m, int n, double a, double r);  // RPrime function from Rossegger

  double Rnk(int n, int k, double r);       //Rnk function from Rossegger
  double Rnk_for_zeroes(int n, double mu);  //Rnk function from Rossegger, as used to find munk zeroes.

  double Limu(double mu, double x);  //Bessel functions of purely imaginary order
  double Kimu(double mu, double x);  //Bessel functions of purely imaginary order

  double Ez(double r, double phi, double z, double r1, double phi1, double z1);
  double Er(double r, double phi, double z, double r1, double phi1, double z1);
  double Ephi(double r, double phi, double z, double r1, double phi1, double z1);

  //alternate versions that don't use precalc constants.
  double Rmn_(int m, int n, double r);  //Rmn function from Rossegger
  //Rmn_for_zeroes doesn't have a way to speed it up with precalcs.
  double Rmn1_(int m, int n, double r);              //Rmn1 function from Rossegger
  double Rmn2_(int m, int n, double r);              //Rmn2 function from Rossegger
  double RPrime_(int m, int n, double a, double r);  // RPrime function from Rossegger
  double Rnk_(int n, int k, double r);               //Rnk function from Rossegger
  double Rnk_for_zeroes_(int n, double mu);          //Rnk function from Rossegger, as used to find munk zeroes.
  double Ez_(double r, double phi, double z, double r1, double phi1, double z1);
  double Er_(double r, double phi, double z, double r1, double phi1, double z1);
  double Ephi_(double r, double phi, double z, double r1, double phi1, double z1);

 protected:
  bool fByFile = false;
  double a = NAN;
  double b = NAN;
  double L = NAN;  //  InnerRadius, OuterRadius, Length of 1/2 the TPC.
  int verbosity = 0;
  double pi = M_PI;
  double epsilon = NAN;  //precision.

  bool tweak = false;

  double MinimumDR = NAN;
  double MinimumDPHI = NAN;
  double MinimumDZ = NAN;

  double FindNextZero(double xstart, double epsilon, int order, double (Rossegger::*func)(int, double));  // Routine to find zeroes of func.
  void FindBetamn(double epsilon);                                                                        // Routine used to fill the Betamn array with resolution epsilon...
  void FindMunk(double epsilon);                                                                          // Routine used to fill the Munk array with resolution epsilon...
  bool CheckZeroes(double epsilon);                                                                       //confirm that the zeroes match to the desired precision.

  void LoadZeroes(const char* destfile);
  void SaveZeroes(const char* destfile);

  double Betamn[NumberOfOrders][NumberOfOrders];  //  Betamn array from Rossegger
  double N2mn[NumberOfOrders][NumberOfOrders];    //  N2mn array from Rossegger
  double Munk[NumberOfOrders][NumberOfOrders];    //  Munk array from Rossegger
  double N2nk[NumberOfOrders][NumberOfOrders];    //  N2nk array from Rossegger

  void PrecalcFreeConstants();  //Routine used to fill the arrays of other values used repeatedly:
  //helpful values to precompute for frequent use:
  double BetaN[NumberOfOrders];                               //  BetaN=(n+1)*pi/L as used in eg 5.32, .46
  double BetaN_a[NumberOfOrders];                             //  BetaN*a as used in eg 5.32, .46
  double BetaN_b[NumberOfOrders];                             //  BetaN*b as used in eg 5.32, .46
  double km_BetaN_a[NumberOfOrders][NumberOfOrders];          //kn(m,BetaN*a) as used in Rossegger 5.32
  double im_BetaN_a[NumberOfOrders][NumberOfOrders];          //in(m,BetaN*a) as used in Rossegger 5.32
  double km_BetaN_b[NumberOfOrders][NumberOfOrders];          //kn(m,BetaN*b) as used in Rossegger 5.33
  double im_BetaN_b[NumberOfOrders][NumberOfOrders];          //in(m,BetaN*b) as used in Rossegger 5.33
  double bessel_denominator[NumberOfOrders][NumberOfOrders];  //BesselI(m,BetaN*a)*BesselK(m,BetaN*b)-BesselI(m,BetaN*b)*BesselK(m,BetaN*a) as in Rossegger 5.65

  void PrecalcDerivedConstants();                         //Routine used to fill repeated values that depend on the Betamn and Munk zeroes:
  double ym_Betamn_a[NumberOfOrders][NumberOfOrders];     //yn(m,Betamn[m][n]*a) as used in Rossegger 5.11
  double jm_Betamn_a[NumberOfOrders][NumberOfOrders];     //jn(m,Betamn[m][n]*a) as used in Rossegger 5.11
  double liMunk_BetaN_a[NumberOfOrders][NumberOfOrders];  //limu(Munk[n][k],BetaN*a) as used in Rossegger 5.45
  double kiMunk_BetaN_a[NumberOfOrders][NumberOfOrders];  //kimu(Munk[n][k],BetaN*a) as used in Rossegger 5.45
  double sinh_Betamn_L[NumberOfOrders][NumberOfOrders];   //sinh(Betamn[m][n]*L)  as in Rossegger 5.64
  double sinh_pi_Munk[NumberOfOrders][NumberOfOrders];    //sinh(pi*Munk[n][k]) as in Rossegger 5.66

  TH2* Tags = nullptr;
  std::map<std::string, TH3*> Grid;
};

#endif /* __SPACECHARGE_H__ */
