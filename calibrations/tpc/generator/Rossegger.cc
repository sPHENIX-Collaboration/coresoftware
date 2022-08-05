#include "Rossegger.h"

#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include <boost/math/special_functions.hpp>  //covers all the special functions.

#include <algorithm>  // for max
#include <cassert>    // for assert
#include <cmath>
#include <cstdlib>  // for exit, abs
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// the limu and kimu terms, that i need to think about a little while longer...
extern "C"
{
  void dkia_(int *IFAC, double *X, double *A, double *DKI, double *DKID, int *IERRO);
  void dlia_(int *IFAC, double *X, double *A, double *DLI, double *DLID, int *IERRO);
}
//

//Bessel Function J_n(x):
#define jn(order, x) boost::math::cyl_bessel_j(order, x)
//Bessel (Neumann) Function Y_n(x):
#define yn(order, x) boost::math::cyl_neumann(order, x)
//Modified Bessel Function of the first kind I_n(x):
#define in(order, x) boost::math::cyl_bessel_i(order, x)
//Modified Bessel Function of the second kind K_n(x):
#define kn(order, x) boost::math::cyl_bessel_k(order, x)
#define limu(im_order, x) Rossegger::Limu(im_order, x)
#define kimu(im_order, x) Rossegger::Kimu(im_order, x)

/*
  This is a modified/renamed copy of Carlos and Tom's "Spacecharge" class, modified to use boost instead of fortran routines, and with phi terms added.
 */

Rossegger::Rossegger(double InnerRadius, double OuterRadius, double Rdo_Z, double precision)
{
  a = InnerRadius;
  b = OuterRadius;
  L = Rdo_Z;

  epsilon = precision;

  verbosity = 0;
  pi = M_PI;

  PrecalcFreeConstants();

  //load the greens functions:
  char zeroesfilename[200];
  sprintf(zeroesfilename, "rosseger_zeroes_eps%1.0E_a%2.2f_b%2.2f_L%2.2f.root", epsilon, a, b, L);
  TFile *fileptr = TFile::Open(zeroesfilename, "READ");
  if (!fileptr)
  {  //generate the lookuptable
    FindMunk(epsilon);
    FindBetamn(epsilon);
    SaveZeroes(zeroesfilename);
  }
  else
  {  //load it from a file
    fileptr->Close();
    LoadZeroes(zeroesfilename);
  }

  PrecalcDerivedConstants();

  std::cout << "Rossegger object initialized as follows:" << std::endl;
  std::cout << "  Inner Radius = " << a << " cm." << std::endl;
  std::cout << "  Outer Radius = " << b << " cm." << std::endl;
  std::cout << "  Half  Length = " << L << " cm." << std::endl;

  if (!CheckZeroes(0.01))
  {
    std::cout << "CheckZeroes(0.01) returned false, exiting" << std::endl;
    exit(1);
  }
  return;
}

double Rossegger::FindNextZero(double xstart, double epsilon, int order, double (Rossegger::*func)(int, double))
{
  double previous = (this->*func)(order, xstart);
  double x = xstart + epsilon;
  double value = previous;

  while (!((value == 0) || (value < 0 && previous > 0) || (value > 0 && previous < 0)))
  {
    //  Rossegger equation 5.12
    value = (this->*func)(order, x);
    if (value == 0) std::cout << "hit it exactly!  Go buy a lottery ticket!" << std::endl;
    if ((value == 0) || (value < 0 && previous > 0) || (value > 0 && previous < 0))
    {
      //when we go from one sign to the other, we have bracketed the zero
      //the following is mathematically equivalent to finding the delta
      //between the x of our previous value and the point where we cross the axis
      //then returning x0=x_old+delta
      double slope = (value - previous) / epsilon;
      double intercept = value - slope * x;
      double x0 = -intercept / slope;
      if (verbosity > 1) std::cout << " " << x0 << "," << std::endl;
      double n0 = (this->*func)(order, x - epsilon);
      double n1 = (this->*func)(order, x + epsilon);
      if ((n0 < 0 && n1 < 0) || (n0 > 0 && n1 > 0))
      {
        printf("neighbors on both sides have the same sign.  Check your function resolution!\n");
      }

      return x0;
    }
    previous = value;
    x += epsilon;
  }
  std::cout << "logic break!\n";
  assert(1 == 2);
  return 0;
}

void Rossegger::FindBetamn(double epsilon)
{
  std::cout << "Now filling the Beta[m][n] Array..." << std::endl;
  if (verbosity > 5) std::cout << "numberOfOrders= " << NumberOfOrders << std::endl;
  for (int m = 0; m < NumberOfOrders; m++)
  {
    if (verbosity) std::cout << "Filling Beta[" << m << "][n]..." << std::endl;

    double x = epsilon;
    for (int n = 0; n < NumberOfOrders; n++)
    {  //  !!!  Off by one from Rossegger convention  !!!
      x = FindNextZero(x, epsilon, m, &Rossegger::Rmn_for_zeroes);
      Betamn[m][n] = x / b;
      x += epsilon;
    }
  }

  //  Now fill in the N2mn array...
  for (int m = 0; m < NumberOfOrders; m++)
  {
    for (int n = 0; n < NumberOfOrders; n++)
    {
      //  Rossegger Equation 5.17
      //  N^2_mn = 2/(pi * beta)^2 [ {Jm(beta a)/Jm(beta b)}^2 - 1 ]
      N2mn[m][n] = 2 / (pi * pi * Betamn[m][n] * Betamn[m][n]);
      //N2mn[m][n] *= (jn(m,Betamn[m][n]*a)*jn(m,Betamn[m][n]*a)/jn(m,Betamn[m][n]*b)/jn(m,Betamn[m][n]*b) - 1.0);
      double jna_over_jnb = jn(m, Betamn[m][n] * a) / jn(m, Betamn[m][n] * b);
      N2mn[m][n] *= (jna_over_jnb * jna_over_jnb - 1.0);
      //rcc note!  in eq 5.17, N2nm is set with betamn[m][n], but from context that looks to be a typo.  The order is mn everywhere else
      if (verbosity > 1) std::cout << "m: " << m << " n: " << n << " N2[m][n]: " << N2mn[m][n];
      double step = 0.01;
      if (verbosity > 1)
      {
        double integral = 0.0;
        for (double r = a; r < b; r += step)
        {
          double rmnval = Rmn(m, n, r);

          integral += rmnval * rmnval * r * step;  //Rmn(m,n,r)*Rmn(m,n,r)*r*step;
        }
        std::cout << " Int: " << integral << std::endl;
      }
      //N2mn[m][n] = integral;
    }
  }

  std::cout << "Done." << std::endl;
}

void Rossegger::FindMunk(double epsilon)
{
  std::cout << "Now filling the Mu[n][k] Array..." << std::endl;
  if (verbosity > 5) std::cout << "numberOfOrders= " << NumberOfOrders << std::endl;
  // We're looking for the zeroes of Rossegger eqn. 5.46:
  // R_nk(mu_nk;a,b)=Limu(Beta_n*a)Kimu(Beta_n*b)-Kimu(Beta_n*a)Limu(Beta_n*b)=0
  // since a and b are fixed, R_nk is a function solely of mu_nk and n.
  // for each 'n' we wish to find the a set of k mu_n's that result in R_nk=0

  //could add an option here to load the munks from file if the dimensions match.

  for (int n = 0; n < NumberOfOrders; n++)  //  !!!  Off by one from Rossegger convention  !!!
  {
    if (verbosity) std::cout << "Filling Mu[" << n << "][k]..." << std::endl;
    double x = epsilon;
    for (int k = 0; k < NumberOfOrders; k++)
    {
      x = FindNextZero(x, epsilon, n, &Rossegger::Rnk_for_zeroes);
      Munk[n][k] = x;
      x += epsilon;
      if (verbosity > 0)
      {
        printf("Mu[%d][%d]=%E\n", n, k, Munk[n][k]);
        printf("adjacent values are Rnk[mu-epsilon]=%E\tRnk[mu+epsilon]=%E\n",
               Rnk_for_zeroes(n, x - epsilon), Rnk_for_zeroes(n, x + epsilon));
        if (verbosity > 100) printf("values of argument to limu and kimu are %f and %f\n",
                                    (n + 1) * pi / L * a, (n + 1) * pi / L * b);
      }
    }
  }

  //  Now fill in the N2nk array...
  for (int n = 0; n < NumberOfOrders; n++)
  {
    for (int k = 0; k < NumberOfOrders; k++)
    {
      //  Rossegger Equation 5.48
      //  Integral of R_nk(r)*R_ns(r) dr/r= delta_ks*N2nk
      //  note that unlike N2mn, there is no convenient shortcut here.
      double integral = 0.0;
      double step = 0.001;

      for (double r = a; r < b; r += step)
      {
        double rnkval = Rnk_(n, k, r);  //must used un-optimized.  We don't have those values yet...

        integral += rnkval * rnkval / r * step;
      }
      if (verbosity > 1)
      {
        std::cout << " Int: " << integral << std::endl;
      }
      N2nk[n][k] = integral;
      if (verbosity > 1) std::cout << "n: " << n << " k: " << k << " N2nk[n][k]: " << N2nk[n][k];
    }
  }

  std::cout << "Done." << std::endl;
  return;
}

bool Rossegger::CheckZeroes(double epsilon)
{
  //confirm that the tabulated zeroes are all zeroes of their respective functions:
  double result;
  for (int m = 0; m < NumberOfOrders; m++)
  {
    for (int n = 0; n < NumberOfOrders; n++)
    {  //  !!!  Off by one from Rossegger convention  !!!
      result = Rmn_for_zeroes(m, Betamn[m][n] * b);
      if (abs(result) > epsilon)
      {
        printf("(m=%d,n=%d) Jm(x)Ym(lx)-Jm(lx)Ym(x) = %f for x=b*%f\n", m, n, result, Betamn[m][n]);
        return false;
      }
    }
  }

  // R_nk(mu_nk;a,b)=Limu(Beta_n*a)Kimu(Beta_n*b)-Kimu(Beta_n*a)Limu(Beta_n*b)=0
  for (int n = 0; n < NumberOfOrders; n++)
  {
    for (int k = 0; k < NumberOfOrders; k++)
    {  //  !!!  Off by one from Rossegger convention  !!!
      result = Rnk_for_zeroes(n, Munk[n][k]);
      if (abs(result) > epsilon * 100)
      {
        printf("(n=%d,k=%d) limu(npi*a/L)kimu(npi*b/L)-kimu(npi*a/L)kimu(npi*b/L) = %f (>eps*100) for mu=%f\n",
               n, k, result, Munk[n][k]);
        return false;
      }
    }
  }

  return true;
}

void Rossegger::PrecalcFreeConstants()
{  //Routine used to fill the arrays of other values used repeatedly
  //these constants depend only on the geometry of the detector
  printf("Precalcing %d geometric constants\n", 3 * NumberOfOrders + 5 * NumberOfOrders * NumberOfOrders);
  for (int n = 0; n < NumberOfOrders; n++)
  {
    BetaN[n] = (n + 1) * pi / L;  //  BetaN=(n+1)*pi/L as used in eg 5.32, .46
    BetaN_a[n] = BetaN[n] * a;    //  BetaN*a as used in eg 5.32, .46
    BetaN_b[n] = BetaN[n] * b;    //  BetaN*b as used in eg 5.32, .46
    for (int m = 0; m < NumberOfOrders; m++)
    {
      km_BetaN_a[m][n] = kn(m, BetaN_a[n]);                                                                                                                      //kn(m,BetaN*a) as used in Rossegger 5.32
      im_BetaN_a[m][n] = in(m, BetaN_a[n]);                                                                                                                      //in(m,BetaN*a) as used in Rossegger 5.32
      km_BetaN_b[m][n] = kn(m, BetaN_b[n]);                                                                                                                      //kn(m,BetaN*b) as used in Rossegger 5.33
      im_BetaN_b[m][n] = in(m, BetaN_b[n]);                                                                                                                      //in(m,BetaN*b) as used in Rossegger 5.33
      bessel_denominator[m][n] = TMath::BesselI(m, BetaN_a[n]) * TMath::BesselK(m, BetaN_b[n]) - TMath::BesselI(m, BetaN_b[n]) * TMath::BesselK(m, BetaN_a[n]);  //TMath::BesselI(m,BetaN*a)*TMath::BesselK(m,BetaN*b)-TMath::BesselI(m,BetaN*b)*TMath::BesselK(m,BetaN*a) as in Rossegger 5.65
    }
  }
  return;
}
void Rossegger::PrecalcDerivedConstants()
{  //Routine used to fill the arrays of other values used repeatedly
   //these constants depend on geometry and the zeroes of special functions
  printf("Precalcing %d geometric constants\n", 6 * NumberOfOrders * NumberOfOrders);

  for (int n = 0; n < NumberOfOrders; n++)
  {
    for (int m = 0; m < NumberOfOrders; m++)
    {
      ym_Betamn_a[m][n] = yn(m, Betamn[m][n] * a);   //yn(m,Betamn[m][n]*a) as used in Rossegger 5.11
      jm_Betamn_a[m][n] = jn(m, Betamn[m][n] * a);   //jn(m,Betamn[m][n]*a) as used in Rossegger 5.11
      sinh_Betamn_L[m][n] = sinh(Betamn[m][n] * L);  //sinh(Betamn[m][n]*L)  as in Rossegger 5.64
    }
    for (int k = 0; k < NumberOfOrders; k++)
    {
      liMunk_BetaN_a[n][k] = limu(Munk[n][k], BetaN_a[n]);  //limu(Munk[n][k],BetaN*a) as used in Rossegger 5.45
      kiMunk_BetaN_a[n][k] = kimu(Munk[n][k], BetaN_a[n]);  //kimu(Munk[n][k],BetaN*a) as used in Rossegger 5.45
      sinh_pi_Munk[n][k] = sinh(pi * Munk[n][k]);           //sinh(pi*Munk[n][k]) as in Rossegger 5.66
    }
  }
  return;
}

double Rossegger::Limu(double mu, double x)
{
  //defined in Rossegger eqn 5.44, also a canonical 'satisfactory companion' to Kimu.
  //could use Griddit?
  //  Rossegger Equation 5.45
  //       Rnk(r) = Limu_nk (BetaN a) Kimu_nk (BetaN r) - Kimu_nk(BetaN a) Limu_nk (BetaN r)

  int IFAC = 1;
  double A = mu;
  double DLI = 0;
  double DERR = 0;
  int IERRO = 0;

  double X = x;
  dlia_(&IFAC, &X, &A, &DLI, &DERR, &IERRO);
  return DLI;
}
double Rossegger::Kimu(double mu, double x)
{
  int IFAC = 1;
  double A = mu;
  double DKI = 0;
  double DERR = 0;
  int IERRO = 0;

  double X = x;
  dkia_(&IFAC, &X, &A, &DKI, &DERR, &IERRO);
  return DKI;
}

double Rossegger::Rmn_for_zeroes(int m, double x)
{
  double lx = a * x / b;
  //  Rossegger Equation 5.12:

  return jn(m, x) * yn(m, lx) - jn(m, lx) * yn(m, x);
}

double Rossegger::Rmn(int m, int n, double r)
{
  if (verbosity > 100) std::cout << "Determine Rmn(" << m << "," << n << "," << r << ") = ";

  //  Check input arguments for sanity...
  int error = 0;
  if (m < 0 || m >= NumberOfOrders) error = 1;
  if (n < 0 || n >= NumberOfOrders) error = 1;
  if (r < a || r > b) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments Rmn(" << m << "," << n << "," << r << ")" << std::endl;
    ;
    return 0;
  }

  //  Calculate the function using C-libraries from boost
  //  Rossegger Equation 5.11:
  //         Rmn(r) = Ym(Beta_mn a)*Jm(Beta_mn r) - Jm(Beta_mn a)*Ym(Beta_mn r)
  double R = 0;
  R = ym_Betamn_a[m][n] * jn(m, Betamn[m][n] * r) - jm_Betamn_a[m][n] * yn(m, Betamn[m][n] * r);

  if (verbosity > 100) std::cout << R << std::endl;
  return R;
}

double Rossegger::Rmn_(int m, int n, double r)
{
  if (verbosity > 100) std::cout << "Determine Rmn(" << m << "," << n << "," << r << ") = ";

  //  Check input arguments for sanity...
  int error = 0;
  if (m < 0 || m >= NumberOfOrders) error = 1;
  if (n < 0 || n >= NumberOfOrders) error = 1;
  if (r < a || r > b) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments Rmn(" << m << "," << n << "," << r << ")" << std::endl;
    ;
    return 0;
  }

  //  Calculate the function using C-libraries from boost
  //  Rossegger Equation 5.11:
  //         Rmn(r) = Ym(Beta_mn a)*Jm(Beta_mn r) - Jm(Beta_mn a)*Ym(Beta_mn r)
  double R = 0;
  R = yn(m, Betamn[m][n] * a) * jn(m, Betamn[m][n] * r) - jn(m, Betamn[m][n] * a) * yn(m, Betamn[m][n] * r);

  if (verbosity > 100) std::cout << R << std::endl;
  return R;
}

double Rossegger::Rmn1(int m, int n, double r)
{
  //  Check input arguments for sanity...
  int error = 0;
  if (m < 0 || m >= NumberOfOrders) error = 1;
  if (n < 0 || n >= NumberOfOrders) error = 1;
  if (r < a || r > b) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments Rmn1(" << m << "," << n << "," << r << ")" << std::endl;
    ;
    return 0;
  }

  //  Calculate using the TMath functions from root.
  //  Rossegger Equation 5.32
  //         Rmn1(r) = Km(BetaN a)Im(BetaN r) - Im(BetaN a) Km(BetaN r)
  double R = 0;
  R = km_BetaN_a[m][n] * in(m, BetaN[n] * r) - im_BetaN_a[m][n] * kn(m, BetaN[n] * r);

  return R;
}

double Rossegger::Rmn1_(int m, int n, double r)
{
  //  Check input arguments for sanity...
  int error = 0;
  if (m < 0 || m >= NumberOfOrders) error = 1;
  if (n < 0 || n >= NumberOfOrders) error = 1;
  if (r < a || r > b) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments Rmn1(" << m << "," << n << "," << r << ")" << std::endl;
    ;
    return 0;
  }

  //  Calculate using the TMath functions from root.
  //  Rossegger Equation 5.32
  //         Rmn1(r) = Km(BetaN a)Im(BetaN r) - Im(BetaN a) Km(BetaN r)
  double R = 0;
  double BetaN_ = (n + 1) * pi / L;
  R = kn(m, BetaN_ * a) * in(m, BetaN_ * r) - in(m, BetaN_ * a) * kn(m, BetaN_ * r);

  return R;
}

double Rossegger::Rmn2(int m, int n, double r)
{
  //  Check input arguments for sanity...
  int error = 0;
  if (m < 0 || m >= NumberOfOrders) error = 1;
  if (n < 0 || n >= NumberOfOrders) error = 1;
  if (r < a || r > b) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments Rmn2(" << m << "," << n << "," << r << ")" << std::endl;
    ;
    return 0;
  }

  //  Calculate using the TMath functions from root.
  //  Rossegger Equation 5.33
  //         Rmn2(r) = Km(BetaN b)Im(BetaN r) - Im(BetaN b) Km(BetaN r)
  double R = 0;
  R = km_BetaN_b[m][n] * in(m, BetaN[n] * r) - im_BetaN_b[m][n] * kn(m, BetaN[n] * r);

  return R;
}

double Rossegger::Rmn2_(int m, int n, double r)
{
  //  Check input arguments for sanity...
  int error = 0;
  if (m < 0 || m >= NumberOfOrders) error = 1;
  if (n < 0 || n >= NumberOfOrders) error = 1;
  if (r < a || r > b) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments Rmn2(" << m << "," << n << "," << r << ")" << std::endl;
    ;
    return 0;
  }

  //  Calculate using the TMath functions from root.
  //  Rossegger Equation 5.33
  //         Rmn2(r) = Km(BetaN b)Im(BetaN r) - Im(BetaN b) Km(BetaN r)
  double R = 0;
  double BetaN_ = (n + 1) * pi / L;
  R = kn(m, BetaN_ * b) * in(m, BetaN_ * r) - in(m, BetaN_ * b) * kn(m, BetaN_ * r);

  return R;
}

double Rossegger::RPrime(int m, int n, double ref, double r)
{
  //  Check input arguments for sanity...
  int error = 0;
  if (m < 0 || m >= NumberOfOrders) error = 1;
  if (n < 0 || n >= NumberOfOrders) error = 1;
  if (ref < a || ref > b) error = 1;
  if (r < a || r > b) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments RPrime(" << m << "," << n << "," << ref << "," << r << ")" << std::endl;
    ;
    return 0;
  }

  double R = 0;
  //  Calculate using the TMath functions from root.
  //  Rossegger Equation 5.65
  //         Rmn2(ref,r) = BetaN/2* [   Km(BetaN ref) {Im-1(BetaN r) + Im+1(BetaN r)}
  //                                  - Im(BetaN ref) {Km-1(BetaN r) + Km+1(BetaN r)}  ]
  //  NOTE:  K-m(z) = Km(z) and I-m(z) = Im(z)... though boost handles negative orders.
  //
  // with: s -> ref,  t -> r,
  double BetaN_ = BetaN[n];
  double term1 = kn(m, BetaN_ * ref) * (in(m - 1, BetaN_ * r) + in(m + 1, BetaN_ * r));
  double term2 = in(m, BetaN_ * ref) * (kn(m - 1, BetaN_ * r) + kn(m + 1, BetaN_ * r));
  R = BetaN_ / 2.0 * (term1 + term2);

  return R;
}

double Rossegger::RPrime_(int m, int n, double ref, double r)
{
  //  Check input arguments for sanity...
  int error = 0;
  if (m < 0 || m >= NumberOfOrders) error = 1;
  if (n < 0 || n >= NumberOfOrders) error = 1;
  if (ref < a || ref > b) error = 1;
  if (r < a || r > b) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments RPrime(" << m << "," << n << "," << ref << "," << r << ")" << std::endl;
    ;
    return 0;
  }

  double R = 0;
  //  Rossegger Equation 5.65
  //         Rmn2(ref,r) = BetaN/2* [   Km(BetaN ref) {Im-1(BetaN r) + Im+1(BetaN r)}
  //                                  - Im(BetaN ref) {Km-1(BetaN r) + Km+1(BetaN r)}  ]
  //  NOTE:  K-m(z) = Km(z) and I-m(z) = Im(z)... though boost handles negative orders.
  //
  // with: s -> ref,  t -> r,
  double BetaN_ = (n + 1) * pi / L;
  double term1 = kn(m, BetaN_ * ref) * (in(m - 1, BetaN_ * r) + in(m + 1, BetaN_ * r));
  double term2 = in(m, BetaN_ * ref) * (kn(m - 1, BetaN_ * r) + kn(m + 1, BetaN_ * r));
  R = BetaN_ / 2.0 * (term1 + term2);

  return R;
}

double Rossegger::Rnk_for_zeroes(int n, double mu)
{
  //unlike Rossegger, we count 'k' and 'n' from zero.
  if (verbosity > 10) printf("Rnk_for_zeroes called with n=%d,mu=%f\n", n, mu);
  double betana = BetaN_a[n];
  double betanb = BetaN_b[n];
  //  Rossegger Equation 5.46
  //       Rnk(r) = Limu_nk (BetaN a) Kimu_nk (BetaN b) - Kimu_nk(BetaN a) Limu_nk (BetaN b)

  return limu(mu, betana) * kimu(mu, betanb) - kimu(mu, betana) * limu(mu, betanb);
}

double Rossegger::Rnk_for_zeroes_(int n, double mu)
{
  //unlike Rossegger, we count 'k' and 'n' from zero.
  if (verbosity > 10) printf("Rnk_for_zeroes called with n=%d,mu=%f\n", n, mu);
  double BetaN_ = (n + 1) * pi / L;  //this is defined in the paragraph before 5.46
  double betana = BetaN_ * a;
  double betanb = BetaN_ * b;
  //  Rossegger Equation 5.46
  //       Rnk(r) = Limu_nk (BetaN a) Kimu_nk (BetaN b) - Kimu_nk(BetaN a) Limu_nk (BetaN b)

  return limu(mu, betana) * kimu(mu, betanb) - kimu(mu, betana) * limu(mu, betanb);
}

double Rossegger::Rnk(int n, int k, double r)
{
  //  Check input arguments for sanity...
  int error = 0;
  if (n < 0 || n >= NumberOfOrders) error = 1;
  if (k < 0 || k >= NumberOfOrders) error = 1;
  if (r < a || r > b) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments Rnk(" << n << "," << k << "," << r << ")" << std::endl;
    ;
    return 0;
  }
  //  Rossegger Equation 5.45
  //       Rnk(r) = Limu_nk (BetaN a) Kimu_nk (BetaN r) - Kimu_nk(BetaN a) Limu_nk (BetaN r)

  return liMunk_BetaN_a[n][k] * kimu(Munk[n][k], BetaN[n] * r) - kiMunk_BetaN_a[n][k] * limu(Munk[n][k], BetaN[n] * r);
}

double Rossegger::Rnk_(int n, int k, double r)
{
  //  Check input arguments for sanity...
  int error = 0;
  if (n < 0 || n >= NumberOfOrders) error = 1;
  if (k < 0 || k >= NumberOfOrders) error = 1;
  if (r < a || r > b) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments Rnk(" << n << "," << k << "," << r << ")" << std::endl;
    ;
    return 0;
  }
  double BetaN_ = (n + 1) * pi / L;
  //  Rossegger Equation 5.45
  //       Rnk(r) = Limu_nk (BetaN a) Kimu_nk (BetaN r) - Kimu_nk(BetaN a) Limu_nk (BetaN r)

  return limu(Munk[n][k], BetaN_ * a) * kimu(Munk[n][k], BetaN_ * r) - kimu(Munk[n][k], BetaN_ * a) * limu(Munk[n][k], BetaN_ * r);
}

double Rossegger::Ez(double r, double phi, double z, double r1, double phi1, double z1)
{
  //rcc streamlined Ez
  int error = 0;
  if (r < a || r > b) error = 1;
  if (phi < 0 || phi > 2 * pi) error = 1;
  if (z < 0 || z > L) error = 1;
  if (r1 < a || r1 > b) error = 1;
  if (phi1 < 0 || phi1 > 2 * pi) error = 1;
  if (z1 < 0 || z1 > L) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments Ez(";
    std::cout << r << ",";
    std::cout << phi << ",";
    std::cout << z << ",";
    std::cout << r1 << ",";
    std::cout << phi1 << ",";
    std::cout << z1;
    std::cout << ")" << std::endl;
    ;
    return 0;
  }
  //Rossegger Equation 5.64
  double G = 0;
  for (int m = 0; m < NumberOfOrders; m++)
  {
    if (verbosity > 10) std::cout << std::endl
                                  << m;
    for (int n = 0; n < NumberOfOrders; n++)
    {
      if (verbosity > 10) std::cout << " " << n;
      double term = 1;  //unitless
      if (verbosity > 10) std::cout << " " << term;
      term *= (2 - ((m == 0) ? 1 : 0)) * cos(m * (phi - phi1));  //unitless
      if (verbosity > 10) std::cout << " " << term;
      term *= Rmn(m, n, r) * Rmn(m, n, r1) / N2mn[m][n];  //units of 1/[L]^2
      if (verbosity > 10) std::cout << " " << term;
      if (z < z1)
      {
        term *= cosh(Betamn[m][n] * z) * sinh(Betamn[m][n] * (L - z1)) / sinh_Betamn_L[m][n];  //unitless
      }
      else
      {
        term *= -cosh(Betamn[m][n] * (L - z)) * sinh(Betamn[m][n] * z1) / sinh_Betamn_L[m][n];  //unitless
      }
      if (verbosity > 10) std::cout << " " << term;
      G += term;
      if (verbosity > 10) std::cout << " " << term << " " << G << std::endl;
    }
  }

  G = G / (2.0 * pi);
  if (verbosity) std::cout << "Ez = " << G << std::endl;

  return G;
}

double Rossegger::Ez_(double r, double phi, double z, double r1, double phi1, double z1)
{
  //if(fByFile && fabs(r-r1)>MinimumDR && fabs(z-z1)>MinimumDZ) return ByFileEZ(r,phi,z,r1,phi1,z1);
  //  Check input arguments for sanity...
  int error = 0;
  if (r < a || r > b) error = 1;
  if (phi < 0 || phi > 2 * pi) error = 1;
  if (z < 0 || z > L) error = 1;
  if (r1 < a || r1 > b) error = 1;
  if (phi1 < 0 || phi1 > 2 * pi) error = 1;
  if (z1 < 0 || z1 > L) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments Ez(";
    std::cout << r << ",";
    std::cout << phi << ",";
    std::cout << z << ",";
    std::cout << r1 << ",";
    std::cout << phi1 << ",";
    std::cout << z1;
    std::cout << ")" << std::endl;
    ;
    return 0;
  }

  double G = 0;
  for (int m = 0; m < NumberOfOrders; m++)
  {
    if (verbosity > 10) std::cout << std::endl
                                  << m;
    for (int n = 0; n < NumberOfOrders; n++)
    {
      if (verbosity > 10) std::cout << " " << n;
      double term = 1 / (2.0 * pi);
      if (verbosity > 10) std::cout << " " << term;
      term *= (2 - ((m == 0) ? 1 : 0)) * cos(m * (phi - phi1));  //unitless
      if (verbosity > 10) std::cout << " " << term;
      term *= Rmn(m, n, r) * Rmn(m, n, r1) / N2mn[m][n];  //units of 1/[L]^2
      if (verbosity > 10) std::cout << " " << term;
      if (z < z1)
      {
        term *= cosh(Betamn[m][n] * z) * sinh(Betamn[m][n] * (L - z1)) / sinh(Betamn[m][n] * L);
      }
      else
      {
        term *= -cosh(Betamn[m][n] * (L - z)) * sinh(Betamn[m][n] * z1) / sinh(Betamn[m][n] * L);
        ;
      }
      if (verbosity > 10) std::cout << " " << term;
      G += term;
      if (verbosity > 10) std::cout << " " << term << " " << G << std::endl;
    }
  }
  if (verbosity) std::cout << "Ez = " << G << std::endl;

  return G;
}

double Rossegger::Er(double r, double phi, double z, double r1, double phi1, double z1)
{
  //rcc streamlined Er

  //as in Rossegger 5.65
  //field at r, phi, z due to unit charge at r1, phi1, z1;
  //if(fByFile && fabs(r-r1)>MinimumDR && fabs(z-z1)>MinimumDZ) return ByFileER(r,phi,z,r1,phi1,z1);
  //  Check input arguments for sanity...
  int error = 0;
  if (r < a || r > b) error = 1;
  if (phi < 0 || phi > 2 * pi) error = 1;
  if (z < 0 || z > L) error = 1;
  if (r1 < a || r1 > b) error = 1;
  if (phi1 < 0 || phi1 > 2 * pi) error = 1;
  if (z1 < 0 || z1 > L) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments Er(";
    std::cout << r << ",";
    std::cout << phi << ",";
    std::cout << z << ",";
    std::cout << r1 << ",";
    std::cout << phi1 << ",";
    std::cout << z1;
    std::cout << ")" << std::endl;
    ;
    return 0;
  }

  double part = 0;
  double G = 0;
  for (int m = 0; m < NumberOfOrders; m++)
  {
    for (int n = 0; n < NumberOfOrders; n++)
    {
      double term = 1;
      part = (2 - ((m == 0) ? 1 : 0)) * cos(m * (phi - phi1));  //unitless
      if (verbosity > 10) printf("(2 - ((m==0)?1:0))*cos(m*(phi-phi1)); = %f\n", part);
      term *= part;
      part = sin(BetaN[n] * z) * sin(BetaN[n] * z1);  //unitless
      if (verbosity > 10) printf("sin(BetaN[n]*z)*sin(BetaN[n]*z1); = %f\n", part);
      term *= part;

      if (r < r1)
      {
        term *= RPrime(m, n, a, r) * Rmn2(m, n, r1);  //units of 1/[L]
      }
      else
      {
        term *= Rmn1(m, n, r1) * RPrime(m, n, b, r);  //units of 1/[L]
      }
      term /= bessel_denominator[m][n];  //unitless
      G += term;
    }
  }

  G = G / (L * pi);  //units of 1/[L] -- net is 1/[L]^2
  if (verbosity) std::cout << "Er = " << G << std::endl;

  return G;
}

double Rossegger::Er_(double r, double phi, double z, double r1, double phi1, double z1)
{
  //doesn't take advantage of precalcs.
  //as in Rossegger 5.65
  //field at r, phi, z due to unit charge at r1, phi1, z1;
  //if(fByFile && fabs(r-r1)>MinimumDR && fabs(z-z1)>MinimumDZ) return ByFileER(r,phi,z,r1,phi1,z1);
  //  Check input arguments for sanity...
  int error = 0;
  if (r < a || r > b) error = 1;
  if (phi < 0 || phi > 2 * pi) error = 1;
  if (z < 0 || z > L) error = 1;
  if (r1 < a || r1 > b) error = 1;
  if (phi1 < 0 || phi1 > 2 * pi) error = 1;
  if (z1 < 0 || z1 > L) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments Er(";
    std::cout << r << ",";
    std::cout << phi << ",";
    std::cout << z << ",";
    std::cout << r1 << ",";
    std::cout << phi1 << ",";
    std::cout << z1;
    std::cout << ")" << std::endl;
    ;
    return 0;
  }

  double G = 0;
  for (int m = 0; m < NumberOfOrders; m++)
  {
    for (int n = 0; n < NumberOfOrders; n++)
    {
      double term = 1 / (L * pi);
      term *= (2 - ((m == 0) ? 1 : 0)) * cos(m * (phi - phi1));
      double BetaN_ = (n + 1) * pi / L;
      term *= sin(BetaN_ * z) * sin(BetaN_ * z1);
      if (r < r1)
      {
        term *= RPrime_(m, n, a, r) * Rmn2_(m, n, r1);
      }
      else
      {
        term *= Rmn1_(m, n, r1) * RPrime_(m, n, b, r);
      }

      term /= TMath::BesselI(m, BetaN_ * a) * TMath::BesselK(m, BetaN_ * b) - TMath::BesselI(m, BetaN_ * b) * TMath::BesselK(m, BetaN_ * a);

      G += term;
    }
  }

  if (verbosity) std::cout << "Er = " << G << std::endl;

  return G;
}

double Rossegger::Ephi(double r, double phi, double z, double r1, double phi1, double z1)
{
  //rcc streamlined Ephi term
  //compute field at rphiz from charge at r1phi1z1
  //  Check input arguments for sanity...
  int error = 0;
  if (r < a || r > b) error = 1;
  if (phi < 0 || phi > 2 * pi) error = 1;
  if (z < 0 || z > L) error = 1;
  if (r1 < a || r1 > b) error = 1;
  if (phi1 < 0 || phi1 > 2 * pi) error = 1;
  if (z1 < 0 || z1 > L) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments Ephi(";
    std::cout << r << ",";
    std::cout << phi << ",";
    std::cout << z << ",";
    std::cout << r1 << ",";
    std::cout << phi1 << ",";
    std::cout << z1;
    std::cout << ")" << std::endl;
    ;
    return 0;
  }

  double G = 0;
  //Rossegger Eqn. 5.66:
  for (int k = 0; k < NumberOfOrders; k++)
  {
    for (int n = 0; n < NumberOfOrders; n++)
    {
      double term = 1;
      term *= sin(BetaN[n] * z) * sin(BetaN[n] * z1);     //unitless
      term *= Rnk(n, k, r) * Rnk(n, k, r1) / N2nk[n][k];  //unitless?

      //the derivative of cosh(munk(pi-|phi-phi1|)
      if (phi > phi1)
      {
        term *= -sinh(Munk[n][k] * (pi - (phi - phi1)));  //unitless
      }
      else
      {
        term *= sinh(Munk[n][k] * (pi - (phi1 - phi)));  //unitless
      }
      term *= 1 / sinh_pi_Munk[n][k];  //unitless
      G += term;
    }
  }

  G = G / (L * r);  //units of 1/[L]^2.  r comes from the phi term in cylindrical gradient expression.
  if (verbosity) std::cout << "Ephi = " << G << std::endl;

  return G;
}

double Rossegger::Ephi_(double r, double phi, double z, double r1, double phi1, double z1)
{
  //compute field at rphiz from charge at r1phi1z1
  //  Check input arguments for sanity...
  int error = 0;
  if (r < a || r > b) error = 1;
  if (phi < 0 || phi > 2 * pi) error = 1;
  if (z < 0 || z > L) error = 1;
  if (r1 < a || r1 > b) error = 1;
  if (phi1 < 0 || phi1 > 2 * pi) error = 1;
  if (z1 < 0 || z1 > L) error = 1;
  if (error)
  {
    std::cout << "Invalid arguments Ephi(";
    std::cout << r << ",";
    std::cout << phi << ",";
    std::cout << z << ",";
    std::cout << r1 << ",";
    std::cout << phi1 << ",";
    std::cout << z1;
    std::cout << ")" << std::endl;
    ;
    return 0;
  }

  verbosity = 1;
  double G = 0;
  //Rossegger Eqn. 5.66:
  for (int k = 0; k < NumberOfOrders; k++)  //off by one from Rossegger convention!
  {
    if (verbosity) std::cout << "\nk=" << k;
    for (int n = 0; n < NumberOfOrders; n++)  //off by one from Rossegger convention!
    {
      if (verbosity) std::cout << " n=" << n;
      double BetaN_ = (n + 1) * pi / L;
      double term = 1 / (L * r);
      if (verbosity) std::cout << " 1/L=" << term;
      term *= sin(BetaN_ * z) * sin(BetaN_ * z1);
      if (verbosity) std::cout << " *sinsin=" << term;
      term *= Rnk_(n, k, r) * Rnk_(n, k, r1);
      if (verbosity) std::cout << " *rnkrnk=" << term;
      term /= N2nk[n][k];
      if (verbosity) std::cout << " */nnknnk=" << term;

      //the derivative of cosh(munk(pi-|phi-phi1|)
      if (phi > phi1)
      {
        term *= -sinh(Munk[n][k] * (pi - (phi - phi1)));
        //term *=  Munk[n][k]*sinh(Munk[n][k]*pi*(phi1-phi));
        //this originally has a factor of Munk in front, but that cancels with one in the denominator
      }
      else
      {
        term *= sinh(Munk[n][k] * (pi - (phi1 - phi)));
        //term *=  -Munk[n][k]*sinh(Munk[n][k]*pi*(phi-phi1));
        //this originally has a factor of Munk in front, but that cancels with one in the denominator
      }
      if (verbosity) std::cout << " *sinh(mu*pi-phi-phi)=" << term;
      term *= 1 / (sinh(pi * Munk[n][k]));
      //term *= 1/(Munk[n][k]*sinh(pi*Munk[n][k]));
      //this originally has a factor of Munk in front, but that cancels with one in the numerator
      G += term;
      if (verbosity) std::cout << "  /sinh=" << term << " G=" << G << std::endl;
    }
  }
  if (verbosity) std::cout << "Ephi = " << G << std::endl;
  verbosity = 0;

  return G;
}

void Rossegger::SaveZeroes(const char *destfile)
{
  TFile *output = TFile::Open(destfile, "RECREATE");
  output->cd();

  TTree *tInfo = new TTree("info", "Mu[n][k] values");
  int ord = NumberOfOrders;
  tInfo->Branch("order", &ord);
  tInfo->Branch("epsilon", &epsilon);
  tInfo->Fill();

  int n, k, m;
  double munk, betamn;
  double n2nk, n2mn;
  TTree *tmunk = new TTree("munk", "Mu[n][k] values");
  tmunk->Branch("n", &n);
  tmunk->Branch("k", &k);
  tmunk->Branch("munk", &munk);
  tmunk->Branch("n2nk", &n2nk);
  for (n = 0; n < ord; n++)
  {
    for (k = 0; k < ord; k++)
    {
      munk = Munk[n][k];
      n2nk = N2nk[n][k];
      tmunk->Fill();
    }
  }

  TTree *tbetamn = new TTree("betamn", "Beta[m][n] values");
  tbetamn->Branch("m", &m);
  tbetamn->Branch("n", &n);
  tbetamn->Branch("betamn", &betamn);
  tbetamn->Branch("n2mn", &n2mn);
  for (m = 0; m < ord; m++)
  {
    for (n = 0; n < ord; n++)
    {
      betamn = Betamn[m][n];
      n2mn = N2mn[m][n];
      tbetamn->Fill();
    }
  }

  tInfo->Write();
  tmunk->Write();
  tbetamn->Write();
  //output->Write();
  output->Close();
  return;
}

void Rossegger::LoadZeroes(const char *destfile)
{
  TFile *f = TFile::Open(destfile, "READ");
  printf("reading rossegger zeroes from %s\n", destfile);
  TTree *tInfo = (TTree *) (f->Get("info"));
  int ord;
  tInfo->SetBranchAddress("order", &ord);
  tInfo->SetBranchAddress("epsilon", &epsilon);
  tInfo->GetEntry(0);
  printf("order=%d,epsilon=%f\n", ord, epsilon);

  int n, k, m;
  double munk, betamn;
  double n2nk, n2mn;
  TTree *tmunk = (TTree *) (f->Get("munk"));
  tmunk->SetBranchAddress("n", &n);
  tmunk->SetBranchAddress("k", &k);
  tmunk->SetBranchAddress("munk", &munk);
  tmunk->SetBranchAddress("n2nk", &n2nk);
  for (int i = 0; i < tmunk->GetEntries(); i++)
  {
    tmunk->GetEntry(i);
    Munk[n][k] = munk;
    N2nk[n][k] = n2nk;
  }

  TTree *tbetamn = (TTree *) (f->Get("betamn"));
  tbetamn->SetBranchAddress("m", &m);
  tbetamn->SetBranchAddress("n", &n);
  tbetamn->SetBranchAddress("betamn", &betamn);
  tbetamn->SetBranchAddress("n2mn", &n2mn);
  for (int i = 0; i < tbetamn->GetEntries(); i++)
  {
    tbetamn->GetEntry(i);
    Betamn[m][n] = betamn;
    N2mn[m][n] = n2mn;
  }

  f->Close();
  return;
}
