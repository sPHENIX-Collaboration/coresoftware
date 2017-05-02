#ifndef __LAPLACESOLUTION_H__
#define __LAPLACESOLUTION_H__

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

#define NumberOfOrders 15  // Convergence problems after 15; Rossegger used 30
#include <string>

class LaplaceSolution
{
 public:
  LaplaceSolution(std::string filename);
  LaplaceSolution(double a=30, double b=80, double L=80);
  virtual ~LaplaceSolution() {}

  void Verbosity(int v) {verbosity=v;}
  double Rmn (int m, int n, double r);  //Rmn function from Rossegger
  double Rmn1(int m, int n, double r);  //Rmn1 function from Rossegger
  double Rmn2(int m, int n, double r);  //Rmn2 function from Rossegger
  double RPrime(int m, int n, double a, double r);  // RPrime function from Rossegger
  
  double Rnk(int n, int k, double r);  //Rnk function from Rossegger

  double Ez  (double r, double phi, double z, double r1, double phi1, double z1);
  double Er  (double r, double phi, double z, double r1, double phi1, double z1);
  double Ephi(double r, double phi, double z, double r1, double phi1, double z1);

 protected:
  bool fByFile;
  double a,b,L;  //  InnerRadius, OuterRadius, Length of 1/2 the TPC.
  int verbosity;
  double pi;

  void FindBetamn(double epsilon);  // Routine used to fill the Betamn array with resolution epsilon...
  void FindMunk(double epsilon);    // Routine used to fill the Munk array with resolution epsilon...

  double Betamn[NumberOfOrders][NumberOfOrders];  //  Betamn array from Rossegger
  double N2mn[NumberOfOrders][NumberOfOrders];    //  N2mn array from Rossegger
  double Munk[NumberOfOrders][NumberOfOrders];    //  Munk array from Rossegger
  double ByFileER(double r, double phi, double z, double r1, double phi1, double z1);
  double ByFileEZ(double r, double phi, double z, double r1, double phi1, double z1);


};

#endif /* __LAPLACESOLUTION_H__ */



