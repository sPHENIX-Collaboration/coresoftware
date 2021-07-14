#ifndef PHG4TPCCENTRALMEMBRANE_H
#define PHG4TPCCENTRALMEMBRANE_H

#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"

//from phg4tpcsteppingaction.cc
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
//R__LOAD_LIBRARY(libphg4hit.so)


// all distances in mm, all angles in rad
// class that generates stripes and dummy hit coordinates
// stripes have width of one mm, length of one pad width, and are centered in middle of sector gaps

using namespace std;

class PHG4TpcCentralMembrane {
public:
  PHG4TpcCentralMembrane(); //default constructor
  int getSearchResult(double xcheck, double ycheck); // check if coords are in a stripe
  int getStripeID(double xcheck, double ycheck);
 
  int fullID;
  double begin_CM, end_CM; // inner and outer radii of central membrane
  
  vector<PHG4Hitv1*> PHG4Hits;
  vector<PHG4Hitv1*> BotVertices;
  vector<PHG4Hitv1*> TopVertices;
  
private:
  static const int nRadii = 8;
  static const int nStripes_R1 = 6;
  static const int nStripes_R2 = 8;
  static const int nStripes_R3 = 12;

  const double mm = 1.0;
  const double cm = 10.0;
  
  int nPads_R1;
  int nPads_R2;
  int nPads_R3;
  
  double padfrac_R1;
  double padfrac_R2;
  double padfrac_R3;
  //double str_width; // width of a stripe
  double arc_r; // radius of arc on end of a stripe
  double R1_e[nRadii], R1[nRadii], R2[nRadii], R3[nRadii];
  
  double str_width_R1_e[nStripes_R1][nRadii];
  double str_width_R1[nStripes_R1][nRadii];
  double str_width_R2[nStripes_R2][nRadii];
  double str_width_R3[nStripes_R3][nRadii];
  double widthmod_R1_e[nRadii];
  double widthmod_R1[nRadii];
  double widthmod_R2[nRadii];
  double widthmod_R3[nRadii];
 
  double spacing_R1_e[nRadii], spacing_R1[nRadii], spacing_R2[nRadii], spacing_R3[nRadii];
  
  //bottom left - 1a
  double x1a_R1_e[nStripes_R1][nRadii], y1a_R1_e[nStripes_R1][nRadii];
  double x1a_R1[nStripes_R1][nRadii], y1a_R1[nStripes_R1][nRadii];
  double x1a_R2[nStripes_R2][nRadii], y1a_R2[nStripes_R2][nRadii];
  double x1a_R3[nStripes_R3][nRadii], y1a_R3[nStripes_R3][nRadii];
  
  //bottom right - 1b
  double x1b_R1_e[nStripes_R1][nRadii], y1b_R1_e[nStripes_R1][nRadii];
  double x1b_R1[nStripes_R1][nRadii], y1b_R1[nStripes_R1][nRadii];
  double x1b_R2[nStripes_R2][nRadii], y1b_R2[nStripes_R2][nRadii];
  double x1b_R3[nStripes_R3][nRadii], y1b_R3[nStripes_R3][nRadii];
  
  //top left - 2a
  double x2a_R1_e[nStripes_R1][nRadii], y2a_R1_e[nStripes_R1][nRadii];
  double x2a_R1[nStripes_R1][nRadii], y2a_R1[nStripes_R1][nRadii];
  double x2a_R2[nStripes_R2][nRadii], y2a_R2[nStripes_R2][nRadii];
  double x2a_R3[nStripes_R3][nRadii], y2a_R3[nStripes_R3][nRadii];
  
  //top right - 2b
  double x2b_R1_e[nStripes_R1][nRadii], y2b_R1_e[nStripes_R1][nRadii];
  double x2b_R1[nStripes_R1][nRadii], y2b_R1[nStripes_R1][nRadii];
  double x2b_R2[nStripes_R2][nRadii], y2b_R2[nStripes_R2][nRadii];
  double x2b_R3[nStripes_R3][nRadii], y2b_R3[nStripes_R3][nRadii];
  
  //left midpoint - 3a
  double x3a_R1_e[nStripes_R1][nRadii], y3a_R1_e[nStripes_R1][nRadii];
  double x3a_R1[nStripes_R1][nRadii], y3a_R1[nStripes_R1][nRadii];
  double x3a_R2[nStripes_R2][nRadii], y3a_R2[nStripes_R2][nRadii];
  double x3a_R3[nStripes_R3][nRadii], y3a_R3[nStripes_R3][nRadii];
  
  //right midpoint - 3b
  double x3b_R1_e[nStripes_R1][nRadii], y3b_R1_e[nStripes_R1][nRadii];
  double x3b_R1[nStripes_R1][nRadii], y3b_R1[nStripes_R1][nRadii];
  double x3b_R2[nStripes_R2][nRadii], y3b_R2[nStripes_R2][nRadii];
  double x3b_R3[nStripes_R3][nRadii], y3b_R3[nStripes_R3][nRadii];
  
  //Check which stripes get removed
  int nGoodStripes_R1_e[nRadii];
  int nGoodStripes_R1[nRadii];
  int nGoodStripes_R2[nRadii];
  int nGoodStripes_R3[nRadii];
  int keepThisAndAfter[nRadii]; //min stripe index
  int keepUntil_R1_e[nRadii]; //max stripe index
  int keepUntil_R1[nRadii];
  int keepUntil_R2[nRadii];
  int keepUntil_R3[nRadii];
  int nStripesIn_R1_e[nRadii];
  int nStripesIn_R1[nRadii];
  int nStripesIn_R2[nRadii];
  int nStripesIn_R3[nRadii];
  int nStripesBefore_R1_e[nRadii];
  int nStripesBefore_R1[nRadii];
  int nStripesBefore_R2[nRadii];
  int nStripesBefore_R3[nRadii];
  int result;

  double botvert, topvert;
  
  int nStripesPerPetal;
  int nPetals;
  int nTotStripes;
  int stripeID;
  
  int nElectrons;
  
  void CalculateVertices(int nStripes, int nPads, double R[], double spacing[], double x1a[][nRadii], double y1a[][nRadii], double x1b[][nRadii], double y1b[][nRadii], double x2a[][nRadii], double y2a[][nRadii], double x2b[][nRadii], double y2b[][nRadii], double x3a[][nRadii], double y3a[][nRadii], double x3b[][nRadii], double y3b[][nRadii], double padfrac, double str_width[][nRadii], double widthmod[], int nGoodStripes[], int keepUntil[], int nStripesIn[], int nStripesBefore[]);

  PHG4Hitv1* GetBotVerticesFromStripe(int moduleID, int radiusID, int stripeID);
  PHG4Hitv1* GetTopVerticesFromStripe(int moduleID, int radiusID, int stripeID);
  
  int SearchModule(int nStripes, double x1a[][nRadii], double x1b[][nRadii], double x2a[][nRadii], double x2b[][nRadii], double y1a[][nRadii], double y1b[][nRadii], double y2a[][nRadii], double y2b[][nRadii], double x3a[][nRadii], double y3a[][nRadii], double x3b[][nRadii], double y3b[][nRadii], double x, double y, int nGoodStripes[]);
  
  PHG4Hitv1* GetPHG4HitFromStripe(int petalID, int moduleID, int radiusID, int stripeID, int nElectrons);
};


#endif
