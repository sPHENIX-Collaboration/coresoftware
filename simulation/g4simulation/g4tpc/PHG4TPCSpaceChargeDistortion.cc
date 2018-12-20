// $Id: $                                                                                             

/*!
 * \file PHG4TPCSpaceChargeDistortion.cc
 * \brief From TKH's SpaceChargeDistortion
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4TPCSpaceChargeDistortion.h"

#include <gsl/gsl_randist.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TFile.h>

#include <iostream>
#include <cassert>

using namespace std;

PHG4TPCSpaceChargeDistortion::PHG4TPCSpaceChargeDistortion(
    const std::string & distortion_map_file, int verbose) :
    PHG4TPCDistortion(verbose)
{
  TFile file(distortion_map_file.c_str());

  if (not file.IsOpen())
    {
      cout
          << "PHG4TPCSpaceChargeDistortion::PHG4TPCSpaceChargeDistortion - Fatal Error - "
          << "Failed to open distortion file " << distortion_map_file << endl;

      exit(13);
    }

  TH3F *rDistortion = dynamic_cast<TH3F *>(file.Get("mapDeltaR"));
  if (not rDistortion)
    {
      cout
          << "PHG4TPCSpaceChargeDistortion::PHG4TPCSpaceChargeDistortion - Fatal Error - "
          << "Failed to find TH3F mapDeltaR in distortion file "
          << distortion_map_file << endl;

      exit(13);
    }
  rDistortion->GetYaxis()->SetRange(1, 1);
  rDistortion2 = dynamic_cast<TH2D *>(rDistortion->Project3D("xz"));
  assert(rDistortion2);
  rDistortion2->SetDirectory(NULL); // make sure to detach from the current TFile

  TH3F *rPhiDistortion = dynamic_cast<TH3F *>(file.Get("mapRDeltaPHI"));
  if (not rPhiDistortion)
    {
      cout
          << "PHG4TPCSpaceChargeDistortion::PHG4TPCSpaceChargeDistortion - Fatal Error - "
          << "Failed to find TH3F mapRDeltaPHI in distortion file "
          << distortion_map_file << endl;

      exit(13);
    }
  rPhiDistortion->GetYaxis()->SetRange(1, 1);
  rPhiDistortion2 = dynamic_cast<TH2D *>(rPhiDistortion->Project3D("xz"));
  assert(rPhiDistortion2);
  rPhiDistortion2->SetDirectory(NULL); // make sure to detach from the current TFile

  //  Default to ALICE values...
  precisionFactor = 0.001;
  accuracyFactor = 0.005;

  if (verbosity > 1)
    {
      cout << "PHG4TPCSpaceChargeDistortion::PHG4TPCSpaceChargeDistortion - "
          << endl;

      double lowX = rDistortion->GetXaxis()->GetBinCenter(1);
      double highX = rDistortion->GetXaxis()->GetBinCenter(
          rDistortion->GetXaxis()->GetNbins());
      double lowY = rDistortion->GetYaxis()->GetBinCenter(1);
      double highY = rDistortion->GetYaxis()->GetBinCenter(
          rDistortion->GetYaxis()->GetNbins());
      double lowZ = rDistortion->GetZaxis()->GetBinCenter(1);
      double highZ = rDistortion->GetZaxis()->GetBinCenter(
          rDistortion->GetZaxis()->GetNbins());

      cout << " " << lowX << " " << highX;
      cout << " " << lowY << " " << highY;
      cout << " " << lowZ << " " << highZ;
      cout << endl;

      double lowX2 = rDistortion2->GetXaxis()->GetBinCenter(1);
      double highX2 = rDistortion2->GetXaxis()->GetBinCenter(
          rDistortion2->GetXaxis()->GetNbins());
      double lowY2 = rDistortion2->GetYaxis()->GetBinCenter(1);
      double highY2 = rDistortion2->GetYaxis()->GetBinCenter(
          rDistortion2->GetYaxis()->GetNbins());

      cout << " " << lowX2 << " " << highX2;
      cout << " " << lowY2 << " " << highY2;
      cout << endl;
    }

  if (verbosity > 1)
    {
      cout << "PHG4TPCSpaceChargeDistortion::PHG4TPCSpaceChargeDistortion - "
          << endl;

      double lowX = rPhiDistortion->GetXaxis()->GetBinCenter(1);
      double highX = rPhiDistortion->GetXaxis()->GetBinCenter(
          rPhiDistortion->GetXaxis()->GetNbins());
      double lowY = rPhiDistortion->GetYaxis()->GetBinCenter(1);
      double highY = rPhiDistortion->GetYaxis()->GetBinCenter(
          rPhiDistortion->GetYaxis()->GetNbins());
      double lowZ = rPhiDistortion->GetZaxis()->GetBinCenter(1);
      double highZ = rPhiDistortion->GetZaxis()->GetBinCenter(
          rPhiDistortion->GetZaxis()->GetNbins());

      cout << " " << lowX << " " << highX;
      cout << " " << lowY << " " << highY;
      cout << " " << lowZ << " " << highZ;
      cout << endl;

      double lowX2 = rPhiDistortion2->GetXaxis()->GetBinCenter(1);
      double highX2 = rPhiDistortion2->GetXaxis()->GetBinCenter(
          rPhiDistortion2->GetXaxis()->GetNbins());
      double lowY2 = rPhiDistortion2->GetYaxis()->GetBinCenter(1);
      double highY2 = rPhiDistortion2->GetYaxis()->GetBinCenter(
          rPhiDistortion2->GetYaxis()->GetNbins());

      cout << " " << lowX2 << " " << highX2;
      cout << " " << lowY2 << " " << highY2;
      cout << endl;
    }

}

PHG4TPCSpaceChargeDistortion::~PHG4TPCSpaceChargeDistortion()
{
  if (rDistortion2)
    delete rDistortion2;
  if (rPhiDistortion2)
    delete rPhiDistortion2;
}

double
PHG4TPCSpaceChargeDistortion::get_r_distortion(double r, double phi, double z)
{
  //  Calculations ONLY for minus z;
  if (z > 0)
    z = -z;

  //double dist = rDistortion->Interpolate(r,lowY,z);
  double dist = rDistortion2->Interpolate(z, r);
  double dist2 = accuracyFactor * dist
      + gsl_ran_gaussian(RandomGenerator, precisionFactor * dist);
  if (verbosity > 0)
    {
      cout <<"PHG4TPCSpaceChargeDistortion::get_r_distortion - input"
          <<" r = "<<r
          <<" phi = "<<phi
          <<" z = "<<z<<endl;
      cout << " Uncorrected R Distortion:" << dist;
      cout << " Corrected R Distortion:" << dist2;
      cout << endl;
    }

  return dist2;
}

double
PHG4TPCSpaceChargeDistortion::get_rphi_distortion(double r, double phi,
    double z)
{
  //  Need to have this histogram in the file...
  //  Calculations ONLY for minus z;
  if (z > 0)
    z = -z;

  //double dist = rPhiDistortion->Interpolate(r,lowY,z);
  double dist = rPhiDistortion2->Interpolate(z, r);
  double dist2 = accuracyFactor * dist
      + gsl_ran_gaussian(RandomGenerator, precisionFactor * dist);
  if (verbosity > 0)
    {
      cout <<"PHG4TPCSpaceChargeDistortion::get_rphi_distortion - input"
          <<" r = "<<r
          <<" phi = "<<phi
          <<" z = "<<z<<endl;
      cout << " Uncorrected rPHI Distortion:" << dist;
      cout << " Corrected rPHI Distortion:" << dist2;
      cout << endl;
    }

  return dist2;
}

