#include "Shifter.h"

#include <TFile.h>
#include <TH3.h>
#include <TString.h>
#include <TVector3.h>

#include <cmath>
#include <cstdlib>  // for getenv

Shifter::Shifter(const std::string &truthfilename, const std::string &correctionfilename)
{
  //load a 'truth' distortion map and, optionally, a map of a measured correction to those distortions
  //this code is currently set up to load a particular correction map that doesn't have distortions
  // in X,Y, and Z components, but rather only in R, R*Phi, and Z components.

  //single event distortion file
  if (!truthfilename.empty())
  {
    forward = TFile::Open(truthfilename.c_str(), "READ");
    if (forward)
    {
      forward->GetObject("hIntDistortionX", hX);
      forward->GetObject("hIntDistortionY", hY);
      forward->GetObject("hIntDistortionZ", hZ);

      //not strictly needed, but handy:
      forward->GetObject("hIntDistortionR", hR);
      forward->GetObject("hIntDistortionP", hPhi);
    }
  }
  if (hX && hY && hZ)
  {
    hasTruth = true;
  }

  //single event distortion file
  if (!correctionfilename.empty())
  {
    //average=TFile::Open(correctionfilename,"READ");
    // hardcoded????????
    std::string correction_filename = std::string(getenv("CALIBRATIONROOT")) + "/distortion_maps/Distortions_full_realistic_micromegas_all-coarse.root";
    average = TFile::Open(correction_filename.c_str(), "READ");
    if (average)
    {
      average->GetObject("hIntDistortionX", hXave);
      average->GetObject("hIntDistortionY", hYave);
      average->GetObject("hIntDistortionZ", hZave);

      average->GetObject("hIntDistortionR", hRave);
      average->GetObject("hIntDistortionP", hPhiave);
    }
  }
  if (hXave && hYave && hZave)
  {
    hasCorrection = true;
  }
}

TVector3 Shifter::ShiftForward(const TVector3 &position)
{
  double x, y, z, xshift, yshift, zshift;
  //const double mm = 1.0;
  //const double cm = 10.0;
  TVector3 shiftposition;

  x = position.X();
  y = position.Y();
  z = position.Z();

  double r = position.Perp();
  double phi = position.Phi();
  if (position.Phi() < 0.0)
  {
    phi = position.Phi() + 2.0 * M_PI;
  }

  //distort coordinate of stripe
  xshift = 0;
  yshift = 0;
  zshift = 0;
  if (hasTruth)
  {
    xshift = hX->Interpolate(phi, r, z);
    yshift = hY->Interpolate(phi, r, z);
    zshift = hZ->Interpolate(phi, r, z);
  }

  //remove average distortion
  if (hasCorrection)
  {
    double raveshift = hRave->Interpolate(phi, r, z);
    double paveshift = hPhiave->Interpolate(phi, r, z);  //hugo confirms the units are cm
    double cosphi = cos(phi);
    double sinphi = sin(phi);
    xshift -= raveshift * cosphi - paveshift * sinphi;
    yshift -= raveshift * sinphi + paveshift * cosphi;

    zshift -= hZave->Interpolate(phi, r, z);
  }

  TVector3 forwardshift(x + xshift, y + yshift, z + zshift);

  return forwardshift;
}

TVector3 Shifter::ShiftBack(const TVector3 &forwardshift)
{
  double x, y, z;
  // const double mm = 1.0;
  //const double cm = 10.0;
  TVector3 shiftposition;

  x = forwardshift.X();
  y = forwardshift.Y();
  z = forwardshift.Z();

  double rforward = forwardshift.Perp();
  double phiforward = forwardshift.Phi();
  if (forwardshift.Phi() < 0.0)
  {
    phiforward += 2.0 * M_PI;
  }

  double xshiftback = -1 * hXBack->Interpolate(phiforward, rforward, z);
  double yshiftback = -1 * hYBack->Interpolate(phiforward, rforward, z);
  double zshiftback = -1 * hZBack->Interpolate(phiforward, rforward, z);

  shiftposition.SetXYZ(x + xshiftback, y + yshiftback, z + zshiftback);

  return shiftposition;
}

TVector3 Shifter::Shift(const TVector3 &position)
{
  return ShiftBack(ShiftForward(position));
}
