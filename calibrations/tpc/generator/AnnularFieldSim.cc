#include "AnnularFieldSim.h"

#include "AnalyticFieldModel.h"
#include "ChargeMapReader.h"
#include "MultiArray.h"  //for TH3 alternative
#include "Rossegger.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVector3.h>

#include <boost/format.hpp>

#include <cassert>  // for assert
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>

#define ALMOST_ZERO 0.00001

AnnularFieldSim::AnnularFieldSim(float in_innerRadius, float in_outerRadius, float in_outerZ,
                                 int r, int roi_r0, int roi_r1, int /*in_rLowSpacing*/, int /*in_rHighSize*/,
                                 int phi, int roi_phi0, int roi_phi1, int /*in_phiLowSpacing*/, int /*in_phiHighSize*/,
                                 int z, int roi_z0, int roi_z1, int /*in_zLowSpacing*/, int /*in_zHighSize*/,
                                 float vdr, LookupCase in_lookupCase, ChargeCase in_chargeCase)
{
  // check well-ordering:
  if (roi_r0 >= r || roi_r1 > r || roi_r0 >= roi_r1)
  {
    exit(1);
  }
  if (roi_phi0 >= phi || roi_phi1 > phi || roi_phi0 >= roi_phi1)
  {
    std::cout << "phi roi is out of range or spans the wrap-around.  Please spare me that math." << std::endl;
    exit(1);
  }
  if (roi_z0 >= z || roi_z1 > z || roi_z0 >= roi_z1)
  {
    exit(1);
  }

  hasTwin = false;  // we never have a twin to start with.
  green_shift = 0;  // and so we don't shift our greens functions unless told to.

  std::cout << boost::str(boost::format("AnnularFieldSim::AnnularFieldSim with (%dx%dx%d) grid") % r % phi % z) << std::endl;
  std::cout << boost::str(boost::format("units m=%1.2E, cm=%1.2E, mm=%1.2E, um=%1.2E,") % m % cm % mm % um) << std::endl;
  std::cout << boost::str(boost::format("units s=%1.2E, us=%1.2E, ns=%1.2E,") % s % us % ns) << std::endl;
  std::cout << boost::str(boost::format("units C=%1.2E, nC=%1.2E, fC=%1.2E, ") % C % nC % fC) << std::endl;
  std::cout << boost::str(boost::format("units Tesla=%1.2E, kGauss=%1.2E") % Tesla % kGauss) << std::endl;

  // debug defaults:
  //
  debug_printActionEveryN = -1;
  debug_printCounter = 0;
  debug_distortionScale.SetXYZ(0, 0, 0);
  debug_npercent = 5;

  // internal states used for the debug flag
  // goes with: if(debugFlag()) print_need_cout("%d: blah\n",__LINE__);

  // load constants of motion, dimensions, etc:
  Enominal = 400 * V / cm;  // v/cm
  Bnominal = 1.4 * Tesla;   // Tesla
  vdrift = vdr * cm / s;    // cm/s
  langevin_T1 = 1.0;
  langevin_T2 = 1.0;
  omegatau_nominal = -10 * (Bnominal / kGauss) * (vdrift / (cm / us)) / (Enominal / (V / cm));  // to match to the familiar formula

  rmin = in_innerRadius * cm;
  rmax = in_outerRadius * cm;
  zmin = 0;
  zmax = in_outerZ * cm;
  if (zmin > zmax)
  {
    zmin = zmax;
    zmax = 0;  // for now, we continue to assume one end of the TPC is at z=0.
  }
  zero_vector.SetXYZ(0, 0, 0);

  // define the size of the volume:
  dim.SetXYZ(1, 0, 0);
  dim.SetPerp(rmax - rmin);
  dim.SetPhi(0);
  phispan = 2 * M_PI;
  dim.SetZ(zmax - zmin);

  // set the green's functions model:
  // UseFreeSpaceGreens();
  // blah
  green = nullptr;

  // load parameters of the whole-volume tiling
  nr = r;
  nphi = phi;
  nz = z;  // number of fundamental bins (f-bins) in each direction
  std::cout << boost::str(boost::format("AnnularFieldSim::AnnularFieldSim set variables nr=%d, nphi=%d, nz=%d") % nr % nphi % nz) << std::endl;
  // and set the default truncation distance:
  truncation_length = -1;  // anything <1 and it won't truncate.

  // calculate the size of an f-bin:
  // note that you have to set a non-zero value to start or perp won't update.
  step.SetXYZ(1, 0, 0);
  step.SetPerp(dim.Perp() / nr);
  step.SetPhi(phispan / nphi);
  step.SetZ(dim.Z() / nz);
  std::cout << boost::str(boost::format("f-bin size:  r=%f,phi=%f, z=%f, wanted %f,%f") % step.Perp() % step.Phi() % ((rmax - rmin) / nr) % (phispan / nphi) % ((zmax - zmin) / nz)) << std::endl;

  // create an array to store the charge in each f-bin
  //  q = new MultiArray<double>(nr, nphi, nz);
  // q->SetAll(0);
  q = new ChargeMapReader(nr, rmin, rmax, nphi, 0, phispan, nz, zmin, zmax);
  chargestring = "No spacecharge present.";

  // load parameters of our region of interest
  rmin_roi = roi_r0;
  phimin_roi = roi_phi0;
  zmin_roi = roi_z0;  // lower edge of our region of interest, measured in f-bins
  rmax_roi = roi_r1;
  phimax_roi = roi_phi1;
  zmax_roi = roi_z1;  // exlcuded upper edge of our region of interest, measured in f-bins
  std::cout << boost::str(boost::format("AnnularFieldSim::AnnularFieldSim set roi variables rmin=%d phimin=%d zmin=%d rmax=%d phimax=%d zmax=%d") % rmin_roi % phimin_roi % zmin_roi % rmax_roi % phimax_roi % zmax_roi) << std::endl;
  // calculate the dimensions, in f-bins in our region of interest
  nr_roi = rmax_roi - rmin_roi;
  nphi_roi = phimax_roi - phimin_roi;
  nz_roi = zmax_roi - zmin_roi;
  std::cout << boost::str(boost::format("AnnularFieldSim::AnnularFieldSim calc'd roi variables nr=%d nphi=%d nz=%d") % nr_roi % nphi_roi % nz_roi) << std::endl;

  // create an array to hold the SC-induced electric field in the roi with the specified dimensions
  Efield = new MultiArray<TVector3>(nr_roi, nphi_roi, nz_roi);
  for (int i = 0; i < Efield->Length(); i++)
  {
    Efield->GetFlat(i)->SetXYZ(0, 0, 0);
  }

  // and to hold the external electric fieldmap over the region of interest
  Eexternal = new MultiArray<TVector3>(nr_roi, nphi_roi, nz_roi);
  for (int i = 0; i < Eexternal->Length(); i++)
  {
    Eexternal->GetFlat(i)->SetXYZ(0, 0, 0);
  }

  // ditto the external magnetic fieldmap
  Bfield = new MultiArray<TVector3>(nr_roi, nphi_roi, nz_roi);
  for (int i = 0; i < Bfield->Length(); i++)
  {
    Bfield->GetFlat(i)->SetXYZ(0, 0, 0);
  }

  // handle the lookup table construction:
  lookupCase = in_lookupCase;
  chargeCase = in_chargeCase;
  if (chargeCase == ChargeCase::NoSpacecharge)
  {
    lookupCase = LookupCase::NoLookup;  // don't build a lookup model if there's no charge.  It just wastes time.
  }

  if (lookupCase == Full3D)
  {
    std::cout << boost::str(boost::format("AnnularFieldSim::AnnularFieldSim building Epartial (full3D) with  nr_roi=%d nphi_roi=%d nz_roi=%d  =~%2.2fM TVector3 objects") % nr_roi % nphi_roi % nz_roi % (nr_roi * nphi_roi * nz_roi * nr * nphi * nz / (1.0e6))) << std::endl;

    Epartial = new MultiArray<TVector3>(nr_roi, nphi_roi, nz_roi, nr, nphi, nz);
    for (int i = 0; i < Epartial->Length(); i++)
    {
      Epartial->GetFlat(i)->SetXYZ(0, 0, 0);
    }
    // and kill the arrays we shouldn't be using:
    Epartial_highres = new MultiArray<TVector3>(1);
    Epartial_highres->GetFlat(0)->SetXYZ(0, 0, 0);

    Epartial_lowres = new MultiArray<TVector3>(1);
    Epartial_lowres->GetFlat(0)->SetXYZ(0, 0, 0);

    Epartial_phislice = new MultiArray<TVector3>(1);
    Epartial_phislice->GetFlat(0)->SetXYZ(0, 0, 0);
    q_lowres = new MultiArray<double>(1);
    *(q_lowres->GetFlat(0)) = 0;
    q_local = new MultiArray<double>(1);
    *(q_local->GetFlat(0)) = 0;
  }
  else if (lookupCase == HybridRes)
  {
    std::cout << "lookupCase==HybridRes" << std::endl;
    // zero out the other two:
    Epartial = new MultiArray<TVector3>(1);
    Epartial->GetFlat(0)->SetXYZ(0, 0, 0);

    Epartial_phislice = new MultiArray<TVector3>(1);
    Epartial_phislice->GetFlat(0)->SetXYZ(0, 0, 0);
  }
  else if (lookupCase == PhiSlice)
  {
    std::cout << "lookupCase==PhiSlice" << std::endl;

    Epartial_phislice = new MultiArray<TVector3>(nr_roi, 1, nz_roi, nr, nphi, nz);
    for (int i = 0; i < Epartial_phislice->Length(); i++)
    {
      Epartial_phislice->GetFlat(i)->SetXYZ(0, 0, 0);
    }
    // zero out the other two:
    Epartial = new MultiArray<TVector3>(1);
    Epartial->GetFlat(0)->SetXYZ(0, 0, 0);
    Epartial_highres = new MultiArray<TVector3>(1);
    Epartial_highres->GetFlat(0)->SetXYZ(0, 0, 0);

    Epartial_lowres = new MultiArray<TVector3>(1);
    Epartial_lowres->GetFlat(0)->SetXYZ(0, 0, 0);

    q_lowres = new MultiArray<double>(1);
    *(q_lowres->GetFlat(0)) = 0;
    q_local = new MultiArray<double>(1);
    *(q_local->GetFlat(0)) = 0;
  }
  else if (lookupCase == Analytic || lookupCase == NoLookup)
  {
    std::cout << "lookupCase==Analytic (or NoLookup)" << std::endl;

    // zero them all out:
    Epartial_phislice = new MultiArray<TVector3>(1);
    Epartial_phislice->GetFlat(0)->SetXYZ(0, 0, 0);

    Epartial = new MultiArray<TVector3>(1);
    Epartial->GetFlat(0)->SetXYZ(0, 0, 0);

    Epartial_highres = new MultiArray<TVector3>(1);
    Epartial_highres->GetFlat(0)->SetXYZ(0, 0, 0);

    Epartial_lowres = new MultiArray<TVector3>(1);
    Epartial_lowres->GetFlat(0)->SetXYZ(0, 0, 0);

    q_lowres = new MultiArray<double>(1);
    *(q_lowres->GetFlat(0)) = 0;
    q_local = new MultiArray<double>(1);
    *(q_local->GetFlat(0)) = 0;
  }
  else
  {
    std::cout << "Ran into wrong lookupCase logic in constructor." << std::endl;
    exit(1);
  }

  return;
}
AnnularFieldSim::AnnularFieldSim(float in_innerRadius, float in_outerRadius, float in_outerZ,
                                 int r, int roi_r0, int roi_r1, int in_rLowSpacing, int in_rHighSize,
                                 int phi, int roi_phi0, int roi_phi1, int in_phiLowSpacing, int in_phiHighSize,
                                 int z, int roi_z0, int roi_z1, int in_zLowSpacing, int in_zHighSize,
                                 float vdr, LookupCase in_lookupCase)
  : AnnularFieldSim(in_innerRadius, in_outerRadius, in_outerZ,
                    r, roi_r0, roi_r1, in_rLowSpacing, in_rHighSize,
                    phi, roi_phi0, roi_phi1, in_phiLowSpacing, in_phiHighSize,
                    z, roi_z0, roi_z1, in_zLowSpacing, in_zHighSize,
                    vdr, in_lookupCase, ChargeCase::FromFile)
{
  std::cout << "AnnularFieldSim::OldConstructor building AnnularFieldSim with ChargeCase::FromFile" << std::endl;
  return;
}
AnnularFieldSim::AnnularFieldSim(float in_innerRadius, float in_outerRadius, float in_outerZ,
                                 int r, int roi_r0, int roi_r1,
                                 int phi, int roi_phi0, int roi_phi1,
                                 int z, int roi_z0, int roi_z1,
                                 float vdr, LookupCase in_lookupCase)
  : AnnularFieldSim(in_innerRadius, in_outerRadius, in_outerZ,
                    r, roi_r0, roi_r1, 1, 3,
                    phi, roi_phi0, roi_phi1, 1, 3,
                    z, roi_z0, roi_z1, 1, 3,
                    vdr, in_lookupCase, ChargeCase::FromFile)
{
  // just passing through for the old version again
  // creates a region with high-res size of 3 (enough to definitely cover the eight local l-bins) and low-res spacing of 1, which ought to match the behavior (with a little more overhead) from when there was no highres-lowres distinction
  std::cout << "AnnularFieldSim::OldConstructor building AnnularFieldSim with local_size=1 in all dimensions, lowres_spacing=1 in all dimensions" << std::endl;
  return;
}
AnnularFieldSim::AnnularFieldSim(float rin, float rout, float dz, int r, int phi, int z, float vdr)
  : AnnularFieldSim(rin, rout, dz,
                    r, 0, r,
                    phi, 0, phi,
                    z, 0, z,
                    vdr, LookupCase::PhiSlice)
{
  // just a pass-through to go from the old style to the more detailed version.
  std::cout << "AnnularFieldSim::OldConstructor building AnnularFieldSim with roi=full in all dimensions" << std::endl;
  return;
}

TVector3 AnnularFieldSim::calc_unit_field(TVector3 at, TVector3 from)
{
  // if(debugFlag()) print_need_cout("%d: AnnularFieldSim::calc_unit_field(at=(r=%f,phi=%f,z=%f))\n",__LINE__,at.Perp(),at.Phi(),at.Z());
  int r_position;
  if (GetRindexAndCheckBounds(at.Perp(), &r_position) != InBounds)
  {
    std::cout << boost::str(boost::format("something's asking for 'at' with r=%f, which is index=%d") % at.Perp() % r_position) << std::endl;
    exit(1);
  }
  // note this is the field due to a fixed point charge in free space.
  // if doing cylindrical calcs with different boundary conditions, this needs to change.

  // this could check roi bounds before returning, if things start acting funny.

  TVector3 field(0, 0, 0);

  // if we're more than truncation_length away, don't bother? rcc food for thought.  right now _length is in bins, so tricky.

  // if the green's function class is not present use free space:
  if (green == nullptr)
  {
    TVector3 delr = at - from;
    field = delr;  // to set the direction.
    if (delr.Mag() < ALMOST_ZERO * ALMOST_ZERO)
    {  // note that this has blurred units -- it should scale with all three dimensions of stepsize.  For lots of phi bins, especially, this might start to read as small before it's really small.
       // do nothing.  the vector is already zero, which will be our approximation.
       // field.SetMag(0);//no contribution if we're in the same cell. -- but root warns if trying to resize something of magnitude zero.
    }
    else
    {
      field.SetMag(k_perm * 1 / (delr * delr));  // scalar product on the bottom. unitful, since we defined k_perm and delr with their correct units. (native units V=1 C=1 cm=1)
    }
    // print_need_cout("calc_unit_field at (%2.2f,%2.2f,%2.2f) from  (%2.2f,%2.2f,%2.2f).  Mag=%2.4fe-9\n",at.x(),at.Y(),at.Z(),from.X(),from.Y(),from.Z(),field.Mag()*1e9);
  }
  else
  {
    double atphi = FilterPhiPos(at.Phi());
    double fromphi = FilterPhiPos(from.Phi());
    double delphi = abs(at.DeltaPhi(from));
    // to allow us to use the same greens function set for both sides of the tpc, we shift into the valid greens region if needed:
    at.SetZ(at.Z() + green_shift);
    from.SetZ(from.Z() + green_shift);
    double Er = green->Er(at.Perp(), atphi, at.Z(), from.Perp(), fromphi, from.Z());
    // RCC manually disabled phi component of green -- actually, a correction to disallow trying to compute phi terms when at the same phi:
    double Ephi = 0;
    if (delphi > ALMOST_ZERO)
    {
      Ephi = green->Ephi(at.Perp(), atphi, at.Z(), from.Perp(), fromphi, from.Z());
    }
    double Ez = green->Ez(at.Perp(), atphi, at.Z(), from.Perp(), fromphi, from.Z());
    field.SetXYZ(-Er, -Ephi, -Ez);  // now these are the correct components if our test point is at y=0 (hence phi=0);
    field = field * epsinv;         // scale field strength, since the greens functions as of Apr 1 2020 do not build-in this factor.
    field.RotateZ(at.Phi());        // rotate to the coordinates of our 'at' point, which is a small rotation for the phislice case.
    // but this does mean we need to be careful! If we are rotating so that this is pointing not in the R direction, then with azimuthally symmetric charge we will have some cancellation we don't want, maybe?  Make sure this matches how we rotate to sum_field!
  }
  return field;
}

double AnnularFieldSim::FilterPhiPos(double phi)
{
  // this primarily takes the region [-pi,0] and maps it to [pi,2pi] by adding 2pi to it.
  // if math has pushed us past 2pi, it also subtracts to try to get us in range.
  double p = phi;
  if (p >= 2 * M_PI)  // phispan)
  {
    p -= 2 * M_PI;
  }
  if (p < 0)
  {
    p += 2 * M_PI;
  }
  if (p >= 2 * M_PI || p < 0)
  {
    std::cout << boost::str(boost::format("AnnularFieldSim::FilterPhiPos asked to filter %f, which is more than range=%f out of bounds.  Check what called this.") % phi % (2 * M_PI)) << std::endl;
    exit(1);
  }
  return p;
}
int AnnularFieldSim::FilterPhiIndex(int phi, int range = -1)
{
  if (range < 0)
  {
    range = nphi;  // default input is range=-1.  in that case, use the intrinsic resolution of the q grid.
  }

  // shifts phi up or down by nphi (=2pi in phi index space), and complains if this doesn't put it in range.
  int p = phi;
  if (p >= range)
  {
    p -= range;
  }
  if (p < 0)
  {
    p += range;
  }
  if (p >= range || p < 0)
  {
    std::cout << boost::str(boost::format("AnnularFieldSim::FilterPhiIndex asked to filter %d, which is more than range=%d out of bounds.  Check what called this.") % phi % range) << std::endl;
    exit(1);
  }
  return p;
}

int AnnularFieldSim::GetRindex(float pos)
{
  float r0f = (pos - rmin) / step.Perp();  // the position in r, in units of step, starting from the low edge of the 0th bin.
  int r0 = floor(r0f);
  return r0;
}

int AnnularFieldSim::GetPhiIndex(float pos)
{
  float p0f = (pos) / step.Phi();  // the position in phi, in units of step, starting from the low edge of the 0th bin.
  int phitemp = floor(p0f);
  int p0 = FilterPhiIndex(phitemp);
  return p0;
}

int AnnularFieldSim::GetZindex(float pos)
{
  float z0f = (pos - zmin) / step.Z();  // the position in z, in units of step, starting from the low edge of the 0th bin.
  int z0 = floor(z0f);
  return z0;
}

AnnularFieldSim::BoundsCase AnnularFieldSim::GetRindexAndCheckBounds(float pos, int *r)
{
  // if(debugFlag()) print_need_cout("%d: AnnularFieldSim::GetRindexAndCheckBounds(r=%f)\n",__LINE__,pos);

  float r0f = (pos - rmin) / step.Perp();  // the position in r, in units of step, starting from the low edge of the 0th bin.
  int r0 = floor(r0f);
  *r = r0;

  int r0lowered_slightly = floor(r0f - ALMOST_ZERO);
  int r0raised_slightly = floor(r0f + ALMOST_ZERO);
  if (r0lowered_slightly >= rmax_roi || r0raised_slightly < rmin_roi)
  {
    return OutOfBounds;
  }

  // now if we are out of bounds, it is because we are on an edge, within range of ALMOST_ZERO of being in fair territory.
  if (r0 >= rmax_roi)
  {
    return OnHighEdge;
  }
  if (r0 < rmin_roi)
  {
    return OnLowEdge;
  }
  // if we're still here, we're in bounds.
  return InBounds;
}
AnnularFieldSim::BoundsCase AnnularFieldSim::GetPhiIndexAndCheckBounds(float pos, int *phi)
{
  // if(debugFlag()) print_need_cout("%d: AnnularFieldSim::GetPhiIndexAndCheckBounds(phi=%f)\n\n",__LINE__,pos);
  float p0f = (pos) / step.Phi();  // the position in phi, in units of step, starting from the low edge of the 0th bin.
  int phitemp = floor(p0f);
  int p0 = FilterPhiIndex(phitemp);
  *phi = p0;

  phitemp = floor(p0f - ALMOST_ZERO);
  int p0lowered_slightly = FilterPhiIndex(phitemp);
  phitemp = floor(p0f + ALMOST_ZERO);
  int p0raised_slightly = FilterPhiIndex(phitemp);
  // annoying detail:  if we are at index 0, we might go above pmax by going down.
  //  and if we are at nphi-1, we might go below pmin by going up.
  //  if we are not at p0=0 or nphi-1, the slight wiggles won't move us.
  //  if we are at p0=0, we are not at or above the max, and lowering slightly won't change that,
  //  and is we are at p0=nphi-1, we are not below the min, and raising slightly won't change that
  if ((p0lowered_slightly >= phimax_roi && p0 != 0) || (p0raised_slightly < phimin_roi && p0 != (nphi - 1)))
  {
    return OutOfBounds;
  }
  // now if we are out of bounds, it is because we are on an edge, within range of ALMOST_ZERO of being in fair territory.
  if (p0 >= phimax_roi)
  {
    return OnHighEdge;
  }
  if (p0 < phimin_roi)
  {
    return OnLowEdge;
  }
  // if we're still here, we're in bounds.
  return InBounds;
}

AnnularFieldSim::BoundsCase AnnularFieldSim::GetZindexAndCheckBounds(float pos, int *z)
{
  // if(debugFlag()) print_need_cout("%d: AnnularFieldSim::GetZindexAndCheckBounds(z=%f)\n\n",__LINE__,pos);
  float z0f = (pos - zmin) / step.Z();  // the position in z, in units of step, starting from the low edge of the 0th bin.
  int z0 = floor(z0f);
  *z = z0;

  int z0lowered_slightly = floor(z0f - ALMOST_ZERO);
  int z0raised_slightly = floor(z0f + ALMOST_ZERO);

  if (z0lowered_slightly >= zmax_roi || z0raised_slightly < zmin_roi)
  {
    return OutOfBounds;
  }
  // now if we are out of bounds, it is because we are on an edge, within range of ALMOST_ZERO of being in fair territory.
  if (z0 >= zmax_roi)
  {
    return OnHighEdge;
  }
  if (z0 < zmin_roi)
  {
    return OnLowEdge;
  }
  // if we're still here, we're in bounds.
  return InBounds;
}

TVector3 AnnularFieldSim::analyticFieldIntegral(float zdest, TVector3 start, MultiArray<TVector3> *field)
{
  // integrates E dz, from the starting point to the selected z position.  The path is assumed to be along z for each step, with adjustments to x and y accumulated after each step.
  // if(debugFlag()) print_need_cout("%d: AnnularFieldSim::fieldIntegral(x=%f,y=%f, z=%f) to z=%f\n\n",__LINE__,start.X(),start.Y(),start.Z(),zdest);
  //   print_need_cout("AnnularFieldSim::analyticFieldIntegral calculating from (%f,%f,%f) (rphiz)=(%f,%f,%f) to z=%f.\n",start.X(),start.Y(),start.Z(),start.Perp(),start.Phi(),start.Z(),zdest);

  // coordinates are assumed to be in native units (cm=1);

  int r, phi;
  bool rOkay = (GetRindexAndCheckBounds(start.Perp(), &r) == InBounds);
  bool phiOkay = (GetPhiIndexAndCheckBounds(FilterPhiPos(start.Phi()), &phi) == InBounds);

  // bool isE=(field==Efield);
  // bool isB=(field==Bfield);

  // print_need_cout("anaFieldInt:  isE=%d, isB=%d, start=(%f,%f,%f)\n",isE,isB,start.X(),start.Y(),start.Z());

  if (!rOkay || !phiOkay)
  {
    std::cout << boost::str(boost::format("AnnularFieldSim::analyticFieldIntegral asked to integrate along (r=%f,phi=%f), index=(%d,%d), which is out of bounds.  returning starting position.") % start.Perp() % start.Phi() % r % phi) << std::endl;
    return start;
  }

  int dir = (start.Z() > zdest ? -1 : 1);  //+1 if going to larger z, -1 if going to smaller;  if they're the same, the sense doesn't matter.

  int zi, zf;
  double startz, endz;
  BoundsCase startBound, endBound;

  // make sure 'zi' is always the smaller of the two numbers, for handling the partial-steps.
  if (dir > 0)
  {
    startBound = GetZindexAndCheckBounds(start.Z(), &zi);  // highest cell with lower bound less than lower bound of integral
    endBound = GetZindexAndCheckBounds(zdest, &zf);        // highest cell with lower bound less than upper lower bound of integral
    startz = start.Z();
    endz = zdest;
  }
  else
  {
    endBound = GetZindexAndCheckBounds(start.Z(), &zf);  // highest cell with lower bound less than lower bound of integral
    startBound = GetZindexAndCheckBounds(zdest, &zi);    // highest cell with lower bound less than upper lower bound of integral
    startz = zdest;
    endz = start.Z();
  }
  bool startOkay = (startBound == InBounds);
  bool endOkay = (endBound == InBounds || endBound == OnHighEdge);  // if we just barely touch out-of-bounds on the high end, we can skip that piece of the integral

  if (!startOkay || !endOkay)
  {
    std::cout << boost::str(boost::format("AnnularFieldSim::analyticFieldIntegral asked to integrate from z=%f to %f, index=%d to %d), which is out of bounds.  returning starting position.") % startz % endz % zi % zf) << std::endl;
    return start;
  }

  start.SetZ(startz);
  TVector3 integral;
  if (field == Efield)
  {
    integral = aliceModel->Eint(endz, start) + Eexternal->Get(r - rmin_roi, phi - phimin_roi, zi - zmin_roi) * (endz - startz);
    return dir * integral;
  }
  else if (field == Bfield)
  {
    return interpolatedFieldIntegral(zdest, start, Bfield);
  }
  return integral;
}

TVector3 AnnularFieldSim::fieldIntegral(float zdest, const TVector3 &start, MultiArray<TVector3> *field)
{
  // integrates E dz, from the starting point to the selected z position.  The path is assumed to be along z for each step, with adjustments to x and y accumulated after each step.
  // if(debugFlag()) print_need_cout("%d: AnnularFieldSim::fieldIntegral(x=%f,y=%f, z=%f) to z=%f\n\n",__LINE__,start.X(),start.Y(),start.Z(),zdest);

  int r, phi;
  bool rOkay = (GetRindexAndCheckBounds(start.Perp(), &r) == InBounds);
  bool phiOkay = (GetPhiIndexAndCheckBounds(FilterPhiPos(start.Phi()), &phi) == InBounds);

  if (!rOkay || !phiOkay)
  {
    std::cout << boost::str(boost::format("AnnularFieldSim::fieldIntegral asked to integrate along (r=%f,phi=%f), index=(%d,%d), which is out of bounds.  returning starting position.") % start.Perp() % start.Phi() % r % phi) << std::endl;
    return start;
  }

  int dir = (start.Z() > zdest ? -1 : 1);  //+1 if going to larger z, -1 if going to smaller;  if they're the same, the sense doesn't matter.

  int zi, zf;
  double startz, endz;
  BoundsCase startBound, endBound;

  // make sure 'zi' is always the smaller of the two numbers, for handling the partial-steps.
  if (dir > 0)
  {
    startBound = GetZindexAndCheckBounds(start.Z(), &zi);  // highest cell with lower bound less than lower bound of integral
    endBound = GetZindexAndCheckBounds(zdest, &zf);        // highest cell with lower bound less than upper lower bound of integral
    startz = start.Z();
    endz = zdest;
  }
  else
  {
    endBound = GetZindexAndCheckBounds(start.Z(), &zf);  // highest cell with lower bound less than lower bound of integral
    startBound = GetZindexAndCheckBounds(zdest, &zi);    // highest cell with lower bound less than upper lower bound of integral
    startz = zdest;
    endz = start.Z();
  }
  bool startOkay = (startBound == InBounds);
  bool endOkay = (endBound == InBounds || endBound == OnHighEdge);  // if we just barely touch out-of-bounds on the high end, we can skip that piece of the integral

  if (!startOkay || !endOkay)
  {
    std::cout << boost::str(boost::format("AnnularFieldSim::fieldIntegral asked to integrate from z=%f to %f, index=%d to %d), which is out of bounds.  returning starting position.") % startz % endz % zi % zf) << std::endl;
    return start;
  }

  TVector3 fieldInt(0, 0, 0);
  // print_need_cout("AnnularFieldSim::fieldIntegral requesting (%d,%d,%d)-(%d,%d,%d) (inclusive) cells\n",r,phi,zi,r,phi,zf-1);
  TVector3 tf;
  for (int i = zi; i < zf; i++)
  {  // count the whole cell of the lower end, and skip the whole cell of the high end.
    tf = field->Get(r - rmin_roi, phi - phimin_roi, i - zmin_roi);
    // print_need_cout("fieldAt (%d,%d,%d)=(%f,%f,%f) step=%f\n",r,phi,i,tf.X(),tf.Y(),tf.Z(),step.Z());
    fieldInt += tf * step.Z();
  }

  // since bins contain their lower bound, but not their upper, I can safely remove the unused portion of the lower cell:
  fieldInt -= field->Get(r - rmin_roi, phi - phimin_roi, zi - zmin_roi) * (startz - zi * step.Z());  // remove the part of the low end cell we didn't travel through

  // but only need to add the used portion of the upper cell if we go past the edge of it meaningfully:
  if (endz / step.Z() - zf > ALMOST_ZERO)
  {
    // print_need_cout("endz/step.Z()=%f, zf=%f\n",endz/step.Z(),zf*1.0);
    // if our final step is actually in the next step.
    fieldInt += field->Get(r - rmin_roi, phi - phimin_roi, zf - zmin_roi) * (endz - zf * step.Z());  // add the part of the high end cell we did travel through
  }

  return dir * fieldInt;
}

TVector3 AnnularFieldSim::GetCellCenter(int r, int phi, int z)
{
  // returns the midpoint of the cell (halfway between each edge, not weighted center)

  TVector3 c(1, 0, 0);
  c.SetPerp((r + 0.5) * step.Perp() + rmin);
  c.SetPhi((phi + 0.5) * step.Phi());
  c.SetZ((z + 0.5) * step.Z() + zmin);

  return c;
}

TVector3 AnnularFieldSim::GetRoiCellCenter(int r, int phi, int z)
{
  // returns the midpoint of the cell relative to the roi (halfway between each edge, not weighted center)

  // return zero if it's out of bounds:
  if (r > nr_roi || r < 0 || phi > nphi_roi || phi < 0 || z > nz_roi || z < 0)
  {
    return zero_vector;
  }

  TVector3 c(1, 0, 0);
  c.SetPerp((r + rmin_roi + 0.5) * step.Perp() + rmin);
  c.SetPhi((phi + phimin_roi + 0.5) * step.Phi());
  c.SetZ((z + zmin_roi + 0.5) * step.Z() + zmin);

  return c;
}

TVector3 AnnularFieldSim::GetGroupCellCenter(int r0, int r1, int phi0, int phi1, int z0, int z1)
{
  // returns the midpoint of the cell (halfway between each edge, not weighted center)
  float ravg = (r0 + r1) / 2.0 + 0.5;
  if (phi0 > phi1)
  {
    phi1 += nphi;
  }
  if (phi0 > phi1)
  {
    std::cout << boost::str(boost::format("phi1(%d)<=phi0(%d) even after boosting phi1.  check what called this!") % phi1 % phi0) << std::endl;
    exit(1);
  }
  float phiavg = (r0 + r1) / 2.0 + 0.5;
  if (phiavg >= nphi)
  {
    phiavg -= nphi;
  }

  float zavg = (z0 + z1) / 2.0 + 0.5;

  TVector3 c(1, 0, 0);
  c.SetPerp((ravg) *step.Perp() + rmin);
  c.SetPhi((phiavg) *step.Phi());
  c.SetZ((zavg) *step.Z() + zmin);

  return c;
}

TVector3 AnnularFieldSim::GetWeightedCellCenter(int r, int phi, int z)
{
  // returns the weighted center of the cell by volume.
  // todo:  this is vaguely hefty, and might be worth storing the result of, if speed is needed
  TVector3 c(1, 0, 0);

  float rin = r * step.Perp() + rmin;
  float rout = rin + step.Perp();

  float rMid = (4 * sin(step.Phi() / 2) * (pow(rout, 3) - pow(rin, 3)) / (3 * step.Phi() * (pow(rout, 2) - pow(rin, 2))));
  c.SetPerp(rMid);
  c.SetPhi((phi + 0.5) * step.Phi());
  c.SetZ((z + 0.5) * step.Z());

  return c;
}

TVector3 AnnularFieldSim::interpolatedFieldIntegral(float zdest, const TVector3 &start, MultiArray<TVector3> *field)
{
  // print_need_cout("AnnularFieldSim::interpolatedFieldIntegral(x=%f,y=%f, z=%f)\n",start.X(),start.Y(),start.Z());

  float r0 = (start.Perp() - rmin) / step.Perp() - 0.5;  // the position in r, in units of step, starting from the center of the 0th bin.
  int r0i = floor(r0);                                   // the integer portion of the position. -- what center is below our position?
  float r0d = r0 - r0i;                                  // the decimal portion of the position. -- how far past center are we?
  int ri[4];                                             // the r position of the four cell centers to consider
  ri[0] = ri[1] = r0i;
  ri[2] = ri[3] = r0i + 1;
  float rw[4];              // the weight of that cell, linear in distance from the center of it
  rw[0] = rw[1] = 1 - r0d;  // 1 when we're on it, 0 when we're at the other one.
  rw[2] = rw[3] = r0d;      // 1 when we're on it, 0 when we're at the other one

  bool skip[] = {false, false, false, false};
  if (ri[0] < rmin_roi)
  {
    skip[0] = true;  // don't bother handling 0 and 1 in the coordinates.
    skip[1] = true;
    rw[2] = rw[3] = 1;  // and weight like we're dead-center on the outer cells.
  }
  if (ri[2] >= rmax_roi)
  {
    skip[2] = true;  // don't bother handling 2 and 3 in the coordinates.
    skip[3] = true;
    rw[0] = rw[1] = 1;  // and weight like we're dead-center on the inner cells.
  }

  // now repeat that structure for phi:
  float p0 = (start.Phi()) / step.Phi() - 0.5;  // the position in phi, in units of step, starting from the center of the 0th bin.
  int p0i = floor(p0);                          // the integer portion of the position. -- what center is below our position?
  float p0d = p0 - p0i;                         // the decimal portion of the position. -- how far past center are we?
  int pi[4];                                    // the phi position of the four cell centers to consider
  pi[0] = pi[2] = FilterPhiIndex(p0i);
  pi[1] = pi[3] = FilterPhiIndex(p0i + 1);
  float pw[4];              // the weight of that cell, linear in distance from the center of it
  pw[0] = pw[2] = 1 - p0d;  // 1 when we're on it, 0 when we're at the other one.
  pw[1] = pw[3] = p0d;      // 1 when we're on it, 0 when we're at the other one

  // to handle wrap-around
  if (pi[0] < phimin_roi || pi[0] >= phimax_roi)
  {
    skip[0] = true;  // don't bother handling 0 and 1 in the coordinates.
    skip[2] = true;
    pw[1] = pw[3] = 1;  // and weight like we're dead-center on the outer cells.
  }
  if (pi[1] >= phimax_roi || pi[1] < phimin_roi)
  {
    skip[1] = true;  // don't bother handling 2 and 3 in the coordinates.
    skip[3] = true;
    pw[0] = pw[2] = 1;  // and weight like we're dead-center on the inner cells.
  }

  int dir = (start.Z() < zdest ? 1 : -1);  //+1 if going to larger z, -1 if going to smaller;  if they're the same, the sense doesn't matter.

  int zi, zf;
  double startz, endz;
  BoundsCase startBound, endBound;

  // make sure 'zi' is always the smaller of the two numbers, for handling the partial-steps.
  if (dir > 0)
  {
    startBound = GetZindexAndCheckBounds(start.Z(), &zi);  // highest cell with lower bound less than lower bound of integral
    endBound = GetZindexAndCheckBounds(zdest, &zf);        // highest cell with lower bound less than upper bound of integral
    startz = start.Z();
    endz = zdest;
  }
  else
  {
    endBound = GetZindexAndCheckBounds(start.Z(), &zf);  // highest cell with lower bound less than lower bound of integral
    startBound = GetZindexAndCheckBounds(zdest, &zi);    // highest cell with lower bound less than upper bound of integral
    startz = zdest;
    endz = start.Z();
  }
  bool startOkay = (startBound == InBounds || startBound == OnLowEdge);  // maybe todo: add handling for being just below the low edge.
  bool endOkay = (endBound == InBounds || endBound == OnHighEdge);       // if we just barely touch out-of-bounds on the high end, we can skip that piece of the integral

  if (!startOkay || !endOkay)
  {
    std::cout << boost::str(boost::format("AnnularFieldSim::InterpolatedFieldIntegral asked to integrate from z=%f to %f, index=%d to %d), which is out of bounds.  returning zero") % startz % endz % zi % zf) << std::endl;
    return zero_vector;
  }

  if (startBound == OnLowEdge)
  {
    // we were just below the low edge, so we will be asked to sample a bin in z we're not actually using
    zi++;  // avoid it.  We weren't integrating anything in it anyway.
  }

  if (false)
  {  // debugging options
    bool isE = (field == Efield);
    bool isB = (field == Bfield);
    std::cout << boost::str(boost::format("interpFieldInt:  isE=%d, isB=%d, start=(%2.4f,%2.4f,%2.4f) dest=%2.4f") % isE % isB % start.X() % start.Y() % start.Z() % zdest) << std::endl;
    std::cout << boost::str(boost::format("interpolating fieldInt for %s at  r=%f,phi=%f") % ((field == Efield) ? "Efield" : "Bfield") % r0 % p0) << std::endl;
    std::cout << boost::str(boost::format("direction=%d, startz=%2.4f, endz=%2.4f, zi=%d, zf=%d") % dir % startz % endz % zi % zf) << std::endl;
    for (int i = 0; i < 4; i++)
    {
      std::cout << boost::str(boost::format("skip[%d]=%d\tw[i]=%2.2f, (pw=%2.2f, rw=%2.2f)") % i % skip[i] % (rw[i] * pw[i]) % pw[i] % rw[i]) << std::endl;
    }
  }

  TVector3 fieldInt(0, 0, 0), partialInt;  // where we'll store integrals as we generate them.

  for (int i = 0; i < 4; i++)
  {
    if (skip[i])
    {
      // print_need_cout("skipping element r=%d,phi=%d\n",ri[i],pi[i]);
      continue;  // we invalidated this one for some reason.
    }
    partialInt.SetXYZ(0, 0, 0);
    for (int j = zi; j < zf; j++)
    {  // count the whole cell of the lower end, and skip the whole cell of the high end.

      partialInt += field->Get(ri[i] - rmin_roi, pi[i] - phimin_roi, j - zmin_roi) * step.Z();
    }
    if (startBound != OnLowEdge)
    {
      partialInt -= field->Get(ri[i] - rmin_roi, pi[i] - phimin_roi, zi - zmin_roi) * (startz - (zi * step.Z() + zmin));  // remove the part of the low end cell we didn't travel through
      // print_need_cout("removing low end of cell we didn't travel through (zi-zmin_roi=%d, length=%f)\n",zi-zmin_roi, startz-(zi*step.Z()+zmin));
    }
    if ((endz - zmin) / step.Z() - zf > ALMOST_ZERO)
    {
      partialInt += field->Get(ri[i] - rmin_roi, pi[i] - phimin_roi, zf - zmin_roi) * (endz - (zf * step.Z() + zmin));  // add the part of the high end cell we did travel through
      // print_need_cout("adding low end of cell we did travel through (zf-zmin_roi=%d, length=%f)\n",zf-zmin_roi, endz-(zf*step.Z()+zmin));
    }
    // print_need_cout("element r=%d,phi=%d, w=%f partialInt=(%2.2E,%2.2E,%2.2E)\n",ri[i],pi[i],rw[i]*pw[i],partialInt.X(),partialInt.Y(),partialInt.Z());

    fieldInt += rw[i] * pw[i] * partialInt;
  }

  return dir * fieldInt;
}

void AnnularFieldSim::load_analytic_spacecharge(float scalefactor = 1)
{
  // scalefactor should be chosen so Rho(pos) returns C/cm^3

  // sphenix:
  // double ifc_radius=20;
  // double ofc_radius=78;
  //   double tpc_halfz=

  double ifc_radius = 83.5 * cm;
  double ofc_radius = 254.5 * cm;
  double tpc_halfz = 250 * cm;

  aliceModel = new AnalyticFieldModel(ifc_radius / cm, ofc_radius / cm, tpc_halfz / cm, scalefactor);
  double totalcharge = 0;
  double localcharge = 0;

  TVector3 pos;
  for (int ifr = 0; ifr < nr; ifr++)
  {
    for (int ifphi = 0; ifphi < nphi; ifphi++)
    {
      for (int ifz = 0; ifz < nz; ifz++)
      {
        pos = GetCellCenter(ifr, ifphi, ifz);
        pos = pos * (1.0 / cm);  // divide by cm so we're in units of cm when we query the charge model.
        double vol = step.Z() * step.Phi() * (2 * (ifr * step.Perp() + rmin) + step.Perp()) * step.Perp();
        // if(debugFlag()) print_need_cout("%d: AnnularFieldSim::load_analytic_spacecharge adding Q=%f into cell (%d,%d,%d)\n",__LINE__,qbin,i,j,k,localr,localphi,localz);
        localcharge = vol * aliceModel->Rho(pos);  // TODO:  figure out what units this is in.
        totalcharge += localcharge;
        // q->Add(ifr, ifphi, ifz, localcharge);  //scalefactor must be applied to charge _and_ field, and so is handled in the aliceModel code.
        q->AddChargeInBin(ifr, ifphi, ifz, localcharge);  // scalefactor must be applied to charge _and_ field, and so is handled in the aliceModel code.
      }
    }
  }
  std::cout << boost::str(boost::format("AnnularFieldSim::load_analytic_spacecharge:  Total charge Q=%E Coulombs") % totalcharge) << std::endl;

  if (lookupCase == HybridRes)
  {
    // go through the q array and build q_lowres.
    for (int i = 0; i < q_lowres->Length(); i++)
    {
      *(q_lowres->GetFlat(i)) = 0;
    }

    // fill our low-res
    // note that this assumes the last bin is short or normal length, not long.
    for (int ifr = 0; ifr < nr; ifr++)
    {
      int r_low = ifr / r_spacing;  // index of our l-bin is just the integer division of the index of our f-bin
      for (int ifphi = 0; ifphi < nphi; ifphi++)
      {
        int phi_low = ifphi / phi_spacing;
        for (int ifz = 0; ifz < nz; ifz++)
        {
          int z_low = ifz / z_spacing;
          q_lowres->Add(r_low, phi_low, z_low, q->GetChargeInBin(ifr, ifphi, ifz));
        }
      }
    }
  }
  return;
}

void AnnularFieldSim::loadEfield(const std::string &filename, const std::string &treename, int zsign)
{
  // prep variables so that loadField can just iterate over the tree entries and fill our selected tree agnostically
  // assumes file stores fields as V/cm.
  TFile fieldFile(filename.c_str(), "READ");
  TTree *fTree;
  fieldFile.GetObject(treename.c_str(), fTree);
  Efieldname = "E:" + filename + ":" + treename;
  float r, z;          // coordinates
  float fr, fz, fphi;  // field components at that coordinate
  fTree->SetBranchAddress("r", &r);
  fTree->SetBranchAddress("er", &fr);
  fTree->SetBranchAddress("z", &z);
  fTree->SetBranchAddress("ez", &fz);
  // phi would go here if we had it.
  fphi = 0;  // no phi components yet.
  // float phi = 0.;
  // phi += 1;
  // phi = 0;  //satisfy picky racf compiler
  loadField(&Eexternal, fTree, &r, nullptr, &z, &fr, &fphi, &fz, V / cm, zsign);
  fieldFile.Close();
  return;
}
void AnnularFieldSim::loadBfield(const std::string &filename, const std::string &treename)
{
  // prep variables so that loadField can just iterate over the tree entries and fill our selected tree agnostically
  // assumes file stores field as Tesla.
  TFile fieldFile(filename.c_str(), "READ");
  TTree *fTree;
  fieldFile.GetObject(treename.c_str(), fTree);
  Bfieldname = "B:" + filename + ":" + treename;
  float r, z;          // coordinates
  float fr, fz, fphi;  // field components at that coordinate
  fTree->SetBranchAddress("r", &r);
  fTree->SetBranchAddress("br", &fr);
  fTree->SetBranchAddress("z", &z);
  fTree->SetBranchAddress("bz", &fz);
  // phi would go here if we had it.
  fphi = 0;  // no phi components yet.
  // float phi = 0.;
  // phi += 1;
  // phi = 0;  //satisfy picky racf compiler
  loadField(&Bfield, fTree, &r, nullptr, &z, &fr, &fphi, &fz, Tesla, 1);
  fieldFile.Close();

  return;
}

void AnnularFieldSim::load3dBfield(const std::string &filename, const std::string &treename, int zsign, float scale)
{
  // prep variables so that loadField can just iterate over the tree entries and fill our selected tree agnostically
  // assumes file stores field as Tesla.
  TFile fieldFile(filename.c_str(), "READ");
  TTree *fTree;
  fieldFile.GetObject(treename.c_str(), fTree);
  Bfieldname = "B(3D):" + filename + ":" + treename + " Scale =" + boost::str(boost::format("%2.2f") % scale);
  float r, z, phi;     // coordinates
  float fr, fz, fphi;  // field components at that coordinate
  fTree->SetBranchAddress("rho", &r);
  fTree->SetBranchAddress("brho", &fr);
  fTree->SetBranchAddress("z", &z);
  fTree->SetBranchAddress("bz", &fz);
  fTree->SetBranchAddress("phi", &phi);
  fTree->SetBranchAddress("bphi", &fphi);
  loadField(&Bfield, fTree, &r, nullptr, &z, &fr, &fphi, &fz, Tesla * scale, zsign);
  return;
}

void AnnularFieldSim::loadField(MultiArray<TVector3> **field, TTree *source, float *rptr, float *phiptr, float *zptr, float *frptr, float *fphiptr, float *fzptr, float fieldunit, int zsign)
{
  // we're loading a tree of unknown size and spacing -- and possibly uneven spacing -- into our local data.
  // formally, we might want to interpolate or otherwise weight, but for now, carve this into our usual bins, and average, similar to the way we load spacecharge.

  bool phiSymmetry = (phiptr == nullptr);  // if the phi pointer is zero, assume phi symmetry.
  int lowres_factor = 10;                  // to fill in gaps, we group together loweres^3 cells into one block and use that average.
  std::cout << boost::str(boost::format("loading field from %f<z<%f") % zmin % zmax) << std::endl;
  TH3F *htEntries = new TH3F("htentries", "num of entries in the field loading", nphi, 0, M_PI * 2.0, nr, rmin, rmax, nz, zmin, zmax);
  TH3F *htSum[3];
  TH3F *htEntriesLow = new TH3F("htentrieslow", "num of lowres entries in the field loading", nphi / lowres_factor + 1, 0, M_PI * 2.0, nr / lowres_factor + 1, rmin, rmax, nz / lowres_factor + 1, zmin, zmax);
  TH3F *htSumLow[3];
  std::string axis[]{"r", "p", "z"};
  for (int i = 0; i < 3; i++)
  {
    htSum[i] = new TH3F(std::string("htsum" + std::to_string(i)).c_str(), std::string("sum of " + axis[i] + "-axis entries in the field loading").c_str(), nphi, 0, M_PI * 2.0, nr, rmin, rmax, nz, zmin, zmax);
    htSumLow[i] = new TH3F(std::string("htsumlow" + std::to_string(i)).c_str(), std::string("sum of low " + axis[i] + "-axis entries in the field loading").c_str(), nphi / lowres_factor + 1, 0, M_PI * 2.0, nr / lowres_factor + 1, rmin, rmax, nz / lowres_factor + 1, zmin, zmax);
  }

  int nEntries = source->GetEntries();
  for (int i = 0; i < nEntries; i++)
  {  // could probably do this with an iterator
    source->GetEntry(i);
    float zval = *zptr * zsign;  // right now, need the ability to flip the sign of the z coordinate.
    // note that the z sign also needs to affect the field sign in that direction, which is handled outside in the z components of the fills
    // if we aren't asking for phi symmetry, build just the one phi strip
    if (!phiSymmetry)
    {
      assert(phiptr);
      htEntries->Fill(*phiptr, *rptr, zval);  // for legacy reasons this histogram, like others, goes phi-r-z.
      htSum[0]->Fill(*phiptr, *rptr, zval, *frptr * fieldunit);
      htSum[1]->Fill(*phiptr, *rptr, zval, *fphiptr * fieldunit);
      htSum[2]->Fill(*phiptr, *rptr, zval, *fzptr * fieldunit * zsign);
      htEntriesLow->Fill(*phiptr, *rptr, zval);  // for legacy reasons this histogram, like others, goes phi-r-z.
      htSumLow[0]->Fill(*phiptr, *rptr, zval, *frptr * fieldunit);
      htSumLow[1]->Fill(*phiptr, *rptr, zval, *fphiptr * fieldunit);
      htSumLow[2]->Fill(*phiptr, *rptr, zval, *fzptr * fieldunit * zsign);
    }
    else
    {  // if we do have phi symmetry, build every phi strip using this one.
      for (int j = 0; j < nphi; j++)
      {
        htEntries->Fill(j * step.Phi(), *rptr, zval);  // for legacy reasons this histogram, like others, goes phi-r-z.
        htSum[0]->Fill(j * step.Phi(), *rptr, zval, *frptr * fieldunit);
        htSum[1]->Fill(j * step.Phi(), *rptr, zval, *fphiptr * fieldunit);
        htSum[2]->Fill(j * step.Phi(), *rptr, zval, *fzptr * fieldunit * zsign);
        htEntriesLow->Fill(j * step.Phi(), *rptr, zval);  // for legacy reasons this histogram, like others, goes phi-r-z.
        htSumLow[0]->Fill(j * step.Phi(), *rptr, zval, *frptr * fieldunit);
        htSumLow[1]->Fill(j * step.Phi(), *rptr, zval, *fphiptr * fieldunit);
        htSumLow[2]->Fill(j * step.Phi(), *rptr, zval, *fzptr * fieldunit * zsign);
      }
    }
  }
  // now we just divide and fill our local plots (which should eventually be stored as histograms, probably) with the values from each hist cell:
  int nemptybins = 0;
  for (int i = 0; i < nphi; i++)
  {
    for (int j = 0; j < nr; j++)
    {
      for (int k = 0; k < nz; k++)
      {
        TVector3 cellcenter = GetCellCenter(j, i, k);
        int bin = htEntries->FindBin(FilterPhiPos(cellcenter.Phi()), cellcenter.Perp(), cellcenter.Z());
        TVector3 fieldvec(htSum[0]->GetBinContent(bin), htSum[1]->GetBinContent(bin), htSum[2]->GetBinContent(bin));
        fieldvec = fieldvec * (1.0 / htEntries->GetBinContent(bin));
        if (htEntries->GetBinContent(bin) < 0.99)
        {
          // no entries here!
          nemptybins++;
        }
        // have to rotate this to the proper direction.
        fieldvec.RotateZ(FilterPhiPos(cellcenter.Phi()));  // rcc caution.  Does this rotation shift the sense of 'up'?
        (*field)->Set(j, i, k, fieldvec);
      }
    }
  }
  if (nemptybins > 0)
  {
    std::cout << boost::str(boost::format("found %d empty bins when constructing %s.  Filling with lower resolution.") % nemptybins % ((*field == Bfield) ? "Bfield" : "Eexternal")) << std::endl;
    for (int i = 0; i < nphi; i++)
    {
      for (int j = 0; j < nr; j++)
      {
        for (int k = 0; k < nz; k++)
        {
          TVector3 cellcenter = GetCellCenter(j, i, k);
          int bin = htEntries->FindBin(FilterPhiPos(cellcenter.Phi()), cellcenter.Perp(), cellcenter.Z());
          if (htEntries->GetBinContent(bin) == 0)
          {
            int lowbin = htEntriesLow->FindBin(FilterPhiPos(cellcenter.Phi()), cellcenter.Perp(), cellcenter.Z());
            TVector3 fieldvec(htSumLow[0]->GetBinContent(lowbin), htSumLow[1]->GetBinContent(lowbin), htSumLow[2]->GetBinContent(lowbin));
            fieldvec = fieldvec * (1.0 / htEntriesLow->GetBinContent(lowbin));
            if (htEntriesLow->GetBinContent(lowbin) < 0.99)
            {
              std::cout << boost::str(boost::format("not enough entries in source to fill fieldmap.  None near r=%f, phi=%f, z=%f. Pick lower granularity!") % cellcenter.Perp() % FilterPhiPos(cellcenter.Phi()) % cellcenter.Z()) << std::endl;
              exit(1);
            }
            // have to rotate this to the proper direction.
            fieldvec.RotateZ(FilterPhiPos(cellcenter.Phi()));  // rcc caution.  Does this rotation shift the sense of 'up'?
            (*field)->Set(j, i, k, fieldvec);
          }
        }
      }
    }
  }
  return;
}

void AnnularFieldSim::load_spacecharge(const std::string &filename, const std::string &histname, float zoffset, float chargescale, float cmscale, bool isChargeDensity)
{
  TFile *f = TFile::Open(filename.c_str());
  TH3 *scmap;
  f->GetObject(histname.c_str(), scmap);
  std::cout << "Loading spacecharge from '" << filename
            << "'.  Seeking histname '" << histname << "'" << std::endl;
  chargesourcename = filename + ":" + histname;
  load_spacecharge(scmap, zoffset, chargescale, cmscale, isChargeDensity, chargesourcename.c_str());
  f->Close();
  return;
}

void AnnularFieldSim::load_and_resample_spacecharge(int new_nphi, int new_nr, int new_nz, const std::string &filename, const std::string &histname, float zoffset, float chargescale, float cmscale, bool isChargeDensity)
{
  TFile *f = TFile::Open(filename.c_str());
  TH3 *scmap;
  f->GetObject(histname.c_str(), scmap);
  std::cout << "Resampling spacecharge from '" << filename
            << "'.  Seeking histname '" << histname << "'" << std::endl;
  chargesourcename = filename + ":" + histname;
  load_and_resample_spacecharge(new_nphi, new_nr, new_nz, scmap, zoffset, chargescale, cmscale, isChargeDensity);
  f->Close();
  return;
}

void AnnularFieldSim::load_and_resample_spacecharge(int new_nphi, int new_nr, int new_nz, TH3 *hist, float zoffset, float chargescale, float cmscale, bool isChargeDensity)
{
  // load spacecharge densities from a histogram, where scalefactor translates into local units of C/cm^3
  // and cmscale translate (hist coord) --> (hist position in cm)
  // noting that the histogram limits may differ from the simulation size, and have different granularity
  // hist is assumed/required to be x=phi, y=r, z=z
  // z offset 'drifts' the charge by that distance toward z=0.

  // Get dimensions of input
  float hrmin = hist->GetYaxis()->GetXmin();
  float hrmax = hist->GetYaxis()->GetXmax();
  float hphimin = hist->GetXaxis()->GetXmin();
  float hphimax = hist->GetXaxis()->GetXmax();
  float hzmin = hist->GetZaxis()->GetXmin();
  float hzmax = hist->GetZaxis()->GetXmax();

  // Get number of bins in each dimension
  int hrn = hist->GetNbinsY();
  int hphin = hist->GetNbinsX();
  int hzn = hist->GetNbinsZ();
  std::cout << boost::str(boost::format("From histogram, native r dimensions: %f<r<%f, hrn (Y axis)=%d") % hrmin % hrmax % hrn) << std::endl;
  std::cout << boost::str(boost::format("From histogram, native phi dimensions: %f<phi<%f, hrn (X axis)=%d") % hphimin % hphimax % hphin) << std::endl;
  std::cout << boost::str(boost::format("From histogram, native z dimensions: %f<z<%f, hrn (Z axis)=%d") % hzmin % hzmax % hzn) << std::endl;

  // do some computation of steps:
  float hrstep = (hrmax - hrmin) / hrn;
  float hphistep = (hphimax - hphimin) / hphin;
  float hzstep = (hzmax - hzmin) / hzn;
  float halfrstep = hrstep * 0.5;
  float halfphistep = hphistep * 0.5;
  float halfzstep = hzstep * 0.5;

  // All we have to do here is resample it so it fits into our expected grid, then pass it to the loader.

  // 1) convert the existing histogram to density if it isn't already:
  if (!isChargeDensity)
  {
    for (int i = 0; i < hphin; i++)
    {
      float phi = hphimin + hphistep * i;
      for (int j = 0; j < hrn; j++)
      {
        float r = hrmin + hrstep * j;
        for (int k = 0; k < hzn; k++)
        {
          float z = hzmin + hzstep * k;
          int bin = hist->FindBin(phi + halfphistep, r + halfrstep, z + halfzstep);
          double vol = hzstep * hphistep * (r + halfrstep) * hrstep;
          hist->SetBinContent(bin, hist->GetBinContent(bin) / vol);
        }
      }
    }
  }
  // I ought to consider a subtler approach that better preserves the overall integral, rather than just skipping the interpolation of the outer edges:
  TH3F *resampled = new TH3F("resampled", "resampled charge", new_nphi, hphimin, hphimax, new_nr, hrmin, hrmax, new_nz, hzmin, hzmax);
  float new_phistep = (hphimax - hphimin) / new_nphi;
  float new_rstep = (hrmax - hrmin) / new_nr;
  float new_zstep = (hzmax - hzmin) / new_nz;

  // 2 resample with interpolation on the interior, and with nearest-bin for the edge bins
  for (int i = 0; i < new_nphi; i++)
  {
    float phi = hphimin + new_phistep * (i + 0.5);                  // bin center of resampled
    float hphi = (phi - hphimin) / hphistep;                        // coord  in source hist
    bool phisafe = ((hphi - 0.75) > 0) && ((hphi + 0.75) < hphin);  // we're past the halfway point of the lowest bin, and short of the halfway point of the highest bin.
    for (int j = 0; j < new_nr; j++)
    {
      float r = hrmin + new_rstep * (j + 0.5);              // bin center of resampled
      float hr = (r - hrmin) / hrstep;                      // coord in source hist
      bool rsafe = ((hr - 0.5) > 0) && ((hr + 0.5) < hrn);  // we're past the halfway point of the lowest bin, and short of the halfway point of the highest bin.
      for (int k = 0; k < new_nz; k++)
      {
        float z = hzmin + new_zstep * (k + 0.5);              // bin center of resampled
        float hz = (z - hzmin) / hzstep;                      // coord in source hist
        bool zsafe = ((hz - 0.5) > 0) && ((hz + 0.5) < hzn);  // we're past the halfway point of the lowest bin, and short of the halfway point of the highest bin.
        // check if we need to do a bin lookup:
        if (phisafe && rsafe && zsafe)
        {
          // print_need_cout("resampling (%d,%d,%d) from (%2.2f/%d,%2.2f/%d,%2.2f/%d) p:(%1.1f<%1.1f<%1.1f) r:(%1.1f<%1.1f<%1.1f) z:(%1.1f<%1.1f<%1.1f)\n",i,j,k,hphi,hphin,hr,hrn,hz,hzn,hphimin,phi,hphimax,hrmin,r,hrmax,hzmin,z,hzmax);
          resampled->Fill(phi, r, z, hist->Interpolate(phi, r, z));  // leave as a density
        }
        else
        {
          int bin = hist->FindBin(phi, r, z);
          resampled->Fill(phi, r, z, hist->GetBinContent(bin));
        }
      }  // k (z loop)
    }    // j (r loop)
  }      // i (phi loop)
  std::cout << "here in load_and_resample_spacecharge, I am going to call load_spacecharge." << std::endl;
  load_spacecharge(resampled, zoffset, chargescale, cmscale, true);
}

void AnnularFieldSim::load_spacecharge(TH3 *hist, float zoffset, float chargescale, float cmscale, bool isChargeDensity, const std::string &inputchargestring)
{
  // new plan:  use ChargeMapReader:
  if (abs(zoffset) > 0.001)
  {
    std::cout << boost::str(boost::format("nonzero zoffset given (%E) but new spacecharge loader can't deal with that.  Failing.") % zoffset) << std::endl;
    assert(false);
  }
  if (isChargeDensity)
  {
    std::cout << "Input dataset is flagged as recording density not total charge, but new loader can't deal with that  Failing.." << std::endl;
    assert(false);
  }
  q->ReadSourceCharge(hist, cmscale, chargescale);

  chargestring = "SC loaded externally: " + inputchargestring + ".";
  return;
}

void AnnularFieldSim::load_digital_current(TH3 *hist, TH2 *gainHist, float chargescale, float cmscale, const std::string &inputchargestring)
{
  q->ReadSourceAdc(hist, gainHist, cmscale, chargescale);

  chargestring = "SC loaded externally: " + inputchargestring + ".";
  return;
}

void AnnularFieldSim::save_spacecharge(const std::string &filename)
{
  // save a histogram giving the spacecharge in our local coordinates, as we'll be using them..
  TH3F *hsc = new TH3F("hInternalSpaceCharge", "Internal Charge Histogram;phi(rad);r(cm);z(cm)", nphi, 0, phispan, nr, rmin, rmax, nz, zmin, zmax);
  for (int i = 0; i < nphi; i++)
  {
    float phi = 0 + step.Phi() * (i + 0.5);
    for (int j = 0; j < nr; j++)
    {
      float r = rmin + step.Perp() * (j + 0.5);
      for (int k = 0; k < nz; k++)
      {
        float z = zmin + step.Z() * (k + 0.5);
        hsc->Fill(phi, r, z, q->GetChargeAtPosition(r, phi, z));
        // old version: hsc->Fill(phi,r,z,q->Get(j,i,k));
      }
    }
  }
  hsc->SaveAs(filename.c_str());
  return;
}

void AnnularFieldSim::add_testcharge(float r, float phi, float z, float coulombs)
{
  q->AddChargeAtPosition(r, phi, z, coulombs * C);
  return;
  /*
  int rcell, phicell, zcell;

  //translate to which cell we're in:
  BoundsCase rb = GetRindexAndCheckBounds(r, &rcell);
  BoundsCase pb = GetPhiIndexAndCheckBounds(phi, &phicell);
  BoundsCase zb = GetZindexAndCheckBounds(z, &zcell);
  //really, we would like to confirm that the cell we find is in the charge map region, not the field map region.
  if (rb == BoundsCase::OutOfBounds || pb == BoundsCase::OutOfBounds || zb == BoundsCase::OutOfBounds)
  {
    print_need_cout("placing %f Coulombs at (r,p,z)=(%f,%f,%f).  Caution that that is out of the roi.\n", coulombs, r, phi, z);
  }
  if (rcell < 0 || phicell < 0 || zcell < 0)
  {
    print_need_cout("Tried placing %f Coulombs at (r,p,z)=(%f,%f,%f).  That is outside of the charge map region.\n", coulombs, r, phi, z);
    return;
  }

  q->Add(rcell, phicell, zcell, coulombs * C);
  if (lookupCase == HybridRes)
  {
    print_need_cout("write the code you didn't write earlier about reloading the lowres map.\n");
    exit(1);
  }

  return;
  */
}

/*
void AnnularFieldSim::populate_analytic_fieldmap(){
  //sum the E field at every point in the region of interest
  // remember that Efield uses relative indices
  print_need_cout("in pop_analytic_fieldmap, n=(%d,%d,%d)\n",nr,nphi,nz);

  TVector3 localF;//holder for the summed field at the current position.
  for (int ir=rmin_roi;ir<rmax_roi;ir++){
    for (int iphi=phimin_roi;iphi<phimax_roi;iphi++){
      for (int iz=zmin_roi;iz<zmax_roi;iz++){
        localF=aliceModel->E(GetCellCenter(ir, iphi, iz))+Eexternal->Get(ir-rmin_roi,iphi-phimin_roi,iz-zmin_roi);
        Efield->Set(ir-rmin_roi,iphi-phimin_roi,iz-zmin_roi,localF);
        //if (localF.Mag()>1e-9)
        if(debugFlag()) print_need_cout("%d: AnnularFieldSim::populate_analytic_fieldmap fieldmap@ (%d,%d,%d) mag=%f\n",__LINE__,ir,iphi,iz,localF.Mag());
      }
    }
  }
  return;
}
*/

void AnnularFieldSim::populate_fieldmap()
{
  // sum the E field at every point in the region of interest
  //  remember that Efield uses relative indices
  std::cout << boost::str(boost::format("in pop_fieldmap, n=(%d,%d,%d)") % nr % nphi % nz) << std::endl;

  std::cout << boost::str(boost::format("populating fieldmap for (%dx%dx%d) grid with (%dx%dx%d) source ") % nr_roi % nphi_roi % nz_roi % nr % nphi % nz) << std::endl;
  if (truncation_length > 0)
  {
    std::cout << boost::str(boost::format(" ==> truncating anything more than %d cells away") % truncation_length) << std::endl;
  }
  unsigned long long totalelements = nr_roi;
  totalelements *= nphi_roi;
  totalelements *= nz_roi;  // breaking up this multiplication prevents a 32bit math overflow
  unsigned long long percent = totalelements / 100 * debug_npercent;
  std::cout << boost::str(boost::format("total elements = %llu") % (totalelements * nr * nphi * nz)) << std::endl;

  int el = 0;

  TVector3 localF;  // holder for the summed field at the current position.
  for (int ir = rmin_roi; ir < rmax_roi; ir++)
  {
    for (int iphi = phimin_roi; iphi < phimax_roi; iphi++)
    {
      for (int iz = zmin_roi; iz < zmax_roi; iz++)
      {
        localF = sum_field_at(ir, iphi, iz);  // asks in global coordinates
        if (!(el % percent))
        {
          std::cout << boost::str(boost::format("populate_fieldmap %llu%%:  ") % ((uint64_t) (debug_npercent) *el / percent));
          std::cout << boost::str(boost::format("sum_field_at (ir=%d,iphi=%d,iz=%d) gives (%E,%E,%E)") % ir % iphi % iz % localF.X() % localF.Y() % localF.Z()) << std::endl;
        }
        el++;

        Efield->Set(ir - rmin_roi, iphi - phimin_roi, iz - zmin_roi, localF);  // sets in roi coordinates.
                                                                               // if (localF.Mag()>1e-9)
                                                                               // if(debugFlag()) print_need_cout("%d: AnnularFieldSim::populate_fieldmap fieldmap@ (%d,%d,%d) mag=%f\n",__LINE__,ir,iphi,iz,localF.Mag());
      }
    }
  }
  return;
}

void AnnularFieldSim::populate_lookup()
{
  // with 'f' being the position the field is being measured at, and 'o' being the position of the charge generating the field.
  // remember the 'f' part of Epartial uses relative indices.
  //   TVector3 (*f)[fx][fy][fz][ox][oy][oz]=field_;
  // print_need_cout("populating lookup for (%dx%dx%d)x(%dx%dx%d) grid\n",fx,fy,fz,ox,oy,oz);

  if (lookupCase == Full3D)
  {
    std::cout << "lookupCase==Full3D" << std::endl;

    populate_full3d_lookup();
  }
  else if (lookupCase == HybridRes)
  {
    std::cout << "lookupCase==HybridRes" << std::endl;
    populate_highres_lookup();
    populate_lowres_lookup();
  }
  else if (lookupCase == PhiSlice)
  {
    std::cout << "Populating lookup:  lookupCase==PhiSlice" << std::endl;
    populate_phislice_lookup();
  }
  else if (lookupCase == Analytic)
  {
    std::cout << "Populating lookup:  lookupCase==Analytic ===> skipping!" << std::endl;
  }
  else if (lookupCase == NoLookup)
  {
    std::cout << "Populating lookup:  lookupCase==NoLookup ===> skipping!" << std::endl;
  }
  else
  {
    exit(1);
  }
  return;
}

void AnnularFieldSim::populate_full3d_lookup()
{
  // with 'f' being the position the field is being measured at, and 'o' being the position of the charge generating the field.
  // remember the 'f' part of Epartial uses relative indices.
  //   TVector3 (*f)[fx][fy][fz][ox][oy][oz]=field_;
  std::cout << boost::str(boost::format("populating full lookup table for (%dx%dx%d)x(%dx%dx%d) grid") % (rmax_roi - rmin_roi) % (phimax_roi - phimin_roi) % (zmax_roi - zmin_roi) % nr % nphi % nz) << std::endl;
  unsigned long long totalelements = (rmax_roi - rmin_roi);
  totalelements *= (phimax_roi - phimin_roi);
  totalelements *= (zmax_roi - zmin_roi);
  totalelements *= nr;
  totalelements *= nphi;
  totalelements *= nz;  // breaking up this multiplication prevents a 32bit math overflow
  unsigned long long percent = totalelements / 100 * debug_npercent;
  std::cout << boost::str(boost::format("total elements = %llu") % totalelements) << std::endl;
  TVector3 at(1, 0, 0);
  TVector3 from(1, 0, 0);
  TVector3 zero(0, 0, 0);

  int el = 0;
  for (int ifr = rmin_roi; ifr < rmax_roi; ifr++)
  {
    for (int ifphi = phimin_roi; ifphi < phimax_roi; ifphi++)
    {
      for (int ifz = zmin_roi; ifz < zmax_roi; ifz++)
      {
        at = GetCellCenter(ifr, ifphi, ifz);
        for (int ior = 0; ior < nr; ior++)
        {
          for (int iophi = 0; iophi < nphi; iophi++)
          {
            for (int ioz = 0; ioz < nz; ioz++)
            {
              el++;
              if (!(el % percent))
              {
                std::cout << boost::str(boost::format("populate_full3d_lookup %d%%") % ((uint64_t) (debug_npercent) *el / percent)) << std::endl;
              }
              from = GetCellCenter(ior, iophi, ioz);

              //*f[ifx][ify][ifz][iox][ioy][ioz]=cacl_unit_field(at,from);
              // print_need_cout("calc_unit_field...\n");
              if (ifr == ior && ifphi == iophi && ifz == ioz)
              {
                Epartial->Set(ifr - rmin_roi, ifphi - phimin_roi, ifz - zmin_roi, ior, iophi, ioz, zero);
              }
              else
              {
                Epartial->Set(ifr - rmin_roi, ifphi - phimin_roi, ifz - zmin_roi, ior, iophi, ioz, calc_unit_field(at, from));
              }
            }
          }
        }
      }
    }
  }
  return;
}

void AnnularFieldSim::populate_highres_lookup()
{
  TVector3 at(1, 0, 0);
  TVector3 from(1, 0, 0);
  TVector3 zero(0, 0, 0);

  // populate_highres_lookup();
  int r_highres_dist = (nr_high - 1) / 2;
  int phi_highres_dist = (nphi_high - 1) / 2;
  int z_highres_dist = (nz_high - 1) / 2;

  // number of fbins contained in the 26 weirdly-shaped edge regions (and one center region which we won't use)
  static int nfbinsin[3][3][3];
  for (auto &i : nfbinsin)
  {
    for (auto &j : i)
    {
      for (int &k : j)
      {
        k = 0;  // we could count total volume, but without knowing the charge prior, it's not clear that'd be /better/
      }
    }
  }

  // todo: if this runs too slowly, I can do geometry instead of looping over all the cells that are possibly in range

  // loop over all the f-bins in the roi:
  TVector3 currentf, newf;  // the averaged field vector without, and then with the new contribution from the f-bin being considered.
  for (int ifr = rmin_roi; ifr < rmax_roi; ifr++)
  {
    int r_parentlow = floor((ifr - r_highres_dist) / (r_spacing * 1.0));       // l-bin partly enclosed in our high-res region
    int r_parenthigh = floor((ifr + r_highres_dist) / (r_spacing * 1.0)) + 1;  // definitely not enclosed in our high-res region
    int r_startpoint = r_parentlow * r_spacing;                                // the first f-bin of the lowest-r f-bin that our h-region touches.  COuld be less than zero.
    int r_endpoint = r_parenthigh * r_spacing;                                 // the first f-bin of the lowest-r l-bin after that that our h-region does not touch.  could be larger than max.

    for (int ifphi = phimin_roi; ifphi < phimax_roi; ifphi++)
    {
      int phi_parentlow = floor(FilterPhiIndex(ifphi - phi_highres_dist) / (phi_spacing * 1.0));  // note this may have wrapped around
      bool phi_parentlow_wrapped = (ifphi - phi_highres_dist < 0);
      int phi_startpoint = phi_parentlow * phi_spacing;  // the first f-bin of the lowest-z f-bin that our h-region touches.
      if (phi_parentlow_wrapped)
      {
        phi_startpoint -= nphi;  // if we wrapped, re-wrap us so we're negative again
      }

      int phi_parenthigh = floor(FilterPhiIndex(ifphi + phi_highres_dist) / (phi_spacing * 1.0)) + 1;  // note that this may have wrapped around
      bool phi_parenthigh_wrapped = (ifphi + phi_highres_dist >= nphi);
      int phi_endpoint = phi_parenthigh * phi_spacing;
      if (phi_parenthigh_wrapped)
      {
        phi_endpoint += nphi;  // if we wrapped, re-wrap us so we're larger than nphi again.  We use these relative coords to determine the position relative to the center of our h-region.
      }

      for (int ifz = zmin_roi; ifz < zmax_roi; ifz++)
      {
        // if(debugFlag()) print_need_cout("%d: AnnularFieldSim::populate_highres_lookup icell=(%d,%d,%d)\n",__LINE__,ifr,ifphi,ifz);

        int z_parentlow = floor((ifz - z_highres_dist) / (z_spacing * 1.0));
        int z_parenthigh = floor((ifz + z_highres_dist) / (z_spacing * 1.0)) + 1;
        int z_startpoint = z_parentlow * z_spacing;  // the first f-bin of the lowest-z f-bin that our h-region touches.
        int z_endpoint = z_parenthigh * z_spacing;   // the first f-bin of the lowest-z l-bin after that that our h-region does not touch.

        // our 'at' position, in global coords:
        at = GetCellCenter(ifr, ifphi, ifz);
        // define the farthest-away parent l-bin cells we're dealing with here:
        // note we're still in absolute coordinates

        // for every f-bin in the l-bins we're dealing with, figure out which relative highres bin it's in, and average the field vector into that bin's vector
        // note most of these relative bins have exactly one f-bin in them.  It's only the edges that can get more.
        // note this is a running average:  Anew=(Aold*Nold+V)/(Nold+1) and so on.
        // note also that we automatically skip f-bins that would've been out of the valid overall volume.
        for (int ir = r_startpoint; ir < r_endpoint; ir++)
        {
          // skip parts that are out of range:
          // could speed this up by moving this into the definition of start and endpoint.
          if (ir < 0)
          {
            ir = 0;
          }
          if (ir >= nr)
          {
            break;
          }

          int rbin = (ir - ifr) + r_highres_dist;  // zeroth bin when we're at max distance below, etc.
          int rcell = 1;
          if (rbin <= 0)
          {
            rbin = 0;
            rcell = 0;
          }
          if (rbin >= nr_high)
          {
            rbin = nr_high - 1;
            rcell = 2;
          }

          for (int iphi = phi_startpoint; iphi < phi_endpoint; iphi++)
          {
            // no phi out-of-range checks since it's circular, but we provide ourselves a filtered version:
            int phiFilt = FilterPhiIndex(iphi);
            int phibin = (iphi - ifphi) + phi_highres_dist;
            int phicell = 1;
            if (phibin <= 0)
            {
              phibin = 0;
              phicell = 0;
            }
            if (phibin >= nphi_high)
            {
              phibin = nphi_high - 1;
              phicell = 2;
            }
            for (int iz = z_startpoint; iz < z_endpoint; iz++)
            {
              if (iz < 0)
              {
                iz = 0;
              }
              if (iz >= nz)
              {
                break;
              }
              int zbin = (iz - ifz) + z_highres_dist;
              int zcell = 1;
              if (zbin <= 0)
              {
                zbin = 0;
                zcell = 0;
              }
              if (zbin >= nz_high)
              {
                zbin = nz_high - 1;
                zcell = 2;
              }
              //'from' is in absolute coordinates
              from = GetCellCenter(ir, phiFilt, iz);

              nfbinsin[rcell][phicell][zcell]++;
              int nf = nfbinsin[rcell][phicell][zcell];
              // coordinates relative to the region of interest:
              int ir_rel = ifr - rmin_roi;
              int iphi_rel = ifphi - phimin_roi;
              int iz_rel = ifz - zmin_roi;
              if (zcell != 1 || rcell != 1 || phicell != 1)
              {
                // we're not in the center, so deal with our weird shapes by averaging:
                // but Epartial is in coordinates relative to the roi
                if (iphi_rel < 0)
                {
                  std::cout << boost::str(boost::format("%d: Getting with phi=%d") % __LINE__ % iphi_rel) << std::endl;
                }
                currentf = Epartial_highres->Get(ir_rel, iphi_rel, iz_rel, rbin, phibin, zbin);
                // to keep this as the average, we multiply what's there back to its initial summed-but-not-divided value
                // then add our new value, and the divide the new sum by the total number of cells
                newf = (currentf * (nf - 1) + calc_unit_field(at, from)) * (1 / (nf * 1.0));
                Epartial_highres->Set(ir_rel, iphi_rel, iz_rel, rbin, phibin, zbin, newf);
              }
              else
              {
                // we're in the center cell, which means any f-bin that's not on the outer edge of our region:
                // calc_unit_field will return zero when at=from, so the center will be automatically zero here.
                if (ifr == rbin && ifphi == phibin && ifz == zbin)
                {
                  Epartial_highres->Set(ir_rel, iphi_rel, iz_rel, rbin, phibin, zbin, zero);
                }
                else
                {  // for extra carefulness, only calc the field if it's not self-to-self.
                  newf = calc_unit_field(at, from);
                  Epartial_highres->Set(ir_rel, iphi_rel, iz_rel, rbin, phibin, zbin, newf);
                }
              }
            }
          }
        }
      }
    }
  }
  return;
}

void AnnularFieldSim::populate_lowres_lookup()
{
  TVector3 at(1, 0, 0);
  TVector3 from(1, 0, 0);
  TVector3 zero(0, 0, 0);
  int fphi_low, fphi_high, fz_low, fz_high;             // edges of the outer l-bin
  int r_low, r_high, phi_low, phi_high, z_low, z_high;  // edges of the inner l-bin

  // todo:  add in handling if roi_low is wrap-around in phi
  for (int ifr = rmin_roi_low; ifr < rmax_roi_low; ifr++)
  {
    int fr_low = ifr * r_spacing;
    int fr_high = fr_low + r_spacing - 1;
    if (fr_high >= nr)
    {
      fr_high = nr - 1;
    }
    for (int ifphi = phimin_roi_low; ifphi < phimax_roi_low; ifphi++)
    {
      fphi_low = ifphi * phi_spacing;
      fphi_high = fphi_low + phi_spacing - 1;
      if (fphi_high >= nphi)
      {
        fphi_high = nphi - 1;  // if our phi l-bins aren't evenly spaced, we need to catch that here.
      }
      for (int ifz = zmin_roi_low; ifz < zmax_roi_low; ifz++)
      {
        fz_low = ifz * z_spacing;
        fz_high = fz_low + z_spacing - 1;
        if (fz_high >= nz)
        {
          fz_high = nz - 1;
        }
        at = GetGroupCellCenter(fr_low, fr_high, fphi_low, fphi_high, fz_low, fz_high);
        // print_need_cout("ifr=%d, rlow=%d,rhigh=%d,r_spacing=%d\n",ifr,r_low,r_high,r_spacing);
        // if(debugFlag())	  print_need_cout("%d: AnnularFieldSim::populate_lowres_lookup icell=(%d,%d,%d)\n",__LINE__,ifr,ifphi,ifz);

        for (int ior = 0; ior < nr_low; ior++)
        {
          r_low = ior * r_spacing;
          r_high = r_low + r_spacing - 1;
          int ir_rel = ifr - rmin_roi_low;

          if (r_high >= nr)
          {
            r_high = nr - 1;
          }
          for (int iophi = 0; iophi < nphi_low; iophi++)
          {
            phi_low = iophi * phi_spacing;
            phi_high = phi_low + phi_spacing - 1;
            if (phi_high >= nphi)
            {
              phi_high = nphi - 1;
            }
            int iphi_rel = ifphi - phimin_roi_low;
            for (int ioz = 0; ioz < nz_low; ioz++)
            {
              z_low = ioz * z_spacing;
              z_high = z_low + z_spacing - 1;
              if (z_high >= nz)
              {
                z_high = nz - 1;
              }
              int iz_rel = ifz - zmin_roi_low;
              from = GetGroupCellCenter(r_low, r_high, phi_low, phi_high, z_low, z_high);

              if (ifr == ior && ifphi == iophi && ifz == ioz)
              {
                Epartial_lowres->Set(ir_rel, iphi_rel, iz_rel, ior, iophi, ioz, zero);
              }
              else
              {  // for extra carefulness, only calc the field if it's not self-to-self.
                Epartial_lowres->Set(ir_rel, iphi_rel, iz_rel, ior, iophi, ioz, calc_unit_field(at, from));
              }

              //*f[ifx][ify][ifz][iox][ioy][ioz]=cacl_unit_field(at,from);
              // print_need_cout("calc_unit_field...\n");
              // this calc's okay.
            }
          }
        }
      }
    }
  }
  return;
}

void AnnularFieldSim::populate_phislice_lookup()
{
  // with 'f' being the position the field is being measured at, and 'o' being the position of the charge generating the field.
  // remember the 'f' part of Epartial uses relative indices.
  //   TVector3 (*f)[fx][fy][fz][ox][oy][oz]=field_;
  std::cout << boost::str(boost::format("populating phislice  lookup for (%dx%dx%d)x(%dx%dx%d) grid") % nr_roi % 1 % nz_roi % nr % nphi % nz) << std::endl;
  unsigned long long totalelements = nr;  // nr*nphi*nz*nr_roi*nz_roi
  totalelements *= nphi;
  totalelements *= nz;
  totalelements *= nr_roi;
  totalelements *= nz_roi;  // breaking up this multiplication prevents a 32bit math overflow
  unsigned long long percent = totalelements / 100 * debug_npercent;
  std::cout << boost::str(boost::format("total elements = %llu") % totalelements) << std::endl;
  TVector3 at(1, 0, 0);
  TVector3 from(1, 0, 0);
  TVector3 zero(0, 0, 0);

  int el = 0;
  for (int ifr = rmin_roi; ifr < rmax_roi; ifr++)
  {
    for (int ifz = zmin_roi; ifz < zmax_roi; ifz++)
    {
      at = GetCellCenter(ifr, 0, ifz);
      for (int ior = 0; ior < nr; ior++)
      {
        for (int iophi = 0; iophi < nphi; iophi++)
        {
          for (int ioz = 0; ioz < nz; ioz++)
          {
            el++;
            from = GetCellCenter(ior, iophi, ioz);
            //*f[ifx][ify][ifz][iox][ioy][ioz]=cacl_unit_field(at,from);
            // print_need_cout("calc_unit_field...\n");
            if (ifr == ior && 0 == iophi && ifz == ioz)
            {
              if (!(el % percent))
              {
                std::cout << boost::str(boost::format("populate_phislice_lookup %llu%%:  ") % ((uint64_t) (debug_npercent) *el / percent));
                std::cout << boost::str(boost::format("self-to-self is zero (ir=%d,iphi=%d,iz=%d) to (or=%d,ophi=0,oz=%d) gives (%E,%E,%E)") % ior % iophi % ioz % ifr % ifz % zero.X() % zero.Y() % zero.Z()) << std::endl;
              }
              Epartial_phislice->Set(ifr - rmin_roi, 0, ifz - zmin_roi, ior, iophi, ioz, zero);
            }
            else
            {
              TVector3 unitf = calc_unit_field(at, from);
              if (true)
              {
                if (!(el % percent))
                {
                  std::cout << boost::str(boost::format("populate_phislice_lookup %llu%%:  ") % ((uint64_t) (debug_npercent) *el / percent));
                  std::cout << boost::str(boost::format("calc_unit_field (ir=%d,iphi=%d,iz=%d) to (or=%d,ophi=0,oz=%d) gives (%E,%E,%E)") % ior % iophi % ioz % ifr % ifz % unitf.X() % unitf.Y() % unitf.Z()) << std::endl;
                }
              }

              Epartial_phislice->Set(ifr - rmin_roi, 0, ifz - zmin_roi, ior, iophi, ioz, unitf);  // the origin phi is relative to zero anyway.
            }
          }
        }
      }
    }
  }
  return;
}

void AnnularFieldSim::load_phislice_lookup(const std::string &sourcefile)
{
  std::cout << boost::str(boost::format("loading phislice  lookup for (%dx%dx%d)x(%dx%dx%d) grid from %s") % nr_roi % 1 % nz_roi % nr % nphi % nz % sourcefile) << std::endl;
  unsigned long long totalelements = nr;  // nr*nphi*nz*nr_roi*nz_roi
  totalelements *= nphi;
  totalelements *= nz;
  totalelements *= nr_roi;
  totalelements *= nz_roi;  // breaking up this multiplication prevents a 32bit math overflow
  unsigned long long percent = totalelements / 100 * debug_npercent;
  std::cout << boost::str(boost::format("total elements = %llu") % totalelements) << std::endl;

  TFile *input = TFile::Open(sourcefile.c_str(), "READ");
  TTree *tInfo;
  input->GetObject("info", tInfo);
  assert(tInfo);

  float file_rmin, file_rmax, file_zmin, file_zmax;
  int file_rmin_roi, file_rmax_roi, file_zmin_roi, file_zmax_roi;
  int file_nr, file_np, file_nz;
  tInfo->SetBranchAddress("rmin", &file_rmin);
  tInfo->SetBranchAddress("rmax", &file_rmax);
  tInfo->SetBranchAddress("zmin", &file_zmin);
  tInfo->SetBranchAddress("zmax", &file_zmax);
  tInfo->SetBranchAddress("rmin_roi_index", &file_rmin_roi);
  tInfo->SetBranchAddress("rmax_roi_index", &file_rmax_roi);
  tInfo->SetBranchAddress("zmin_roi_index", &file_zmin_roi);
  tInfo->SetBranchAddress("zmax_roi_index", &file_zmax_roi);
  tInfo->SetBranchAddress("nr", &file_nr);
  tInfo->SetBranchAddress("nphi", &file_np);
  tInfo->SetBranchAddress("nz", &file_nz);
  tInfo->GetEntry(0);

  std::cout << "param\tobj\tfile" << std::endl;
  ;
  std::cout << boost::str(boost::format("nr\t%d\t%d") % nr % file_nr) << std::endl;
  std::cout << boost::str(boost::format("np\t%d\t%d") % nphi % file_np) << std::endl;
  std::cout << boost::str(boost::format("nz\t%d\t%d") % nz % file_nz) << std::endl;
  std::cout << boost::str(boost::format("rmin\t%2.2f\t%2.2f") % rmin % file_rmin) << std::endl;
  std::cout << boost::str(boost::format("rmax\t%2.2f\t%2.2f") % rmax % file_rmax) << std::endl;
  std::cout << boost::str(boost::format("zmin\t%2.2f\t%2.2f") % zmin % file_zmin) << std::endl;
  std::cout << boost::str(boost::format("zmax\t%2.2f\t%2.2f") % zmax % file_zmax) << std::endl;
  std::cout << boost::str(boost::format("rmin_roi\t%d\t%d") % rmin_roi % file_rmin_roi) << std::endl;
  std::cout << boost::str(boost::format("rmax_roi\t%d\t%d") % rmax_roi % file_rmax_roi) << std::endl;
  std::cout << boost::str(boost::format("zmin_roi\t%d\t%d") % zmin_roi % file_zmin_roi) << std::endl;
  std::cout << boost::str(boost::format("zmax_roi\t%d\t%d") % zmax_roi % file_zmax_roi) << std::endl;

  if (file_rmin != rmin || file_rmax != rmax ||
      file_zmin != zmin || file_zmax != zmax ||
      file_rmin_roi != rmin_roi || file_rmax_roi != rmax_roi ||
      file_zmin_roi != zmin_roi || file_zmax_roi != zmax_roi ||
      file_nr != nr || file_np != nphi || file_nz != nz)
  {
    std::cout << "file parameters do not match fieldsim parameters:" << std::endl;
    ;

    exit(1);
  }

  TTree *tLookup;
  input->GetObject("phislice", tLookup);
  assert(tLookup);
  int ior, ifr, iophi, ioz, ifz;
  TVector3 *unitf = nullptr;
  tLookup->SetBranchAddress("ir_source", &ior);
  tLookup->SetBranchAddress("ir_target", &ifr);
  tLookup->SetBranchAddress("ip_source", &iophi);
  // always zero: tLookup->SetBranchAddress("ip_target",&ifphi);
  tLookup->SetBranchAddress("iz_source", &ioz);
  tLookup->SetBranchAddress("iz_target", &ifz);
  tLookup->SetBranchAddress("Evec", &unitf);  // assume field has units V/(C*cm)

  int el = 0;
  std::cout << boost::str(boost::format("%s has %lld entries") % sourcefile % tLookup->GetEntries());
  for (int i = 0; i < (int) totalelements; i++)
  {
    el++;
    tLookup->GetEntry(i);
    // print_need_cout("loading i=%d\n",i);
    Epartial_phislice->Set(ifr - rmin_roi, 0, ifz - zmin_roi, ior, iophi, ioz, (*unitf) * (-1.0) * (V / (cm * C)));  // load assuming field has units V/(C*cm), which is how we save it.
    // note that we save the gradient terms, not the field, hence we need to multiply by (-1.0)
    if (!(el % percent))
    {
      std::cout << boost::str(boost::format("load_phislice_lookup %llu%%:  ") % ((uint64_t) (debug_npercent) *el / percent));
      std::cout << boost::str(boost::format("field from (ir=%d,iphi=%d,iz=%d) to (or=%d,ophi=0,oz=%d) is (%E,%E,%E)") % ior % iophi % ioz % ifr % ifz % unitf->X() % unitf->Y() % unitf->Z()) << std::endl;
    }
  }

  input->Close();
  return;
}

void AnnularFieldSim::save_phislice_lookup(const std::string &destfile)
{
  std::cout << boost::str(boost::format("saving phislice  lookup for (%dx%dx%d)x(%dx%dx%d) grid to %s") % nr_roi % 1 % nz_roi % nr % nphi % nz % destfile) << std::endl;
  unsigned long long totalelements = nr;  // nr*nphi*nz*nr_roi*nz_roi
  totalelements *= nphi;
  totalelements *= nz;
  totalelements *= nr_roi;
  totalelements *= nz_roi;  // breaking up this multiplication prevents a 32bit math overflow
  unsigned long long percent = totalelements / 100 * debug_npercent;
  std::cout << boost::str(boost::format("total elements = %llu") % totalelements) << std::endl;

  TFile *output = TFile::Open(destfile.c_str(), "RECREATE");
  output->cd();

  TTree *tInfo = new TTree("info", "Information about Lookup Table");
  tInfo->Branch("rmin", &rmin);
  tInfo->Branch("rmax", &rmax);
  tInfo->Branch("zmin", &zmin);
  tInfo->Branch("zmax", &zmax);
  tInfo->Branch("rmin_roi_index", &rmin_roi);
  tInfo->Branch("rmax_roi_index", &rmax_roi);
  tInfo->Branch("zmin_roi_index", &zmin_roi);
  tInfo->Branch("zmax_roi_index", &zmax_roi);
  tInfo->Branch("nr", &nr);
  tInfo->Branch("nphi", &nphi);
  tInfo->Branch("nz", &nz);
  tInfo->Fill();
  std::cout << "info tree built." << std::endl;
  ;

  TTree *tLookup = new TTree("phislice", "Phislice Lookup Table");
  int ior, ifr, iophi, ioz, ifz;
  TVector3 unitf;
  tLookup->Branch("ir_source", &ior);
  tLookup->Branch("ir_target", &ifr);
  tLookup->Branch("ip_source", &iophi);
  // always zero: tLookup->Branch("ip_target",&ifphi);
  tLookup->Branch("iz_source", &ioz);
  tLookup->Branch("iz_target", &ifz);
  tLookup->Branch("Evec", &unitf);
  std::cout << "lookup tree built." << std::endl;
  ;

  int el = 0;
  for (ifr = rmin_roi; ifr < rmax_roi; ifr++)
  {
    for (ifz = zmin_roi; ifz < zmax_roi; ifz++)
    {
      for (ior = 0; ior < nr; ior++)
      {
        for (iophi = 0; iophi < nphi; iophi++)
        {
          for (ioz = 0; ioz < nz; ioz++)
          {
            el++;
            unitf = Epartial_phislice->Get(ifr - rmin_roi, 0, ifz - zmin_roi, ior, iophi, ioz) * (-1 / (V / (C * cm)));  // save in units of V/(C*cm) note that we introduce a -1 here for legcy reasons.
            if (true)
            {
              if (!(el % percent))
              {
                std::cout << boost::str(boost::format("save_phislice_lookup %llu%%:  ") % ((uint64_t) (debug_npercent) *el / percent));
                std::cout << boost::str(boost::format("field from (ir=%d,iphi=%d,iz=%d) to (or=%d,ophi=0,oz=%d) is (%E,%E,%E)") % ior % iophi % ioz % ifr % ifz % unitf.X() % unitf.Y() % unitf.Z()) << std::endl;
              }
            }

            tLookup->Fill();
          }
        }
      }
    }
  }
  output->cd();
  tInfo->Write();
  tLookup->Write();
  // output->Write();
  output->Close();
  return;
}

void AnnularFieldSim::setFlatFields(float B, float E)
{
  // these only cover the roi, but since we address them flat, we don't need to know that here.
  std::cout << boost::str(boost::format("AnnularFieldSim::setFlatFields(B=%f Tesla,E=%f V/cm)") % B % E) << std::endl;
  std::cout << boost::str(boost::format("lengths:  Eext=%d, Bfie=%d") % Eexternal->Length() % Bfield->Length()) << std::endl;
  Efieldname = "E:Flat:" + std::to_string(E);
  Bfieldname = "B:Flat:" + std::to_string(B);

  Enominal = E * (V / cm);
  Bnominal = B * Tesla;
  for (int i = 0; i < Eexternal->Length(); i++)
  {
    Eexternal->GetFlat(i)->SetXYZ(0, 0, Enominal);
  }
  for (int i = 0; i < Bfield->Length(); i++)
  {
    Bfield->GetFlat(i)->SetXYZ(0, 0, Bnominal);
  }
  UpdateOmegaTau();
  return;
}

TVector3 AnnularFieldSim::sum_field_at(int r, int phi, int z)
{
  // sum the E field over all nr by ny by nz cells of sources, at the global coordinate position r,phi,z.
  // note the specific position in Epartial is in relative coordinates.
  //  if(debugFlag()) print_need_cout("%d: AnnularFieldSim::sum_field_at(r=%d,phi=%d, z=%d)\n",__LINE__,r,phi,z);

  /*
    for the near future:
  TVector3 sum=(sum_local_field_at(r, phi, z)
                + sum_nonlocal_field_at(r,phi,z)
                + Escale*Eexternal->Get(r-rmin_roi,phi-phimin_roi,z-zmin_roi));
  return sum;
  */

  TVector3 sum(0, 0, 0);
  if (lookupCase == Full3D)
  {
    sum += sum_full3d_field_at(r, phi, z);
  }
  else if (lookupCase == HybridRes)
  {
    sum += sum_local_field_at(r, phi, z);
    sum += sum_nonlocal_field_at(r, phi, z);
  }
  else if (lookupCase == PhiSlice)
  {
    sum += sum_phislice_field_at(r, phi, z);
  }
  else if (lookupCase == Analytic)
  {
    sum += aliceModel->E(GetCellCenter(r, phi, z));
  }
  else if (lookupCase == NoLookup)
  {
    // do nothing.  We are forcibly assuming E from spacecharge is zero everywhere.
  }
  sum += Eexternal->Get(r - rmin_roi, phi - phimin_roi, z - zmin_roi);
  if (debugFlag())
  {
    std::cout << boost::str(boost::format("summed field at (%d,%d,%d)=(%f,%f,%f)") % r % phi % z % sum.X() % sum.Y() % sum.Z()) << std::endl;
  }

  return sum;
}

TVector3 AnnularFieldSim::sum_full3d_field_at(int r, int phi, int z)
{
  // sum the E field over all nr by ny by nz cells of sources, at the specific position r,phi,z.
  // note the specific position in Epartial is in relative coordinates.
  // print_need_cout("AnnularFieldSim::sum_field_at(r=%d,phi=%d, z=%d)\n",r,phi,z);
  TVector3 sum(0, 0, 0);
  float rdist, phidist, zdist, remdist;
  for (int ir = 0; ir < nr; ir++)
  {
    if (truncation_length > 0)
    {
      rdist = abs(ir - r);
      remdist = sqrt(truncation_length * truncation_length - rdist * rdist);
      if (remdist < 0)
      {
        continue;  // skip if we're too far away
      }
    }
    for (int iphi = 0; iphi < nphi; iphi++)
    {
      if (truncation_length > 0)
      {
        phidist = fmin(abs(iphi - phi), abs(abs(iphi - phi) - nphi));  // think about this in phi... rcc food for thought.
        remdist = sqrt(truncation_length * truncation_length - phidist * phidist);
        if (remdist < 0)
        {
          continue;  // skip if we're too far away
        }
      }
      for (int iz = 0; iz < nz; iz++)
      {
        if (truncation_length > 0)
        {
          zdist = abs(iz - z);
          remdist = sqrt(truncation_length * truncation_length - zdist * zdist);
          if (remdist < 0)
          {
            continue;  // skip if we're too far away
          }
        }
        // sum+=*partial[x][phi][z][ix][iphi][iz] * *q[ix][iphi][iz];
        if (r == ir && phi == iphi && z == iz)
        {
          continue;  // dont' compute self-to-self field.
        }
        sum += Epartial->Get(r - rmin_roi, phi - phimin_roi, z - zmin_roi, ir, iphi, iz) * q->GetChargeInBin(ir, iphi, iz);
      }
    }
  }
  // print_need_cout("summed field at (%d,%d,%d)=(%f,%f,%f)\n",x,y,z,sum.X(),sum.Y(),sum.Z());
  return sum;
}

TVector3 AnnularFieldSim::sum_local_field_at(int r, int phi, int z)
{
  // do the summation of the inner high-resolution region charges:
  //
  //  bin 0  1 2 ...  n-2 n-1
  //   . . .|.|.|.|.|.|.|. . .
  //        | | | | | | |
  //   . . .|.|.|.|.|.|.|. . .
  //   -----+-+-+-+-+-+-+-----
  //   . . .|.|.|.|.|.|.|. . .
  //   -----+-+-+-+-+-+-+-----
  //   . . .|.|.|.|.|.|.|. . .
  //   -----+-+-+-+-+-+-+-----
  //   . . .|.|.|.|.|.|.|. . .
  //   -----+-+-+-+-+-+-+-----
  //   . . .|.|.|.|.|.|.|. . .
  //        | | | | | | |
  //   . . .|.|.|.|.|.|.|. . .
  //
  //

  int r_highres_dist = (nr_high - 1) / 2;
  int phi_highres_dist = (nphi_high - 1) / 2;
  int z_highres_dist = (nz_high - 1) / 2;

  int r_parentlow = floor((r - r_highres_dist) / (r_spacing * 1.0));       // partly enclosed in our high-res region
  int r_parenthigh = floor((r + r_highres_dist) / (r_spacing * 1.0)) + 1;  // definitely not enclosed in our high-res region
  int phi_parentlow = floor((phi - phi_highres_dist) / (phi_spacing * 1.0));
  int phi_parenthigh = floor((phi + phi_highres_dist) / (phi_spacing * 1.0)) + 1;  // note that this can be bigger than nphi!  We keep track of that.
  int z_parentlow = floor((z - z_highres_dist) / (z_spacing * 1.0));
  int z_parenthigh = floor((z + z_highres_dist) / (z_spacing * 1.0)) + 1;
  std::cout << boost::str(boost::format("AnnularFieldSim::sum_local_field_at parents: rlow=%d,philow=%d,zlow=%d,rhigh=%d,phihigh=%d,zhigh=%d") % r_parentlow % phi_parentlow % z_parentlow % r_parenthigh % phi_parenthigh % z_parenthigh) << std::endl;

  // zero our current qlocal holder:
  for (int i = 0; i < q_local->Length(); i++)
  {
    *(q_local->GetFlat(i)) = 0;
  }

  // get the charge involved in the local highres block:
  for (int ir = r_parentlow * r_spacing; ir < r_parenthigh * r_spacing; ir++)
  {
    if (ir < 0)
    {
      ir = 0;
    }
    if (ir >= nr)
    {
      break;
    }
    int rbin = (ir - r) + r_highres_dist;  // index in our highres locale.  zeroth bin when we're at max distance below, etc.
    if (rbin < 0)
    {
      rbin = 0;
    }
    if (rbin >= nr_high)
    {
      rbin = nr_high - 1;
    }
    for (int iphi = phi_parentlow * phi_spacing; iphi < phi_parenthigh * phi_spacing; iphi++)
    {
      // no phi range checks since it's circular.
      int phiFilt = FilterPhiIndex(iphi);
      int phibin = (iphi - phi) + phi_highres_dist;
      if (phibin < 0)
      {
        phibin = 0;
      }
      if (phibin >= nphi_high)
      {
        phibin = nphi_high - 1;
      }
      for (int iz = z_parentlow * z_spacing; iz < z_parenthigh * z_spacing; iz++)
      {
        if (iz < 0)
        {
          iz = 0;
        }
        if (iz >= nz)
        {
          break;
        }
        // if(debugFlag()) print_need_cout("%d: AnnularFieldSim::sum_field_at, reading q in f-bin at(r=%d,phi=%d, z=%d)\n",__LINE__,ir,iphi,iz);

        int zbin = (iz - z) + z_highres_dist;
        if (zbin < 0)
        {
          zbin = 0;
        }
        if (zbin >= nz_high)
        {
          zbin = nz_high - 1;
        }
        // print_need_cout("filtering in local highres block\n");
        q_local->Add(rbin, phibin, zbin, q->GetChargeInBin(ir, phiFilt, iz));
        // print_need_cout("done filtering in local highres block\n");
      }
    }
  }

  // now q_highres is up to date for our current cell of interest.
  // start building our full sum by scaling the local lookup table by q.
  // note that the lookup table needs to have already accounted for cell centers.

  TVector3 sum(0, 0, 0);

  // note that Epartial_highres returns zero if we're outside of the global region.  q_local will also be zero there.
  // these are loops over the position in the epartial highres grid, so relative to the point in question:
  // reminder: variables r, phi, and z are global f-bin indices.

  // assuming highres correctly gives a zero when asked for the middle element in each of the last three indices,
  // and assuming q_local has no contribution from regions outside the global bounds, this is correct:
  for (int ir = 0; ir < nr_high; ir++)
  {
    for (int iphi = 0; iphi < nphi_high; iphi++)
    {
      for (int iz = 0; iz < nz_high; iz++)
      {
        // first three are relative to the roi, last three are relative to the point in the first three.  ooph.
        if (phi - phimin_roi < 0)
        {
          std::cout << boost::str(boost::format("%d: Getting with phi=%d") % __LINE__ % (phi - phimin_roi)) << std::endl;
        }
        sum += Epartial_highres->Get(r - rmin_roi, phi - phimin_roi, z - zmin_roi, ir, iphi, iz) * q_local->Get(ir, iphi, iz);
      }
    }
  }

  return sum;
}

TVector3 AnnularFieldSim::sum_nonlocal_field_at(int r, int phi, int z)
{
  // now we look for our low_res contribution, which will be the interpolated summed field from the eight blocks closest to our f-bin of interest:

  // find our interpolated position between the eight nearby lowres cells:
  // this is similar to the interpolated integral stuff, except we're at integer steps inside integer blocks
  // and we have z as well, now.
  bool skip[] = {false, false, false, false, false, false, false, false};

  float r0 = r / (1.0 * r_spacing) - 0.5 - rmin_roi_low;  // the position in r, in units of r_spacing, starting from the center of the 0th l-bin in the roi.
  int r0i = floor(r0);                                    // the integer portion of the position. -- what center is below our position?
  float r0d = r0 - r0i;                                   // the decimal portion of the position. -- how far past center are we?
  // instead of listing all eight, I'll address these as i%2, (i/2)%2 and (i/4)%2 to avoid typos
  int ri[2];  // the r position of the eight cell centers to consider.
  ri[0] = r0i;
  ri[1] = r0i + 1;
  float rw[2];      // the weight of that cell, linear in distance from the center of it
  rw[0] = 1 - r0d;  // 1 when we're on it, 0 when we're at the other one.
  rw[1] = r0d;      // 1 when we're on it, 0 when we're at the other one

  // zero out if the requested l-bin is out of the roi (happens if we're close to the outer than the inner edge of the outermost l-bin)
  if (ri[0] < 0)
  {
    for (int i = 0; i < 8; i++)
    {
      if ((i / 4) % 2 == 0)
      {
        skip[i] = true;  // don't handle contributions from ri[0].
      }
    }
    rw[1] = 1;  // and weight like we're dead-center on the outer cells.
  }
  if (ri[1] >= nr_roi_low)
  {
    for (int i = 0; i < 8; i++)
    {
      if ((i / 4) % 2 == 1)
      {
        skip[i] = true;  // don't bother handling ri[1].
      }
    }
    rw[0] = 1;  // and weight like we're dead-center on the inner cells.
  }

  // now repeat that structure for phi:

  float p0 = phi / (1.0 * phi_spacing) - 0.5 - phimin_roi_low;  // the position 'phi' in  units of phi_spacing, starting from the center of the 0th bin.
  // print_need_cout("prepping to filter low, p0=%f, phi=%d, phi_spacing=%d, phimin_roi_low=%d\n",p0,phi,phi_spacing,phimin_roi_low);

  int p0i = floor(p0);   // the integer portion of the position. -- what center is below our position?
  float p0d = p0 - p0i;  // the decimal portion of the position. -- how far past center are we?
  int pi[4];             // the local phi position of the four l-bin centers to consider
  // print_need_cout("filtering low, p0i=%d, nphi_low=%d\n",p0i,nphi_low);
  // Q: why do I need to filter here?
  // A: Worst case scenario, this produces -1<p0<0.  If that wraps around and is still in the roi, I should include it
  // A: Likewise, if I get p0>nphi_low, then that means it ought to wrap around and be considered against the other limit.
  pi[0] = FilterPhiIndex(p0i, nphi_low);
  pi[1] = FilterPhiIndex(p0i + 1, nphi_low);
  // having filtered, we will always have a positive number.

  // print_need_cout("done filtering low\n");
  float pw[2];      // the weight of that cell, linear in distance from the center of it
  pw[0] = 1 - p0d;  // 1 when we're on it, 0 when we're at the other one.
  pw[1] = p0d;      // 1 when we're on it, 0 when we're at the other one

  // note that since phi wraps around itself, we have the possibility of being outside by being above/below the opposite end of where we thought we were
  // remember these are coordinates wrt the beginning of the roi, and are positive because we filtered.  If that rotation didn't put them back in the roi, then they weren't before, either.
  if (pi[0] >= nphi_roi_low)
  {
    for (int i = 0; i < 8; i++)
    {
      if ((i / 2) % 2 == 0)
      {
        skip[i] = true;  // don't bother handling pi[0].
      }
    }
    pw[1] = 1;  // and weight like we're dead-center on the high cells.
  }
  if (pi[1] >= nphi_roi_low)
  {
    for (int i = 0; i < 8; i++)
    {
      if ((i / 2) % 2 == 1)
      {
        skip[i] = true;  // don't bother handling pi[1].
      }
    }
    pw[0] = 1;  // and weight like we're dead-center on the high cells.
  }

  // and once more for z.  ooph.

  float z0 = z / (1.0 * z_spacing) - 0.5 - zmin_roi_low;  // the position in r, in units of r_spacing, starting from the center of the 0th bin.
  int z0i = floor(z0);                                    // the integer portion of the position. -- what center is below our position?
  float z0d = z0 - z0i;                                   // the decimal portion of the position. -- how far past center are we?
  // instead of listing all eight, I'll address these as i%2, (i/2)%2 and (i/4)%2 to avoid typos
  int zi[2];  // the r position of the eight cell centers to consider.
  zi[0] = z0i;
  zi[1] = z0i + 1;
  float zw[2];      // the weight of that cell, linear in distance from the center of it
  zw[0] = 1 - z0d;  // 1 when we're on it, 0 when we're at the other one.
  zw[1] = z0d;      // 1 when we're on it, 0 when we're at the other one

  if (zi[0] < 0)
  {
    for (int i = 0; i < 8; i++)
    {
      if ((i) % 2 == 0)
      {
        skip[i] = true;  // don't bother handling zi[0].
      }
    }
    zw[1] = 1;  // and weight like we're dead-center on the higher cells.
  }
  if (zi[1] >= nz_roi_low)
  {
    for (int i = 0; i < 8; i++)
    {
      if ((i) % 2 == 1)
      {
        skip[i] = true;  // don't bother handling zi[1].
      }
    }
    zw[0] = 1;  // and weight like we're dead-center on the lower cells.
  }

  TVector3 sum(0, 0, 0);
  // at this point, we should be skipping all destination l-bins that are out-of-bounds.
  // note that if any out-of-bounds ones survive somehow, the call to Epartial_lowres will fail loudly.
  int lBinEdge[2];              // lower and upper (included) edges of the low-res bin, measured in f-bins, reused per-dimension
  int hRegionEdge[2];           // lower and upper edge of the high-res region, measured in f-bins, reused per-dimension.
  bool overlapsPhi, overlapsZ;  // whether we overlap in R, phi, and z.

  int r_highres_dist = (nr_high - 1) / 2;
  int phi_highres_dist = (nphi_high - 1) / 2;
  int z_highres_dist = (nz_high - 1) / 2;

  for (int ir = 0; ir < nr_low; ir++)
  {
    lBinEdge[0] = ir * r_spacing;
    lBinEdge[1] = (ir + 1) * r_spacing - 1;
    hRegionEdge[0] = r - r_highres_dist;
    hRegionEdge[1] = r + r_highres_dist;
    bool overlapsR = (lBinEdge[0] <= hRegionEdge[1] && hRegionEdge[0] <= lBinEdge[1]);
    for (int iphi = 0; iphi < nphi_low; iphi++)
    {
      lBinEdge[0] = iphi * phi_spacing;
      lBinEdge[1] = (iphi + 1) * phi_spacing - 1;
      hRegionEdge[0] = phi - phi_highres_dist;
      hRegionEdge[1] = phi + phi_highres_dist;
      overlapsPhi = (lBinEdge[0] <= hRegionEdge[1] && hRegionEdge[0] <= lBinEdge[1]);
      for (int iz = 0; iz < nz_low; iz++)
      {
        lBinEdge[0] = iz * z_spacing;
        lBinEdge[1] = (iz + 1) * z_spacing - 1;
        hRegionEdge[0] = z - z_highres_dist;
        hRegionEdge[1] = z + z_highres_dist;
        overlapsZ = (lBinEdge[0] <= hRegionEdge[1] && hRegionEdge[0] <= lBinEdge[1]);
        // conceptually: see if the l-bin overlaps with the high-res region:
        // the high-res region includes all indices from r-(nr_high-1)/2 to r+(nr_high-1)/2.
        // each low-res region includes all indices from ir*r_spacing to (ir+1)*r_spacing-1.
        if (overlapsR && overlapsPhi && overlapsZ)
        {
          // if their bounds are interleaved in all dimensions, there is overlap, and we've already summed this region.
          continue;
        }
        else
        {
          // if(debugFlag()) print_need_cout("%d: AnnularFieldSim::sum_field_at, considering l-bin at(r=%d,phi=%d, z=%d)\n",__LINE__,ir,iphi,iz);

          for (int i = 0; i < 8; i++)
          {
            if (skip[i])
            {
              continue;
            }
            if (ri[(i / 4) % 2] + rmin_roi_low == ir && pi[(i / 2) % 2] + phimin_roi_low == iphi && zi[(i) % 2] + zmin_roi_low == iz)
            {
              std::cout << boost::str(boost::format("considering an l-bins effect on itself, r=%d,phi=%d,z=%d (matches i=%d, not skipped), means we're not interpolating fairly") % ir % iphi % iz % i) << std::endl;
              exit(1);
            }
            // the ri, pi, and zi elements are relative to the roi, as needed for Epartial.
            // the ir, iphi, and iz are all absolute, as needed for q_lowres
            if (pi[(i / 2) % 2] < 0)
            {
              std::cout << boost::str(boost::format("%d: Getting with phi=%d") % __LINE__ % (pi[(i / 2) % 2])) << std::endl;
            }
            sum += (Epartial_lowres->Get(ri[(i / 4) % 2], pi[(i / 2) % 2], zi[(i) % 2], ir, iphi, iz) * q_lowres->Get(ir, iphi, iz)) * zw[(i) % 2] * pw[(i / 2) % 2] * rw[(i / 4) % 2];
          }
        }
      }
    }
  }

  return sum;
}

TVector3 AnnularFieldSim::sum_phislice_field_at(int r, int phi, int z)
{
  // sum the E field over all nr by ny by nz cells of sources, at the specific position r,phi,z.
  // note the specific position in Epartial is in relative coordinates.
  // print_need_cout("AnnularFieldSim::sum_field_at(r=%d,phi=%d, z=%d)\n",r,phi,z);
  TVector3 pos = GetRoiCellCenter(r - rmin_roi, phi - phimin_roi, z - zmin_roi);
  TVector3 slicepos = GetRoiCellCenter(r - rmin_roi, 0, z - zmin_roi);
  float rotphi = pos.Phi() - slicepos.Phi();  // probably this is phi*step.Phi();

  // caution:  if you print progress of each of these, it will print a great deal of data, since this is called per-bin
  // unsigned long long totalelements=nr*nphi*nz;
  //  unsigned long long percent=totalelements/100;

  // unsigned long long el=0;

  TVector3 sum(0, 0, 0);
  TVector3 unrotatedField(0, 0, 0);
  TVector3 unitField(0, 0, 0);
  int phirel;
  for (int ir = 0; ir < nr; ir++)
  {
    for (int iphi = 0; iphi < nphi; iphi++)
    {
      for (int iz = 0; iz < nz; iz++)
      {
        // sum+=*partial[x][phi][z][ix][iphi][iz] * *q[ix][iphi][iz];
        if (r == ir && phi == iphi && z == iz)
        {
          continue;  // dont' compute self-to-self field.
        }
        phirel = FilterPhiIndex(iphi - phi);
        unitField = Epartial_phislice->Get(r - rmin_roi, 0, z - zmin_roi, ir, phirel, iz);
        unitField.RotateZ(rotphi);  // previously was rotate by the step.Phi()*phi.    //annoying that I can't rename this to 'rotated field' here without unnecessary overhead.

        sum += unitField * q->GetChargeInBin(ir, iphi, iz);
        ;

        /*
        if(!(el%percent)) {print_need_cout("summing phislices %d%%:  ",(int)(el/percent));
          print_need_cout("unit field at (r=%d,p=%d,z=%d) from  (ir=%d,ip=%d,iz=%d) is (%E,%E,%E) (xyz), q=%E\n",
                 r,phi,z,ir,iphi,iz,unitField.X(),unitField.Y(),unitField.Z(),q->GetChargeInBin(ir,iphi,iz));
        }
        el++;
        */
      }
    }
  }
  // print_need_cout("summed field at (%d,%d,%d)=(%f,%f,%f)\n",x,y,z,sum.X(),sum.Y(),sum.Z());
  return sum;
}

TVector3 AnnularFieldSim::swimToInAnalyticSteps(float zdest, TVector3 start, int steps = 1, int *goodToStep = nullptr)
{
  // assume coordinates are given in native units (cm=1 unless that changed!).
  double zdist = (zdest - start.Z()) * cm;
  double zstep = zdist / steps;
  start = start * cm;  // scale us to cm.

  TVector3 ret = start;
  TVector3 accumulated_distortion(0, 0, 0);
  TVector3 accumulated_drift(0, 0, 0);
  TVector3 drift_step(0, 0, zstep);

  int rt, pt, zt;  // just placeholders for the bounds-checking.
  BoundsCase zBound;
  for (int i = 0; i < steps; i++)
  {
    zBound = GetZindexAndCheckBounds(ret.Z(), &zt);
    if (zBound == OnLowEdge)
    {
      // print_need_cout("AnnularFieldSIm::swimToInAnalyticSteps requests z-nudge from z=%f to %f\n", ret.Z(), ret.Z()+ALMOST_ZERO);//nudge it in z:
      ret.SetZ(ret.Z() + ALMOST_ZERO);
    }
    if (GetRindexAndCheckBounds(ret.Perp(), &rt) != InBounds || GetPhiIndexAndCheckBounds(FilterPhiPos(ret.Phi()), &pt) != InBounds || (zBound == OutOfBounds))
    {
      std::cout << boost::str(boost::format(
                                  "AnnularFieldSim::swimToInAnalyticSteps at step %d asked to swim particle from (%f,%f,%f)(cm) (rphiz)=(%fcm,%frad,%fcm)which is outside the ROI.") %
                              i % (ret.X() / cm) % (ret.Y() / cm) % (ret.Z() / cm) % (ret.Perp() / cm) % ret.Phi() % (ret.Z() / cm))
                << std::endl;
      std::cout << boost::str(boost::format(" -- %f <= r < %f \t%f <= phi < %f \t%f <= z < %f ") % (rmin_roi * step.Perp() + rmin) % (rmax_roi * step.Perp() + rmin) % (phimin_roi * step.Phi()) % (phimax_roi * step.Phi()) % (zmin_roi * step.Z()) % (zmax_roi * step.Z())) << std::endl;
      std::cout << "Returning last valid position." << std::endl;
      ;
      if (!(goodToStep == nullptr))
      {
        *goodToStep = i - 1;
      }
      return ret * (1.0 / cm);
    }
    // print_need_cout("AnnularFieldSim::swimToInAnalyticSteps at step %d, asked to swim particle from (%f,%f,%f) (rphiz)=(%f,%f,%f).\n",i,ret.X(),ret.Y(),ret.Z(),ret.Perp(),ret.Phi(),ret.Z());
    // rcc note: once I put the z distoriton back in, I need to check that ret.Z+zstep is still in bounds:
    accumulated_distortion += GetStepDistortion(ret.Z() + zstep, ret, true, true);

    accumulated_drift += drift_step;

    // this seems redundant, but if the distortions are small they may lose precision and stop actually changing the position when step size is small.  This allows them to accumulate separately so they can grow properly:
    ret = start + accumulated_distortion + accumulated_drift;
  }

  return ret * (1.0 / cm);
}

TVector3 AnnularFieldSim::swimToInSteps(float zdest, const TVector3 &start, int steps = 1, bool interpolate = false, int *goodToStep = nullptr)
{
  int *success = nullptr;
  TVector3 straightline(start.X(), start.Y(), zdest);
  TVector3 distortion = GetTotalDistortion(zdest, start, steps, interpolate, goodToStep, success);
  return straightline + distortion;
}

// NOLINTNEXTLINE(misc-no-recursion)
TVector3 AnnularFieldSim::GetTotalDistortion(float zdest, const TVector3 &start, int steps, bool interpolate, int *goodToStep, int *success)
{
  // work in native units is automatic.
  // double zdist=(zdest-start.Z())*cm;
  // start=start*cm;

  // check the z bounds:
  int rt, pt, zt;  // just placeholders for the bounds-checking.
  BoundsCase zBound;
  zBound = GetZindexAndCheckBounds(zdest, &zt);
  bool rdrswitch = RdeltaRswitch;
  *success = 0;
  if (zBound == OutOfBounds)
  {
    if (hasTwin)
    {  // check if this destination is valid for our twin, if we have one:
      zBound = GetZindexAndCheckBounds(-zdest, &zt);
      if (zBound == InBounds || zBound == OnLowEdge)
      {
        // destination is in the twin
        return twin->GetTotalDistortion(zdest, start, steps, interpolate, goodToStep, success);
      }
    }
    // otherwise, we're not in the twin, and default to our usual gripe:
    std::cout << boost::str(boost::format("AnnularFieldSim::GetTotalDistortion starting at (%f,%f,%f)=(r%f,p%f,z%f) asked to drift to z=%f, which is outside the ROI.  hasTwin= %d.  Returning zero_vector.") % start.X() % start.Y() % start.Z() % start.Perp() % FilterPhiPos(start.Phi()) % start.Z() % zdest % (int) hasTwin) << std::endl;
    std::cout << boost::str(boost::format(" -- %f <= r < %f \t%f <= phi < %f \t%f <= z < %f ") % (rmin_roi * step.Perp() + rmin) % (rmax_roi * step.Perp() + rmin) % (phimin_roi * step.Phi()) % (phimax_roi * step.Phi()) % (zmin_roi * step.Z()) % (zmax_roi * step.Z())) << std::endl;
    *goodToStep = 0;
    *success = 0;
    return zero_vector;  // rcchere
  }
  else if (zBound == OnLowEdge)
  {
    // nudge it in z:
    zdest += ALMOST_ZERO;
  }
  zBound = GetZindexAndCheckBounds(start.Z(), &zt);
  if (zBound == OutOfBounds)
  {
    std::cout << boost::str(boost::format("AnnularFieldSim::GetTotalDistortion starting at (%f,%f,%f)=(r%f,p%f,z%f) asked to drift from z=%f, which is outside the ROI.  Returning zero_vector.") % start.X() % start.Y() % start.Z() % start.Perp() % FilterPhiPos(start.Phi()) % start.Z() % start.Z()) << std::endl;
    std::cout << boost::str(boost::format(" -- %f <= r < %f \t%f <= phi < %f \t%f <= z < %f ") % (rmin_roi * step.Perp() + rmin) % (rmax_roi * step.Perp() + rmin) % (phimin_roi * step.Phi()) % (phimax_roi * step.Phi()) % (zmin_roi * step.Z()) % (zmax_roi * step.Z())) << std::endl;
    *goodToStep = 0;
    *success = 0;
    return zero_vector;
  }
  else if (zBound == OnLowEdge)
  {
    // nudge it in z:
    zdest += ALMOST_ZERO;
  }

  // now we are guaranteed the z limits are in range, and don't need to check them again.

  double zstep = (zmax - zmin) / steps;  // always a positive number
  // print_need_cout("Lets do some handy math here, the value of the first zstep: %f \n",zstep);
  // print_need_cout("Lets do some handy math here, the value of zmax is: %f \n",zmax);
  // print_need_cout("Lets do some handy math here, the value of zmin is: %f \n",zmin);
  // print_need_cout("Lets do some handy math here, the value of steps is: %d \n",steps);
  if ((zdest - start.Z()) < 0)
  {
    zstep = -zstep;  // It goes two directions, so there will be positive and negative value
  }
  int integernumbersteps = (zdest - start.Z()) / zstep;
  // print_need_cout("Lets do some handy math here, the value of integernumbersteps is: %d \n",integernumbersteps);
  // print_need_cout("Lets do some handy math here, the value of zdest is: %f \n",zdest);
  // print_need_cout("Lets do some handy math here, the value of start.Z() is: %f \n",start.Z());
  // print_need_cout("Lets do some handy math here, the value of zstep is: %f \n",zstep);

  // print_need_cout("Hello, this is the zmax-zmin length: %f\n", (zmax-zmin));
  // print_need_cout("Hello, this is the zdest-start.Z(): %f\n", (zdest-start.Z()));
  // print_need_cout("Hello, this is the zstep number: %f\n", zstep);
  // print_need_cout("Hello, this is the integernumbersteps number: %d\n", integernumbersteps);

  TVector3 position = start;
  TVector3 accumulated_distortion(0, 0, 0);
  TVector3 accumulated_drift(0, 0, 0);
  TVector3 drift_step(0, 0, zstep);

  // the conceptual approach here is to get the vector distortion in each z step, and use the transverse component of that to update the position of the particle for the next step, while accumulating the total distortion separately from the position.  This allows a small residual to accumulate, rather than being lost.  We do not correct the z position, so that the stepping does not 'skip' parts of the trajectory.

  // describe what you're doing here:
  // print_need_cout("Hello again, this is the zdest: %f\n", (zdest));
  // print_need_cout("Hello again, this is the start.Z(): %f\n", (start.Z()));

  accumulated_distortion += GetStepDistortion(zdest - (integernumbersteps * zstep), position, true, false);

  // short-circuit if there's no travel length:

  position.SetX(start.X() + accumulated_distortion.X());
  position.SetY(start.Y() + accumulated_distortion.Y());
  position.SetZ(zdest - (integernumbersteps * zstep));

  float StartR, FinalR;
  float DeltaR;
  for (int i = 0; i < integernumbersteps; i++)
  {
    if (GetRindexAndCheckBounds(position.Perp(), &rt) != InBounds || GetPhiIndexAndCheckBounds(FilterPhiPos(position.Phi()), &pt) != InBounds || (zBound == OutOfBounds))
    {
      // print_need_cout("AnnularFieldSim::GetTotalDistortion starting at (%f,%f,%f)=(r%f,p%f,z%f) with drift_step=%f, at step %d, asked to swim particle from (%f,%f,%f) (rphiz)=(%f,%f,%f)which is outside the ROI.\n", start.X(), start.Y(), start.Z(), start.Perp(), FilterPhiPos(start.Phi()), start.Z(), zstep, i, position.X(), position.Y(), position.Z(), position.Perp(), position.Phi(), position.Z());
      // print_need_cout(" -- %f <= r < %f \t%f <= phi < %f \t%f <= z < %f \n", rmin_roi * step.Perp() + rmin, rmax_roi * step.Perp() + rmin, phimin_roi * step.Phi(), phimax_roi * step.Phi(), zmin_roi * step.Z(), zmax_roi * step.Z());
      // print_need_cout("Returning last good position.\n");
      if (!(goodToStep == nullptr))
      {
        *goodToStep = i - 1, *success = 0;
      }
      // assert (1==2);
      return (accumulated_distortion);
    }

    StartR = position.Perp();
    // as we accumulate distortions, add these to the x and y positions, but do not change the z position, otherwise we will 'skip' parts of the drift, which is not the intended behavior.
    accumulated_distortion += GetStepDistortion(position.Z() + zstep, position, true, false);
    position.SetX(start.X() + accumulated_distortion.X());
    position.SetY(start.Y() + accumulated_distortion.Y());
    FinalR = position.Perp();
    // PositionXAfter=start.X() + accumulated_distortion.X();
    // PositionYAfter=start.Y() + accumulated_distortion.Y();//

    // StartR=sqrt(PositionXBefore*PositionXBefore+PositionYBefore*PositionYBefore);
    DeltaR = FinalR - StartR;  // sqrt(PositionXAfter*PositionXAfter+PositionYAfter*PositionYAfter)-StartR;
    if (rdrswitch)
    {
      hRdeltaRComponent->Fill(StartR, DeltaR);
    }
    position.SetZ(position.Z() + zstep);
  }
  if (GetRindexAndCheckBounds(position.Perp(), &rt) != InBounds || GetPhiIndexAndCheckBounds(FilterPhiPos(position.Phi()), &pt) != InBounds || (zBound == OutOfBounds))
  {
    *goodToStep = integernumbersteps - 1;
    *success = 0;
    return (accumulated_distortion);
  }
  *goodToStep = integernumbersteps;
  *success = 1;
  return accumulated_distortion;
}

void AnnularFieldSim::PlotFieldSlices(const std::string &filebase, const TVector3 &pos, char which)
{
  bool mapEfield = true;
  if (which == 'B')
  {
    mapEfield = false;
  }
  which = mapEfield ? 'E' : 'B';
  std::string units;
  if (mapEfield)
  {
    units = "V/cm";
  }
  else
  {
    units = "T";
  }

  std::cout << boost::str(boost::format("plotting field slices for %c field, slicing at (%1.2F,%1.2f,%1.2f)...") % which % pos.Perp() % FilterPhiPos(pos.Phi()) % pos.Z()) << std::endl;
  std::cout << "file=" << filebase << std::endl;
  ;
  std::string plotfilename = filebase + "." + which + "field_slices.pdf";
  TVector3 inner = GetInnerEdge();
  TVector3 outer = GetOuterEdge();
  TVector3 steplocal = GetFieldStep();

  TCanvas *c;

  TH2F *hEfield[3][3];
  TH2F *hCharge[3];
  TH1F *hEfieldComp[3][3];
  std::string axis[] = {"r", "p", "z", "r", "p", "z"};
  float axisval[] = {(float) pos.Perp(), (float) FilterPhiPos(pos.Phi()), (float) pos.Z(), (float) pos.Perp(), (float) FilterPhiPos(pos.Phi()), (float) pos.Z()};
  int axn[] = {nr_roi, nphi_roi, nz_roi, nr_roi, nphi_roi, nz_roi};
  float axtop[] = {(float) outer.Perp(), 2 * M_PI, (float) outer.Z(), (float) outer.Perp(), 2 * M_PI, (float) outer.Z()};
  float axbot[] = {(float) inner.Perp(), 0, (float) inner.Z(), (float) inner.Perp(), 0, (float) inner.Z()};

  // if we are in charge of a twin, extend our axes:
  if (hasTwin)
  {
    // axtop[2]=axtop[5]=(float)(twin->GetOuterEdge().Z());
    axbot[2] = axbot[5] = (float) (twin->GetInnerEdge().Z());
  }
  std::cout << boost::str(boost::format("rpz bounds are %f<r%f\t %f<phi%f\t %f<z%f") % axbot[0] % axtop[0] % axbot[1] % axtop[1] % axbot[2] % axtop[2]) << std::endl;
  float axstep[6];
  for (int i = 0; i < 6; i++)
  {
    axstep[i] = (axtop[i] - axbot[i]) / (1.0 * axn[i]);
  }
  TVector3 field;
  TVector3 lpos;

  // it's a bit meta, but this loop over axes is a compact way to generate all the histogram titles.
  for (int ax = 0; ax < 3; ax++)
  {
    // loop over which axis slice to take
    hCharge[ax] = new TH2F(std::string("hCharge" + axis[ax]).c_str(),
                           boost::str(boost::format("Spacecharge Distribution in the %s%s plane at %s=%2.3f (C/cm^3);%s;%s") % axis[ax + 1] % axis[ax + 2] % axis[ax] % axisval[ax] % axis[ax + 1] % axis[ax + 2]).c_str(),
                           axn[ax + 1], axbot[ax + 1], axtop[ax + 1],
                           axn[ax + 2], axbot[ax + 2], axtop[ax + 2]);
    for (int i = 0; i < 3; i++)
    {
      // loop over which axis of the field to read
      hEfield[ax][i] = new TH2F(boost::str(boost::format("h%sfield%s_%s%s") % which % axis[i] % axis[ax + 1] % axis[ax + 2]).c_str(),
                                boost::str(boost::format("%s component of %s Field in the %s%s plane at %s=%2.3f (%s);%s;%s") % axis[i] % which % axis[ax + 1] % axis[ax + 2] % axis[ax] % axisval[ax] % units.c_str() % axis[ax + 1] % axis[ax + 2]).c_str(),
                                axn[ax + 1], axbot[ax + 1], axtop[ax + 1],
                                axn[ax + 2], axbot[ax + 2], axtop[ax + 2]);
      hEfieldComp[ax][i] = new TH1F(boost::str(boost::format("h%sfieldComp%s_%s%s") % which % axis[i] % axis[ax + 1] % axis[ax + 2]).c_str(),
                                    boost::str(boost::format("Log Magnitude of %s component of %s Field in the %s%s plane at %s=%2.3f;log10(mag)") % axis[i] % which % axis[ax + 1] % axis[ax + 2] % axis[ax] % axisval[ax]).c_str(),
                                    200, -5, 5);
    }
  }

  float rpz_coord[3];
  for (int ax = 0; ax < 3; ax++)
  {  // we have three sets of 'slices'.  the R slice is a 2d plot in phi-z, etc.
    rpz_coord[ax] = axisval[ax] + axstep[ax] / 2;
    for (int i = 0; i < axn[ax + 1]; i++)
    {  // for each slice, loop over the bins of the 2d plot:
      rpz_coord[(ax + 1) % 3] = axbot[ax + 1] + (i + 0.5) * axstep[ax + 1];
      for (int j = 0; j < axn[ax + 2]; j++)
      {
        rpz_coord[(ax + 2) % 3] = axbot[ax + 2] + (j + 0.5) * axstep[ax + 2];
        lpos.SetXYZ(rpz_coord[0], 0, rpz_coord[2]);
        lpos.SetPhi(rpz_coord[1]);
        if (false && ax == 0)
        {
          std::cout << boost::str(boost::format("sampling rpz=(%f,%f,%f)=(%f,%f,%f) after conversion to xyz=(%f,%f,%f)") % rpz_coord[0] % rpz_coord[1] % rpz_coord[2] % lpos.Perp() % (FilterPhiPos(lpos.Phi())) % lpos.Z() % lpos.X() % lpos.Y() % lpos.Z()) << std::endl;
        }
        if (mapEfield)
        {
          // GetFieldAt automatically asks the twin if we are out of bounds here.
          field = GetFieldAt(lpos) * (1.0 * cm / V);  // get units so we're drawing in V/cm when we draw.
        }
        else
        {
          field = GetBFieldAt(lpos) * (1.0 / Tesla);  // get units so we're drawing in Tesla when we draw.
          // if (hasTwin && lpos.Z()<0) field=twin->GetBFieldAt(lpos)*1.0/Tesla;
        }
        field.RotateZ(-rpz_coord[1]);  // rotate us so we can read the y component as the phi component
        // if (field.Mag()>0) continue; //nothing has mag zero because of the drift field.
        hEfield[ax][0]->Fill(rpz_coord[(ax + 1) % 3], rpz_coord[(ax + 2) % 3], field.X());
        hEfield[ax][1]->Fill(rpz_coord[(ax + 1) % 3], rpz_coord[(ax + 2) % 3], field.Y());
        hEfield[ax][2]->Fill(rpz_coord[(ax + 1) % 3], rpz_coord[(ax + 2) % 3], field.Z());
        hCharge[ax]->Fill(rpz_coord[(ax + 1) % 3], rpz_coord[(ax + 2) % 3], GetChargeAt(lpos));
        hEfieldComp[ax][0]->Fill((abs(field.X())));
        hEfieldComp[ax][1]->Fill((abs(field.Y())));
        hEfieldComp[ax][2]->Fill((abs(field.Z())));
      }
    }
  }

  c = new TCanvas("cfieldslices", "electric field", 1200, 800);
  c->Divide(4, 3);
  gStyle->SetOptStat();
  for (int ax = 0; ax < 3; ax++)
  {
    for (int i = 0; i < 3; i++)
    {
      c->cd(ax * 4 + i + 1);
      gPad->SetRightMargin(0.15);
      hEfield[ax][i]->SetStats(false);
      hEfield[ax][i]->Draw("colz");
      // hEfield[ax][i]->GetListOfFunctions()->Print();
      gPad->Modified();
      // hEfieldComp[ax][i]->Draw();//"colz");
    }
    c->cd(ax * 4 + 4);
    gPad->SetRightMargin(0.15);
    hCharge[ax]->SetStats(false);
    hCharge[ax]->Draw("colz");
    // pal=dynamic_cast<TPaletteAxis*>(hCharge[ax]->GetListOfFunctions()->FindObject("palette"));
    if (false)
    {
      // pal->SetX1NDC(0.86);
      // pal->SetX2NDC(0.91);
      gPad->Modified();
    }
  }
  c->SaveAs(plotfilename.c_str());
  std::cout << "after plotting field slices..." << std::endl;
  std::cout << "file=" << filebase << std::endl;

  return;
}

void AnnularFieldSim::GenerateSeparateDistortionMaps(const std::string &filebase, int nSteps, int r_subsamples, int p_subsamples, int z_subsamples, int /*z_substeps*/, bool andCartesian)
{
  // generates the distortion map for one or both sides of the detector, separating them so
  // they do not try to interpolate across the CM.

  // 1) pick a map spacing ('s')
  TVector3 steplocal(step.Perp() / r_subsamples, 0, step.Z() / z_subsamples);
  steplocal.SetPhi(step.Phi() / p_subsamples);
  float deltar = steplocal.Perp();  //(rf-ri)/nr;
  float deltap = steplocal.Phi();   //(pf-pi)/np;
  float deltaz = steplocal.Z();     //(zf-zi)/nz;
  bool rdrswitch = RdeltaRswitch;
  TVector3 stepzvec(0, 0, deltaz);

  //  int nSteps = 500;  //how many steps to take in the particle path.  hardcoded for now.  Think about this later.

  int nSides = 1;
  if (hasTwin)
  {
    nSides = 2;
  }
  // idea for a faster way to build a map:

  // 2) generate the distortions steplocal.Z() away from the readout
  // 3) generate the distortion from (i)*steplocal.Z() away to (i-1) away, then add the interpolated value from the (i-1) layer

  // for interpolation, Henry needs one extra buffer bin on each side.

  // so we define the histogram bounds (the 'h' suffix) to be the full range
  // plus an additional step in each direction so interpolation can work at the edges
  TVector3 lowerEdge = GetRoiCellCenter(rmin_roi, phimin_roi, zmin_roi);
  TVector3 upperEdge = GetRoiCellCenter(rmax_roi - 1, phimax_roi - 1, zmax_roi - 1);
  int nph = nphi * p_subsamples + 2;  // number of phibins in the histogram
  int nrh = nr * r_subsamples + 2;    // number of r bins in the histogram
  int nzh = nz * z_subsamples + 2;    // number of z you get the idea.

  float rih = lowerEdge.Perp() - 0.5 * step.Perp() - steplocal.Perp();             // lower bound of roi, minus one
  float rfh = upperEdge.Perp() + 0.5 * step.Perp() + steplocal.Perp();             // upper bound of roi, plus one
  float pih = FilterPhiPos(lowerEdge.Phi()) - 0.5 * step.Phi() - steplocal.Phi();  // can't automate this or we'll run afoul of phi looping.
  float pfh = FilterPhiPos(upperEdge.Phi()) + 0.5 * step.Phi() + steplocal.Phi();  // can't automate this or we'll run afoul of phi looping.
  float zih = lowerEdge.Z() - 0.5 * step.Z() - steplocal.Z();                      // lower bound of roi, minus one
  float zfh = upperEdge.Z() + 0.5 * step.Z() + steplocal.Z();                      // upper bound of roi, plus one
  float z_readout = upperEdge.Z() + 0.5 * step.Z();                                // readout plane.  Note we assume this is positive.

  std::cout << "generating distortion map..." << std::endl;
  std::cout << "file=" << filebase << std::endl;
  std::cout << boost::str(boost::format("Phi:  %d steps from %f to %f (field has %d steps)") % nph % pih % pfh % GetFieldStepsPhi()) << std::endl;
  std::cout << boost::str(boost::format("R:  %d steps from %f to %f (field has %d steps)") % nrh % rih % rfh % GetFieldStepsR()) << std::endl;
  std::cout << boost::str(boost::format("Z:  %d steps from %f to %f (field has %d steps)") % nzh % zih % zfh % GetFieldStepsZ()) << std::endl;
  std::string distortionFilename = filebase + ".distortion_map.hist.root";
  std::string summaryFilename = filebase + "distortion_summary.pdf";
  std::string diffSummaryFilename = filebase + "differential_summary.pdf";

  TFile *outf = TFile::Open(distortionFilename.c_str(), "RECREATE");
  outf->cd();

  // actual output maps:

  TH3F *hDistortionR = new TH3F("hDistortionR", "Per-z-bin Distortion in the R direction as a function of (phi,r,z) (centered in r,phi, z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);
  TH3F *hDistortionP = new TH3F("hDistortionP", "Per-z-bin Distortion in the RPhi direction as a function of (phi,r,z)  (centered in r,phi, z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);
  TH3F *hDistortionZ = new TH3F("hDistortionZ", "Per-z-bin Distortion in the Z direction as a function of (phi,r,z)  (centered in r,phi, z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);
  TH3F *hIntDistortionR = new TH3F("hIntDistortionR", "Integrated R Distortion from (phi,r,z) to z=0 (centered in r,phi, and z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);
  TH3F *hIntDistortionP = new TH3F("hIntDistortionP", "Integrated R Distortion from (phi,r,z) to z=0 (centered in r,phi, and z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);
  TH3F *hIntDistortionZ = new TH3F("hIntDistortionZ", "Integrated R Distortion from (phi,r,z) to z=0  (centered in r,phi, and z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);

  TH3F *hIntDistortionX = new TH3F("hIntDistortionX", "Integrated X Distortion from (phi,r,z) to z=0 (centered in r,phi, and z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);
  TH3F *hIntDistortionY = new TH3F("hIntDistortionY", "Integrated Y Distortion from (phi,r,z) to z=0 (centered in r,phi, and z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);

  const int nMapComponents = 6;
  TH3F *hSeparatedMapComponent[2][6];  // side, then xyzrp
  TH3F *hValidToStepComponent[2];
  TH3F *hReachesReadout[2];
  std::string side[2];
  side[0] = "soloz";
  if (hasTwin)
  {
    side[0] = "posz";
    side[1] = "negz";
  }
  std::string sepAxis[] = {"X", "Y", "Z", "R", "P", "RPhi"};
  float zlower, zupper;
  for (int i = 0; i < nSides; i++)
  {
    if (i == 0)
    {  // doing the positive side
      zlower = fmin(zih, zfh);
      zupper = fmax(zih, zfh);
    }
    else
    {
      zlower = -1 * fmax(zih, zfh);
      zupper = -1 * fmin(zih, zfh);
    }
    for (int j = 0; j < nMapComponents; j++)
    {
      hSeparatedMapComponent[i][j] = new TH3F(std::string("hIntDistortion" + sepAxis[j] + "_" + side[i]).c_str(),
                                              std::string("Integrated" + sepAxis[j] + " Deflection drifting from (phi,r,z) to z=endcap;phi;r;z (" + side[i] + " side)").c_str(),
                                              nph, pih, pfh, nrh, rih, rfh, nzh, zlower, zupper);
    }
    hValidToStepComponent[i] = new TH3F(std::string("hValidToStep_" + side[i]).c_str(),
                                        std::string("Step before out of bounds drifting from (phi,r,z) to z=endcap;phi;r;z (" + side[i] + " side)").c_str(),
                                        nph, pih, pfh, nrh, rih, rfh, nzh, zlower, zupper);
    hReachesReadout[i] = new TH3F(std::string("hReachesReadout_" + side[i]).c_str(),
                                  std::string("Bool value for checking if successfuly drifting from (phi,r,z) to z=endcap and reaches the readout;phi;r;z (" + side[i] + " side)").c_str(),
                                  nph, pih, pfh, nrh, rih, rfh, nzh, zlower, zupper);
  }
  if (rdrswitch)
  {
    hRdeltaRComponent = new TH2F("hRdeltaR", "Bool value of input for a 2D R_versus_deltaR plot from Zstart to Zend;R;deltaR", nrh, rih, rfh, 10 * nrh, rih - 70, rfh + 170);
    // hRdeltaRComponent = new TH2F("hRdeltaR","Bool value of input for a 2D R_versus_deltaR plot from Zstart to Zend;R;deltaR", nrh, rih, rfh, 10*nrh, rih-20, rfh-77);
    twin->hRdeltaRComponent = hRdeltaRComponent;
  }
  // monitor plots, and the position that that plot monitors at:

  // TVector3 pos((nrh/2+0.5)*s.Perp()+rih,0,(nzh/2+0.5)*s.Z()+zih);
  TVector3 pos(((int) (nrh / 2) + 0.5) * steplocal.Perp() + rih, 0, zmin + ((int) (nz / 2) + 0.5) * step.Z());
  float posphi = ((int) (nph / 2) + 0.5) * steplocal.Phi() + pih;
  pos.SetPhi(posphi);
  // int xi[3]={nrh/2,nph/2,nzh/2};
  int xi[3] = {(int) floor((pos.Perp() - rih) / steplocal.Perp()), (int) floor((posphi - pih) / steplocal.Phi()), (int) floor((pos.Z() - zih) / steplocal.Z())};
  if (!hasTwin)
  {
    std::cout << boost::str(boost::format("rpz slice indices= (%d,%d,%d) (no twin)") % xi[0] % xi[1] % xi[2]) << std::endl;
  }
  int twinz = (-pos.Z() - zih) / steplocal.Z();
  if (hasTwin)
  {
    std::cout << boost::str(boost::format("rpz slice indices= (%d,%d,%d) twinz=%d") % xi[0] % xi[1] % xi[2] % twinz) << std::endl;
  }

  const std::string axname[] = {"r", "p", "z", "r", "p", "z"};
  int axn[] = {nrh, nph, nzh, nrh, nph, nzh};
  float axval[] = {(float) pos.Perp(), (float) pos.Phi(), (float) pos.Z(), (float) pos.Perp(), (float) pos.Phi(), (float) pos.Z()};
  float axbot[] = {rih, pih, zih, rih, pih, zih};
  float axtop[] = {rfh, pfh, zfh, rfh, pfh, zfh};
  TH2F *hIntDist[3][3];
  TH1F *hRDist[2][3];  // now with a paired friend for the negative side
  TH2F *hDiffDist[3][3];
  TH1F *hRDiffDist[2][3];
  for (int i = 0; i < 3; i++)
  {
    // loop over which axis of the distortion to read
    for (int ax = 0; ax < 3; ax++)
    {
      // loop over which plane to work in
      hDiffDist[ax][i] = new TH2F(std::string("hDiffDist" + axname[i] + "_" + axname[ax + 1] + axname[ax + 2]).c_str(),
                                  boost::str(boost::format("%s component of diff. distortion in  %s%s plane at %s=%2.3f;%s;%s") % axname[i] % axname[ax + 1] % axname[ax + 2] % axname[ax] % axval[ax] % axname[ax + 1] % axname[ax + 2]).c_str(),
                                  axn[ax + 1], axbot[ax + 1], axtop[ax + 1],
                                  axn[ax + 2], axbot[ax + 2], axtop[ax + 2]);
      hIntDist[ax][i] = new TH2F(std::string("hIntDist" + axname[i] + "_" + axname[ax + 1] + axname[ax + 2]).c_str(),
                                 boost::str(boost::format("%s component of int. distortion in  %s%s plane at %s=%2.3f;%s;%s") % axname[i] % axname[ax + 1] % axname[ax + 2] % axname[ax] % axval[ax] % axname[ax + 1] % axname[ax + 2]).c_str(),
                                 axn[ax + 1], axbot[ax + 1], axtop[ax + 1],
                                 axn[ax + 2], axbot[ax + 2], axtop[ax + 2]);
    }
    hRDist[0][i] = new TH1F(std::string("hRDist" + axname[i]).c_str(),
                            boost::str(boost::format("%s component of int. distortion vs r with %s=%2.3f and %s=%2.3f;r(cm);#delta (cm)") % axname[i] % axname[1] % axval[1] % axname[2] % axval[2]).c_str(),
                            axn[0], axbot[0], axtop[0]);
    hRDist[1][i] = new TH1F(std::string("hRDist" + axname[i] + "Neg").c_str(),
                            boost::str(boost::format("%s component of int. distortion vs r with %s=%2.3f and %s=-%2.3f;r(cm);#delta (cm)") % axname[i] % axname[1] % axval[1] % axname[2] % axval[2]).c_str(),
                            axn[0], axbot[0], axtop[0]);
    hRDiffDist[0][i] = new TH1F(std::string("hRDiffDist" + axname[i]).c_str(),
                                boost::str(boost::format("%s component of diff. distortion vs r with %s=%2.3f and %s=%2.3f;r(cm);#delta (cm)") % axname[i] % axname[1] % axval[1] % axname[2] % axval[2]).c_str(),
                                axn[0], axbot[0], axtop[0]);
    hRDiffDist[1][i] = new TH1F(std::string("hRDiffDist" + axname[i] + "Neg").c_str(),
                                boost::str(boost::format("%s component of diff. distortion vs r with %s=%2.3f and %s=-%2.3f;r(cm);#delta (cm)") % axname[i] % axname[1] % axval[1] % axname[2] % axval[2]).c_str(),
                                axn[0], axbot[0], axtop[0]);
  }

  TVector3 inpart, outpart;
  TVector3 diffdistort, distort;
  int validToStep;
  int successCheck;

  // TTree version:
  float partR, partP, partZ;
  int ir, ip, iz;
  float distortR, distortP, distortZ;
  float distortX, distortY;
  float diffdistR, diffdistP, diffdistZ;
  TTree *dTree = new TTree("dTree", "Distortion per step z");
  dTree->Branch("r", &partR);
  dTree->Branch("p", &partP);
  dTree->Branch("z", &partZ);
  dTree->Branch("ir", &ir);
  dTree->Branch("ip", &ip);
  dTree->Branch("iz", &iz);
  dTree->Branch("dr", &distortR);
  dTree->Branch("drp", &distortP);
  dTree->Branch("dz", &distortZ);

  std::cout << boost::str(boost::format("generating separated distortion map with (%dx%dx%d) grid ") % nrh % nph % nzh) << std::endl;
  unsigned long long totalelements = nrh;
  totalelements *= nph;
  totalelements *= nzh;  // breaking up this multiplication prevents a 32bit math overflow
  if (hasTwin)
  {
    totalelements *= 2;  // if we have a twin, we have twice as many z bins as we thought.
  }

  unsigned long long percent = totalelements / 100;
  unsigned long long waypoint = percent * debug_npercent;
  std::cout << boost::str(boost::format("total elements = %llu") % totalelements) << std::endl;

  int el = 0;

  // we want to loop over the entire region to be mapped, but we also need to include
  // one additional bin at each edge, to allow the mc drift code to interpolate properly.
  // hence we count from -1 to n+1, and manually adjust the position in those edge cases
  // to avoid sampling nonphysical regions in r and z.  the phi case is free to wrap as
  //  normal.

  // note that we apply the adjustment to the particle position (inpart) and not the plotted position (partR etc)
  inpart.SetXYZ(1, 0, 0);
  for (ir = 0; ir < nrh; ir++)
  {
    partR = (ir + 0.5) * deltar + rih;
    if (ir == 0)
    {
      inpart.SetPerp(partR + deltar);
    }
    else if (ir == nrh - 1)
    {
      inpart.SetPerp(partR - deltar);
    }
    else
    {
      inpart.SetPerp(partR);
    }
    for (ip = 0; ip < nph; ip++)
    {
      partP = (ip + 0.5) * deltap + pih;
      inpart.SetPhi(partP);
      // since phi loops, there's no need to adjust phis that are out of bounds.
      for (iz = 0; iz < nzh; iz++)
      {
        partZ = (iz) *deltaz + zih;  // start us at the EDGE of the bin,
        if (iz == 0)
        {
          inpart.SetZ(partZ + deltaz);
        }
        else if (iz == nzh - 1)
        {
          inpart.SetZ(partZ - deltaz);
        }
        else
        {
          inpart.SetZ(partZ);
        }
        partZ += 0.5 * deltaz;  // move to center of histogram bin.
        for (int localside = 0; localside < nSides; localside++)
        {
          if (localside == 0)
          {
            diffdistort = zero_vector;  // GetTotalDistortion(inpart.Z() + deltaz, inpart, nSteps, true, &validToStep, &successCheck);
            distort = GetTotalDistortion(z_readout, inpart, nSteps, true, &validToStep, &successCheck);
          }
          else
          {
            // if we have more than one side,
            // flip z coords and do the twin instead:
            partZ *= -1;                   // position to place in histogram
            inpart.SetZ(-1 * inpart.Z());  // position to seek in sim
            diffdistort = zero_vector;     // twin->GetTotalDistortion(inpart.Z() - deltaz, inpart, nSteps, true, &validToStep, &successCheck);
            distort = twin->GetTotalDistortion(-z_readout, inpart, nSteps, true, &validToStep, &successCheck);
          }

          diffdistort.RotateZ(-inpart.Phi());  // rotate so that distortion components are wrt the x axis
          diffdistP = diffdistort.Y();         // the phi component is now the y component.
          diffdistR = diffdistort.X();         // and the r component is the x component
          diffdistZ = diffdistort.Z();

          distortX = distort.X();
          distortY = distort.Y();
          distort.RotateZ(-inpart.Phi());  // rotate so that distortion components are wrt the x axis
          distortP = distort.Y();          // the phi component is now the y component.
          distortR = distort.X();          // and the r component is the x component
          distortZ = distort.Z();

          float distComp[nMapComponents];  // by components
          distComp[0] = distortX;
          distComp[1] = distortY;
          distComp[2] = distortZ;
          distComp[3] = distortR;
          distComp[4] = distortP / partR;  // 'P' now refers to phi in radians, rather than unitful
          distComp[5] = distortP;          // 'RPhi' is the one that correponds to the []meter unitful phi-hat value

          for (int c = 0; c < nMapComponents; c++)
          {
            hSeparatedMapComponent[localside][c]->Fill(partP, partR, partZ, distComp[c]);
          }

          hValidToStepComponent[localside]->Fill(partP, partR, partZ, validToStep);
          hReachesReadout[localside]->Fill(partP, partR, partZ, successCheck);

          // recursive integral distortion:
          // get others working first!

          // print_need_cout("part=(rpz)(%f,%f,%f),distortP=%f\n",partP,partR,partZ,distortP);
          hIntDistortionR->Fill(partP, partR, partZ, distortR);
          hIntDistortionP->Fill(partP, partR, partZ, distortP);
          hIntDistortionZ->Fill(partP, partR, partZ, distortZ);

          if (andCartesian)
          {
            hIntDistortionX->Fill(partP, partR, partZ, distortX);
            hIntDistortionY->Fill(partP, partR, partZ, distortY);
          }

          // now we fill particular slices for integral visualizations:
          if (ir == xi[0] && localside == 0)
          {  // r slice
            // print_need_cout("ir=%d, r=%f (pz)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",ir,partR,ip,iz,distortR,distortP);
            hIntDist[0][0]->Fill(partP, partZ, distortR);
            hIntDist[0][1]->Fill(partP, partZ, distortP);
            hIntDist[0][2]->Fill(partP, partZ, distortZ);
            hDiffDist[0][0]->Fill(partP, partZ, diffdistR);
            hDiffDist[0][1]->Fill(partP, partZ, diffdistP);
            hDiffDist[0][2]->Fill(partP, partZ, diffdistZ);
          }
          if (ip == xi[1] && localside == 0)
          {  // phi slice
            // print_need_cout("ip=%d, p=%f (rz)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",ip,partP,ir,iz,distortR,distortP);
            hIntDist[1][0]->Fill(partZ, partR, distortR);
            hIntDist[1][1]->Fill(partZ, partR, distortP);
            hIntDist[1][2]->Fill(partZ, partR, distortZ);
            hDiffDist[1][0]->Fill(partZ, partR, diffdistR);
            hDiffDist[1][1]->Fill(partZ, partR, diffdistP);
            hDiffDist[1][2]->Fill(partZ, partR, diffdistZ);

            if (iz == xi[2] && localside == 0)
            {  // z slices of phi slices= r line at mid phi, mid z:
              hRDist[0][0]->Fill(partR, distortR);
              hRDist[0][1]->Fill(partR, distortP);
              hRDist[0][2]->Fill(partR, distortZ);
              hRDiffDist[0][0]->Fill(partR, diffdistR);
              hRDiffDist[0][1]->Fill(partR, diffdistP);
              hRDiffDist[0][2]->Fill(partR, diffdistZ);
            }
            if (hasTwin && iz == twinz && localside == 1)
            {  // z slices of phi slices= r line at mid phi, mid z:
              hRDist[1][0]->Fill(partR, distortR);
              hRDist[1][1]->Fill(partR, distortP);
              hRDist[1][2]->Fill(partR, distortZ);
              hRDiffDist[1][0]->Fill(partR, diffdistR);
              hRDiffDist[1][1]->Fill(partR, diffdistP);
              hRDiffDist[1][2]->Fill(partR, diffdistZ);
            }
          }
          if (iz == xi[2] && localside == 0)
          {  // z slice
            // print_need_cout("iz=%d, z=%f (rp)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",iz,partZ,ir,ip,distortR,distortP);

            hIntDist[2][0]->Fill(partR, partP, distortR);
            hIntDist[2][1]->Fill(partR, partP, distortP);
            hIntDist[2][2]->Fill(partR, partP, distortZ);
            hDiffDist[2][0]->Fill(partR, partP, diffdistR);
            hDiffDist[2][1]->Fill(partR, partP, diffdistP);
            hDiffDist[2][2]->Fill(partR, partP, diffdistZ);
          }

          if (!(el % waypoint))
          {
            std::cout << boost::str(boost::format("generating distortions %d%%:  ") % ((int) (el / percent)));
            std::cout << boost::str(boost::format("distortion at (ir=%d,ip=%d,iz=%d) is (%E,%E,%E)") % ir % ip % iz % distortR % distortP % distortZ) << std::endl;
          }
          el++;
        }
      }
    }
  }
  std::cout << "Completed distortion generation.  Saving outputs..." << std::endl;

  TCanvas *canvas = new TCanvas("cdistort", "distortion integrals", 1200, 800);
  // take 10 of the bottom of this for data?
  std::cout << "was able to make a tcanvas" << std::endl;
  canvas->cd();
  TPad *c = new TPad("cplots", "distortion integral plots", 0, 0.2, 1, 1);
  canvas->cd();
  TPad *textpad = new TPad("ctext", "distortion integral plots", 0, 0.0, 1, 0.2);
  std::cout << "was able to make some tpads" << std::endl;

  c->Divide(4, 3);
  gStyle->SetOptStat();
  std::cout << "was able to interact with gStyle" << std::endl;

  for (int i = 0; i < 3; i++)
  {
    // component
    for (int ax = 0; ax < 3; ax++)
    {
      std::cout << boost::str(boost::format("looping over components i=%d ax=%d") % i % ax) << std::endl;

      // plane
      c->cd(i * 4 + ax + 1);
      gPad->SetRightMargin(0.15);
      hIntDist[ax][i]->SetStats(false);
      hIntDist[ax][i]->Draw("colz");
    }

    std::cout << "drawing R profile " << i << std::endl;

    c->cd(i * 4 + 4);
    hRDist[0][i]->SetStats(false);
    hRDist[0][i]->SetFillColor(kRed);
    hRDist[0][i]->Draw("hist");
    if (hasTwin)
    {
      std::cout << "drawing R profile twin " << i << std::endl;
      hRDist[1][i]->SetStats(false);
      hRDist[1][i]->SetLineColor(kBlue);
      hRDist[1][i]->Draw("hist,same");
    }
  }
  std::cout << "switching to textpad" << std::endl;

  textpad->cd();
  float texpos = 0.9;
  float texshift = 0.12;
  TLatex *tex = new TLatex(0.0, texpos, "Fill Me In");
  std::cout << "built TLatex" << std::endl;

  tex->SetTextSize(texshift * 0.8);
  tex->DrawLatex(0.05, texpos, GetFieldString().c_str());
  texpos -= texshift;
  tex->DrawLatex(0.05, texpos, GetChargeString().c_str());
  texpos -= texshift;
  // tex->DrawLatex(0.05,texpos,Form("Drift Field = %2.2f V/cm",GetNominalE()));texpos-=texshift;
  tex->DrawLatex(0.05, texpos, boost::str(boost::format("Drifting grid of (rp)=(%d x %d) electrons with %d steps") % nrh % nph % nSteps).c_str());
  texpos -= texshift;
  tex->DrawLatex(0.05, texpos, GetLookupString().c_str());
  texpos -= texshift;
  tex->DrawLatex(0.05, texpos, GetGasString().c_str());
  texpos -= texshift;
  if (debug_distortionScale.Mag() > 0)
  {
    tex->DrawLatex(0.05, texpos, boost::str(boost::format("Distortion scaled by (r,p,z)=(%2.2f,%2.2f,%2.2f)") % debug_distortionScale.X() % debug_distortionScale.Y() % debug_distortionScale.Z()).c_str());
    texpos -= texshift;
  }
  texpos = 0.9;
  std::cout << "cd'ing to canvas:" << std::endl;

  canvas->cd();
  std::cout << "draw1" << std::endl;

  c->Draw();
  canvas->cd();
  std::cout << "draw2" << std::endl;
  textpad->Draw();
  std::cout << "was able to complete drawing on both pads" << std::endl;

  canvas->SaveAs(summaryFilename.c_str());

  // canvas->cd();
  // already done TPad *c=new TPad("cplots","distortion differential plots",0,0.2,1,1);
  // canvas->cd();
  // already done TPad *textpad=new TPad("ctext","distortion differential plots",0,0.0,1,0.2);
  // already done c->Divide(4,3);
  // gStyle->SetOptStat();
  for (int i = 0; i < 3; i++)
  {
    // component
    for (int ax = 0; ax < 3; ax++)
    {
      // plane
      c->cd(i * 4 + ax + 1);
      gPad->SetRightMargin(0.15);
      hDiffDist[ax][i]->SetStats(false);
      hDiffDist[ax][i]->Draw("colz");
    }
    c->cd(i * 4 + 4);
    hRDiffDist[0][i]->SetStats(false);
    hRDiffDist[0][i]->SetFillColor(kRed);
    hRDiffDist[0][i]->Draw("hist");
    if (hasTwin)
    {
      hRDiffDist[1][i]->SetStats(false);
      hRDiffDist[1][i]->SetLineColor(kBlue);
      hRDiffDist[1][i]->Draw("hist same");
    }
  }
  textpad->cd();
  texpos = 0.9;
  texshift = 0.12;
  tex->SetTextSize(texshift * 0.8);
  tex->DrawLatex(0.05, texpos, GetFieldString().c_str());
  texpos -= texshift;
  tex->DrawLatex(0.05, texpos, GetChargeString().c_str());
  texpos -= texshift;
  // tex->DrawLatex(0.05,texpos,Form("Drift Field = %2.2f V/cm",GetNominalE()));texpos-=texshift;
  tex->DrawLatex(0.05, texpos, boost::str(boost::format("Drifting grid of (rp)=(%d x %d) electrons with %d steps") % nrh % nph % nSteps).c_str());
  texpos -= texshift;
  tex->DrawLatex(0.05, texpos, GetLookupString().c_str());
  texpos -= texshift;
  tex->DrawLatex(0.05, texpos, GetGasString().c_str());
  texpos -= texshift;
  tex->DrawLatex(0.05, texpos, "Differential Plots");
  texpos -= texshift;
  if (debug_distortionScale.Mag() > 0)
  {
    tex->DrawLatex(0.05, texpos, boost::str(boost::format("Distortion scaled by (r,p,z)=(%2.2f,%2.2f,%2.2f)") % debug_distortionScale.X() % debug_distortionScale.Y() % debug_distortionScale.Z()).c_str());
    texpos -= texshift;
  }
  texpos = 0.9;

  canvas->cd();
  c->Draw();
  canvas->cd();
  textpad->Draw();
  canvas->SaveAs(diffSummaryFilename.c_str());

  std::cout << "saving map histograms to:" << distortionFilename << std::endl;

  outf->cd();
  for (int i = 0; i < nSides; i++)
  {
    std::cout << "Saving side " << side[i] << std::endl;
    for (int j = 0; j < nMapComponents; j++)
    {
      hSeparatedMapComponent[i][j]->GetSumw2()->Set(0);
      hSeparatedMapComponent[i][j]->Write();
    }
    hValidToStepComponent[i]->GetSumw2()->Set(0);
    hValidToStepComponent[i]->Write();

    hReachesReadout[i]->GetSumw2()->Set(0);
    hReachesReadout[i]->Write();
  }
  if (rdrswitch)
  {
    hRdeltaRComponent->GetSumw2()->Set(0);
    hRdeltaRComponent->Write();
  }

  hDistortionR->GetSumw2()->Set(0);
  hDistortionP->GetSumw2()->Set(0);
  hDistortionZ->GetSumw2()->Set(0);
  hIntDistortionR->GetSumw2()->Set(0);
  hIntDistortionP->GetSumw2()->Set(0);
  hIntDistortionZ->GetSumw2()->Set(0);
  if (andCartesian)
  {
    hIntDistortionX->GetSumw2()->Set(0);
    hIntDistortionY->GetSumw2()->Set(0);
  }
  hDistortionR->Write();
  hDistortionP->Write();
  hDistortionZ->Write();
  hIntDistortionR->Write();
  hIntDistortionP->Write();
  hIntDistortionZ->Write();
  if (false && andCartesian)
  {
    hIntDistortionX->Write();
    hIntDistortionY->Write();
  }
  std::cout << "finished writing histograms" << std::endl;
  // dTree->Write();
  //  print_need_cout("wrote dTree\n");
  outf->Close();
  // print_need_cout("map:%s.closed\n",distortionFilename.Data());

  std::cout << "wrote separated map and summary to " << filebase << std::endl;
  std::cout << "map: " << distortionFilename << std::endl;
  std::cout << "summary: " << summaryFilename << std::endl;
  // print_need_cout("wrote map and summary to %s and %s.\n",distortionFilename.Data(),summaryFilename.Data());
  return;
}

void AnnularFieldSim::GenerateDistortionMaps(const std::string &filebase, int r_subsamples, int p_subsamples, int z_subsamples, int /*z_substeps*/, bool andCartesian)
{
  // generates the distortion map for the full detector instead of one map per side.
  // This produces wrong behavior when interpolating across the CM and should not be used
  // unless you're doing some debugging/backward compatibility checking.  Eventually this should be removed  (said in early 2022)
  std::cout << "WARNING:  You called the version of the distortion generator that generates a unified map.  Are you sure you meant to do this?" << std::endl;

  // 1) pick a map spacing ('s')
  bool makeUnifiedMap = true;  // if true, generate a single map across the whole TPC.  if false, save two maps, one for each side.

  TVector3 steplocal(step.Perp() / r_subsamples, 0, step.Z() / z_subsamples);
  steplocal.SetPhi(step.Phi() / p_subsamples);
  float deltar = steplocal.Perp();  //(rf-ri)/nr;
  float deltap = steplocal.Phi();   //(pf-pi)/np;
  float deltaz = steplocal.Z();     //(zf-zi)/nz;
  TVector3 stepzvec(0, 0, deltaz);
  int nSteps = 500;  // how many steps to take in the particle path.  hardcoded for now.  Think about this later.

  // idea for a faster way to build a map:

  // 2) generate the distortions steplocal.Z() away from the readout
  // 3) generate the distortion from (i)*steplocal.Z() away to (i-1) away, then add the interpolated value from the (i-1) layer

  // for interpolation, Henry needs one extra buffer bin on each side.

  // so we define the histogram bounds (the 'h' suffix) to be the full range
  // plus an additional step in each direction so interpolation can work at the edges
  TVector3 lowerEdge = GetRoiCellCenter(rmin_roi, phimin_roi, zmin_roi);
  TVector3 upperEdge = GetRoiCellCenter(rmax_roi - 1, phimax_roi - 1, zmax_roi - 1);
  int nph = nphi * p_subsamples + 2;  // nuber of phibins in the histogram
  int nrh = nr * r_subsamples + 2;    // number of r bins in the histogram
  int nzh = nz * z_subsamples + 2;    // number of z you get the idea.
  int successCheck;

  if (hasTwin && makeUnifiedMap)
  {  // double the z range if we have a twin.  r and phi are the same, unless we had a phi roi...
    lowerEdge.SetZ(-1 * upperEdge.Z());
    nzh += nz * z_subsamples;
  }

  float rih = lowerEdge.Perp() - 0.5 * step.Perp() - steplocal.Perp();             // lower bound of roi, minus one
  float rfh = upperEdge.Perp() + 0.5 * step.Perp() + steplocal.Perp();             // upper bound of roi, plus one
  float pih = FilterPhiPos(lowerEdge.Phi()) - 0.5 * step.Phi() - steplocal.Phi();  // can't automate this or we'll run afoul of phi looping.
  float pfh = FilterPhiPos(upperEdge.Phi()) + 0.5 * step.Phi() + steplocal.Phi();  // can't automate this or we'll run afoul of phi looping.
  float zih = lowerEdge.Z() - 0.5 * step.Z() - steplocal.Z();                      // lower bound of roi, minus one
  float zfh = upperEdge.Z() + 0.5 * step.Z() + steplocal.Z();                      // upper bound of roi, plus one
  float z_readout = upperEdge.Z() + 0.5 * step.Z();                                // readout plane.  Note we assume this is positive.

  std::cout << "generating distortion map..." << std::endl;
  std::cout << "file=" << filebase << std::endl;
  std::cout << boost::str(boost::format("Phi:  %d steps from %f to %f (field has %d steps)") % nph % pih % pfh % GetFieldStepsPhi()) << std::endl;
  std::cout << boost::str(boost::format("R:  %d steps from %f to %f (field has %d steps)") % nrh % rih % rfh % GetFieldStepsR()) << std::endl;
  std::cout << boost::str(boost::format("Z:  %d steps from %f to %f (field has %d steps)") % nzh % zih % zfh % GetFieldStepsZ()) << std::endl;
  std::string distortionFilename = filebase + ".distortion_summary.pdf";
  std::string summaryFilename = filebase + ".distortion_summary.pdf";
  std::string diffSummaryFilename = filebase + ".differential_summary.pdf";

  TFile *outf = TFile::Open(distortionFilename.c_str(), "RECREATE");
  outf->cd();

  // actual output maps:

  TH3F *hDistortionR = new TH3F("hDistortionR", "Per-z-bin Distortion in the R direction as a function of (phi,r,z) (centered in r,phi, z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);
  TH3F *hDistortionP = new TH3F("hDistortionP", "Per-z-bin Distortion in the RPhi direction as a function of (phi,r,z)  (centered in r,phi, z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);
  TH3F *hDistortionZ = new TH3F("hDistortionZ", "Per-z-bin Distortion in the Z direction as a function of (phi,r,z)  (centered in r,phi, z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);
  TH3F *hIntDistortionR = new TH3F("hIntDistortionR", "Integrated R Distortion from (phi,r,z) to z=0 (centered in r,phi, and z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);
  TH3F *hIntDistortionP = new TH3F("hIntDistortionP", "Integrated R Distortion from (phi,r,z) to z=0 (centered in r,phi, and z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);
  TH3F *hIntDistortionZ = new TH3F("hIntDistortionZ", "Integrated R Distortion from (phi,r,z) to z=0  (centered in r,phi, and z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);

  TH3F *hIntDistortionX = new TH3F("hIntDistortionX", "Integrated X Distortion from (phi,r,z) to z=0 (centered in r,phi, and z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);
  TH3F *hIntDistortionY = new TH3F("hIntDistortionY", "Integrated Y Distortion from (phi,r,z) to z=0 (centered in r,phi, and z);phi;r;z", nph, pih, pfh, nrh, rih, rfh, nzh, zih, zfh);

  /*
  TH3F* hNewIntDistortionR=new TH3F("hNewIntDistortionR","Recursively Integrated R Distortion from (r,phi,z) to z=0 (centered in r,phi, and z);phi;r;z",nph,pih,pfh,nrh,rih,rfh,nzh,zih,zfh);
  TH3F* hNewIntDistortionP=new TH3F("hNewIntDistortionP","Recursively Integrated R Distortion from (r,phi,z) to z=0 (centered in r,phi, and z);phi;r;z",nph,pih,pfh,nrh,rih,rfh,nzh,zih,zfh);
  TH3F* hNewIntDistortionZ=new TH3F("hNewIntDistortionZ","Recursively Integrated R Distortion from (r,phi,z) to z=0  (centered in r,phi, and z);phi;r;z",nph,pih,pfh,nrh,rih,rfh,nzh,zih,zfh);
//for the interchanging distortion maps.
    TH2F* hNewIntDistortionR=new TH3F("hNewIntDistortionR","Recursively Integrated R Distortion from (r,phi,z) to z=0 (centered in r,phi, and z);phi;r;z",nph,pih,pfh,nrh,rih,rfh,nzh,zih,zfh);
  TH3F* hNewIntDistortionP=new TH3F("hNewIntDistortionP","Recursively Integrated R Distortion from (r,phi,z) to z=0 (centered in r,phi, and z);phi;r;z",nph,pih,pfh,nrh,rih,rfh,nzh,zih,zfh);
  TH3F* hNewIntDistortionZ=new TH3F("hNewIntDistortionZ","Recursively Integrated R Distortion from (r,phi,z) to z=0  (centered in r,phi, and z);phi;r;z",nph,pih,pfh,nrh,rih,rfh,nzh,zih,zfh);
  */

  // monitor plots, and the position that plot monitors at:

  // TVector3 pos((nrh/2+0.5)*steplocal.Perp()+rih,0,(nzh/2+0.5)*steplocal.Z()+zih);
  TVector3 pos(((int) (nrh / 2) + 0.5) * steplocal.Perp() + rih, 0, zmin + ((int) (nz / 2) + 0.5) * step.Z());
  float posphi = ((int) (nph / 2) + 0.5) * steplocal.Phi() + pih;
  pos.SetPhi(posphi);
  // int xi[3]={nrh/2,nph/2,nzh/2};
  int xi[3] = {(int) floor((pos.Perp() - rih) / steplocal.Perp()), (int) floor((posphi - pih) / steplocal.Phi()), (int) floor((pos.Z() - zih) / steplocal.Z())};
  if (!hasTwin)
  {
    std::cout << boost::str(boost::format("rpz slice indices= (%d,%d,%d) (no twin)") % xi[0] % xi[1] % xi[2]) << std::endl;
  }
  int twinz = 0;  // this is meant to be the matching position to xi[2] in the twin, hence generally -1*pos.  Better to just ask the twin rather than trying to calculate it ourselves...
  if (hasTwin)
  {
    twinz = twin->GetZindex(-1 * pos.Z());
    std::cout << boost::str(boost::format("rpz slice indices= (%d,%d,%d) twinz=%d") % xi[0] % xi[1] % xi[2] % twinz) << std::endl;
  }

  const std::string axname[] = {"r", "p", "z", "r", "p", "z"};
  int axn[] = {nrh, nph, nzh, nrh, nph, nzh};
  float axval[] = {(float) pos.Perp(), (float) pos.Phi(), (float) pos.Z(), (float) pos.Perp(), (float) pos.Phi(), (float) pos.Z()};
  float axbot[] = {rih, pih, zih, rih, pih, zih};
  float axtop[] = {rfh, pfh, zfh, rfh, pfh, zfh};
  TH2F *hIntDist[3][3];
  TH1F *hRDist[2][3];  // now with a paired friend for the negative side
  TH2F *hDiffDist[3][3];
  TH1F *hRDiffDist[2][3];
  for (int i = 0; i < 3; i++)
  {
    // loop over which axis of the distortion to read
    for (int ax = 0; ax < 3; ax++)
    {
      // loop over which plane to work in
      hDiffDist[ax][i] = new TH2F(std::string("hDiffDist" + axname[i] + "_" + axname[ax + 1] + axname[ax + 2]).c_str(),
                                  boost::str(boost::format("%s component of diff. distortion in  %s%s plane at %s=%2.3f;%s;%s") % axname[i] % axname[ax + 1] % axname[ax + 2] % axname[ax] % axval[ax] % axname[ax + 1] % axname[ax + 2]).c_str(),
                                  axn[ax + 1], axbot[ax + 1], axtop[ax + 1],
                                  axn[ax + 2], axbot[ax + 2], axtop[ax + 2]);
      hIntDist[ax][i] = new TH2F(std::string("hDIntDist" + axname[i] + "_" + axname[ax + 1] + axname[ax + 2]).c_str(),
                                 boost::str(boost::format("%s component of int. distortion in  %s%s plane at %s=%2.3f;%s;%s") % axname[i] % axname[ax + 1] % axname[ax + 2] % axname[ax] % axval[ax] % axname[ax + 1] % axname[ax + 2]).c_str(),
                                 axn[ax + 1], axbot[ax + 1], axtop[ax + 1],
                                 axn[ax + 2], axbot[ax + 2], axtop[ax + 2]);
    }
    hRDist[0][i] = new TH1F(std::string("hRDist" + axname[i]).c_str(),
                            boost::str(boost::format("%s component of int. distortion vs r with %s=%2.3f and %s=%2.3f;r(cm);#delta (cm)") % axname[i] % axname[1] % axval[1] % axname[2] % axval[2]).c_str(),
                            axn[0], axbot[0], axtop[0]);
    hRDist[1][i] = new TH1F(std::string("hRDist" + axname[i] + "Neg").c_str(),
                            boost::str(boost::format("%s component of int. distortion vs r with %s=%2.3f and %s=-%2.3f;r(cm);#delta (cm)") % axname[i] % axname[1] % axval[1] % axname[2] % axval[2]).c_str(),
                            axn[0], axbot[0], axtop[0]);
    hRDiffDist[0][i] = new TH1F(std::string("hRDist" + axname[i]).c_str(),
                                boost::str(boost::format("%s component of diff. distortion vs r with %s=%2.3f and %s=%2.3f;r(cm);#delta (cm)") % axname[i] % axname[1] % axval[1] % axname[2] % axval[2]).c_str(),
                                axn[0], axbot[0], axtop[0]);
    hRDiffDist[1][i] = new TH1F(std::string("hRDist" + axname[i] + "Neg").c_str(),
                                boost::str(boost::format("%s component of diff. distortion vs r with %s=%2.3f and %s=-%2.3f;r(cm);#delta (cm)") % axname[i] % axname[1] % axval[1] % axname[2] % axval[2]).c_str(),
                                axn[0], axbot[0], axtop[0]);
  }

  TVector3 inpart, outpart;
  TVector3 distort;
  int validToStep;

  // TTree version:
  float partR, partP, partZ;
  int ir, ip, iz;
  float distortR, distortP, distortZ;
  float distortX, distortY;
  float diffdistR, diffdistP, diffdistZ;
  TTree *dTree = new TTree("dTree", "Distortion per step z");
  dTree->Branch("r", &partR);
  dTree->Branch("p", &partP);
  dTree->Branch("z", &partZ);
  dTree->Branch("ir", &ir);
  dTree->Branch("ip", &ip);
  dTree->Branch("iz", &iz);
  dTree->Branch("dr", &distortR);
  dTree->Branch("drp", &distortP);
  dTree->Branch("dz", &distortZ);

  std::cout << "generating distortion map with (" << nrh << "x" << nph << "x" << nzh << " grid" << std::endl;
  unsigned long long totalelements = nrh;
  totalelements *= nph;
  totalelements *= nzh;  // breaking up this multiplication prevents a 32bit math overflow
  unsigned long long percent = totalelements / 100 * debug_npercent;
  std::cout << boost::str(boost::format("total elements = %llu") % totalelements) << std::endl;

  int el = 0;

  // we want to loop over the entire region to be mapped, but we also need to include
  // one additional bin at each edge, to allow the mc drift code to interpolate properly.
  // hence we count from -1 to n+1, and manually adjust the position in those edge cases
  // to avoid sampling nonphysical regions in r and z.  the phi case is free to wrap as
  //  normal.

  // note that we apply the adjustment to the particle position (inpart) and not the plotted position (partR etc)
  inpart.SetXYZ(1, 0, 0);
  for (ir = 0; ir < nrh; ir++)
  {
    partR = (ir + 0.5) * deltar + rih;
    if (ir == 0)
    {
      inpart.SetPerp(partR + deltar);
    }
    else if (ir == nrh - 1)
    {
      inpart.SetPerp(partR - deltar);
    }
    else
    {
      inpart.SetPerp(partR);
    }
    for (ip = 0; ip < nph; ip++)
    {
      partP = (ip + 0.5) * deltap + pih;
      inpart.SetPhi(partP);
      // since phi loops, there's no need to adjust phis that are out of bounds.
      for (iz = 0; iz < nzh; iz++)
      {
        partZ = (iz) *deltaz + zih;  // start us at the EDGE of the bin, maybe has problems at the CM when twinned.
        if (iz == 0)
        {
          inpart.SetZ(partZ + deltaz);
        }
        else if (iz == nzh - 1)
        {
          inpart.SetZ(partZ - deltaz);
        }
        else
        {
          inpart.SetZ(partZ);
        }
        partZ += 0.5 * deltaz;  // move to center of histogram bin.

        // print_need_cout("iz=%d, zcoord=%2.2f, bin=%d\n",iz,partZ,  hIntDist[0][0]->GetYaxis()->FindBin(partZ));

        // differential distortion:
        // be careful with the math of a distortion.  The R distortion is NOT the perp() component of outpart-inpart -- that's the transverse magnitude of the distortion!
        if (hasTwin && inpart.Z() < 0)
        {
          distort = twin->GetTotalDistortion(inpart.Z(), inpart + stepzvec, nSteps, true, &validToStep, &successCheck);  // step across the cell in the opposite direction, starting at the high side and going to the low side..
        }
        else
        {
          distort = GetTotalDistortion(inpart.Z() + deltaz, inpart, nSteps, true, &validToStep, &successCheck);
        }
        distort.RotateZ(-inpart.Phi());  // rotate so that that is on the x axis
        diffdistP = distort.Y();         // the phi component is now the y component.
        diffdistR = distort.X();         // and the r component is the x component
        diffdistZ = distort.Z();
        hDistortionR->Fill(partP, partR, partZ, diffdistR);
        hDistortionP->Fill(partP, partR, partZ, diffdistP);
        hDistortionZ->Fill(partP, partR, partZ, diffdistZ);
        dTree->Fill();

        // integral distortion:
        if (hasTwin && makeUnifiedMap && inpart.Z() < 0)
        {
          distort = twin->GetTotalDistortion(-z_readout, inpart + stepzvec, nSteps, true, &validToStep, &successCheck);
        }
        else
        {
          distort = GetTotalDistortion(z_readout, inpart, nSteps, true, &validToStep, &successCheck);
        }
        distortX = distort.X();
        distortY = distort.Y();
        distort.RotateZ(-inpart.Phi());  // rotate so that that is on the x axis
        distortP = distort.Y();          // the phi component is now the y component.
        distortR = distort.X();          // and the r component is the x component
        distortZ = distort.Z();

        // recursive integral distortion:
        // get others working first!

        // print_need_cout("part=(rpz)(%f,%f,%f),distortP=%f\n",partP,partR,partZ,distortP);
        hIntDistortionR->Fill(partP, partR, partZ, distortR);
        hIntDistortionP->Fill(partP, partR, partZ, distortP);
        hIntDistortionZ->Fill(partP, partR, partZ, distortZ);

        if (andCartesian)
        {
          hIntDistortionX->Fill(partP, partR, partZ, distortX);
          hIntDistortionY->Fill(partP, partR, partZ, distortY);
        }

        // now we fill particular slices for integral visualizations:
        if (ir == xi[0])
        {  // r slice
          // print_need_cout("ir=%d, r=%f (pz)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",ir,partR,ip,iz,distortR,distortP);
          hIntDist[0][0]->Fill(partP, partZ, distortR);
          hIntDist[0][1]->Fill(partP, partZ, distortP);
          hIntDist[0][2]->Fill(partP, partZ, distortZ);
          hDiffDist[0][0]->Fill(partP, partZ, diffdistR);
          hDiffDist[0][1]->Fill(partP, partZ, diffdistP);
          hDiffDist[0][2]->Fill(partP, partZ, diffdistZ);
        }
        if (ip == xi[1])
        {  // phi slice
          // print_need_cout("ip=%d, p=%f (rz)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",ip,partP,ir,iz,distortR,distortP);
          hIntDist[1][0]->Fill(partZ, partR, distortR);
          hIntDist[1][1]->Fill(partZ, partR, distortP);
          hIntDist[1][2]->Fill(partZ, partR, distortZ);
          hDiffDist[1][0]->Fill(partZ, partR, diffdistR);
          hDiffDist[1][1]->Fill(partZ, partR, diffdistP);
          hDiffDist[1][2]->Fill(partZ, partR, diffdistZ);

          if (iz == xi[2])
          {  // z slices of phi slices= r line at mid phi, mid z:
            hRDist[0][0]->Fill(partR, distortR);
            hRDist[0][1]->Fill(partR, distortP);
            hRDist[0][2]->Fill(partR, distortZ);
            hRDiffDist[0][0]->Fill(partR, diffdistR);
            hRDiffDist[0][1]->Fill(partR, diffdistP);
            hRDiffDist[0][2]->Fill(partR, diffdistZ);
          }
          if (hasTwin && iz == twinz)
          {  // z slices of phi slices= r line at mid phi, mid z:
            hRDist[1][0]->Fill(partR, distortR);
            hRDist[1][1]->Fill(partR, distortP);
            hRDist[1][2]->Fill(partR, distortZ);
            hRDiffDist[1][0]->Fill(partR, diffdistR);
            hRDiffDist[1][1]->Fill(partR, diffdistP);
            hRDiffDist[1][2]->Fill(partR, diffdistZ);
          }
        }
        if (iz == xi[2])
        {  // z slice
          // print_need_cout("iz=%d, z=%f (rp)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",iz,partZ,ir,ip,distortR,distortP);

          hIntDist[2][0]->Fill(partR, partP, distortR);
          hIntDist[2][1]->Fill(partR, partP, distortP);
          hIntDist[2][2]->Fill(partR, partP, distortZ);
          hDiffDist[2][0]->Fill(partR, partP, diffdistR);
          hDiffDist[2][1]->Fill(partR, partP, diffdistP);
          hDiffDist[2][2]->Fill(partR, partP, diffdistZ);
        }

        if (!(el % percent))
        {
          std::cout << boost::str(boost::format("generating distortions %d%%:  ") % ((int) (debug_npercent * (el / percent))));
          std::cout << boost::str(boost::format("distortion at (ir=%d,ip=%d,iz=%d) is (%E,%E,%E)") % ir % ip % iz % distortR % distortP % distortZ) << std::endl;
        }
        el++;
      }
    }
  }

  TCanvas *canvas = new TCanvas("cdistort", "distortion integrals", 1200, 800);
  // take 10 of the bottom of this for data?
  canvas->cd();
  TPad *c = new TPad("cplots", "distortion integral plots", 0, 0.2, 1, 1);
  canvas->cd();
  TPad *textpad = new TPad("ctext", "distortion integral plots", 0, 0.0, 1, 0.2);
  c->Divide(4, 3);
  gStyle->SetOptStat();
  for (int i = 0; i < 3; i++)
  {
    // component
    for (int ax = 0; ax < 3; ax++)
    {
      // plane
      c->cd(i * 4 + ax + 1);
      gPad->SetRightMargin(0.15);
      hIntDist[ax][i]->SetStats(false);
      hIntDist[ax][i]->Draw("colz");
    }
    c->cd(i * 4 + 4);
    hRDist[0][i]->SetStats(false);
    hRDist[0][i]->SetFillColor(kRed);
    hRDist[0][i]->Draw("hist");
    if (hasTwin)
    {
      hRDist[1][i]->SetStats(false);
      hRDist[1][i]->SetLineColor(kBlue);
      hRDist[1][i]->Draw("hist,same");
    }
    /*
      double Cut = 40;
      h->SetFillColor(kRed);
      TH1F *hNeg = (TH1F*)hRDist[i]->Clone(Form("hNegRDist%d",i));
      hNeg->SetFillColor(kGreen);
      for (int n = 1; n <= hNeg->GetNbinsX(); n++) {
      hNeg->SetBinContent(n,Cut);
      }
      h3->Draw(); h.Draw("same");
      TH1F *h2 = (TH1F*)h->Clone("h2");
      h2->SetFillColor(kGray-4);
      h2->SetMaximum(Cut);
      h2->Draw("same");
    */
  }
  textpad->cd();
  float texpos = 0.9;
  float texshift = 0.12;
  TLatex *tex = new TLatex(0.0, texpos, "Fill Me In");
  tex->SetTextSize(texshift * 0.8);
  tex->DrawLatex(0.05, texpos, GetFieldString().c_str());
  texpos -= texshift;
  tex->DrawLatex(0.05, texpos, GetChargeString().c_str());
  texpos -= texshift;
  // tex->DrawLatex(0.05,texpos,Form("Drift Field = %2.2f V/cm",GetNominalE()));texpos-=texshift;
  tex->DrawLatex(0.05, texpos, boost::str(boost::format("Drifting grid of (rp)=(%d x %d) electrons with %d steps") % nrh % nph % nSteps).c_str());
  texpos -= texshift;
  tex->DrawLatex(0.05, texpos, GetLookupString().c_str());
  texpos -= texshift;
  tex->DrawLatex(0.05, texpos, GetGasString().c_str());
  texpos -= texshift;
  if (debug_distortionScale.Mag() > 0)
  {
    tex->DrawLatex(0.05, texpos, boost::str(boost::format("Distortion scaled by (r,p,z)=(%2.2f,%2.2f,%2.2f)") % debug_distortionScale.X() % debug_distortionScale.Y() % debug_distortionScale.Z()).c_str());
    texpos -= texshift;
  }
  texpos = 0.9;

  canvas->cd();
  c->Draw();
  canvas->cd();
  textpad->Draw();
  canvas->SaveAs(summaryFilename.c_str());

  // canvas->cd();
  // already done TPad *c=new TPad("cplots","distortion differential plots",0,0.2,1,1);
  // canvas->cd();
  // already done TPad *textpad=new TPad("ctext","distortion differential plots",0,0.0,1,0.2);
  // already done c->Divide(4,3);
  // gStyle->SetOptStat();
  for (int i = 0; i < 3; i++)
  {
    // component
    for (int ax = 0; ax < 3; ax++)
    {
      // plane
      c->cd(i * 4 + ax + 1);
      gPad->SetRightMargin(0.15);
      hDiffDist[ax][i]->SetStats(false);
      hDiffDist[ax][i]->Draw("colz");
    }
    c->cd(i * 4 + 4);
    hRDiffDist[0][i]->SetStats(false);
    hRDiffDist[0][i]->SetFillColor(kRed);
    hRDiffDist[0][i]->Draw("hist");
    if (hasTwin)
    {
      hRDiffDist[1][i]->SetStats(false);
      hRDiffDist[1][i]->SetLineColor(kBlue);
      hRDiffDist[1][i]->Draw("hist same");
    }
  }
  textpad->cd();
  texpos = 0.9;
  texshift = 0.12;
  tex->SetTextSize(texshift * 0.8);
  tex->DrawLatex(0.05, texpos, GetFieldString().c_str());
  texpos -= texshift;
  tex->DrawLatex(0.05, texpos, GetChargeString().c_str());
  texpos -= texshift;
  // tex->DrawLatex(0.05,texpos,Form("Drift Field = %2.2f V/cm",GetNominalE()));texpos-=texshift;
  tex->DrawLatex(0.05, texpos, boost::str(boost::format("Drifting grid of (rp)=(%d x %d) electrons with %d steps") % nrh % nph % nSteps).c_str());
  texpos -= texshift;
  tex->DrawLatex(0.05, texpos, GetLookupString().c_str());
  texpos -= texshift;
  tex->DrawLatex(0.05, texpos, GetGasString().c_str());
  texpos -= texshift;
  tex->DrawLatex(0.05, texpos, "Differential Plots");
  texpos -= texshift;
  if (debug_distortionScale.Mag() > 0)
  {
    tex->DrawLatex(0.05, texpos, boost::str(boost::format("Distortion scaled by (r,p,z)=(%2.2f,%2.2f,%2.2f)") % debug_distortionScale.X() % debug_distortionScale.Y() % debug_distortionScale.Z()).c_str());
    texpos -= texshift;
  }
  texpos = 0.9;

  canvas->cd();
  c->Draw();
  canvas->cd();
  textpad->Draw();
  canvas->SaveAs(diffSummaryFilename.c_str());

  //  print_need_cout("map:%s.\n",distortionFilename.c_str());

  outf->cd();

  hDistortionR->Write();
  hDistortionP->Write();
  hDistortionZ->Write();
  hIntDistortionR->Write();
  hIntDistortionP->Write();
  hIntDistortionZ->Write();
  if (andCartesian)
  {
    hIntDistortionX->Write();
    hIntDistortionY->Write();
  }
  dTree->Write();
  outf->Close();
  // print_need_cout("map:%s.closed\n",distortionFilename.c_str());

  /*
  //all histograms associated with this file should be deleted when we closed the file.
  //delete the histograms we made:
  TH3F *m;
  m = (TH3F*)gROOT­>FindObject("hDistortionR");  if(m) { print_need_cout("found in memory still.\n"); m­>Delete();}
  m = (TH3F*)gROOT­>FindObject("hDistortionP");  if(m) m­>Delete();
  m = (TH3F*)gROOT­>FindObject("hDistortionZ");  if(m) m­>Delete();
  m = (TH3F*)gROOT­>FindObject("hIntDistortionR");  if(m) m­>Delete();
  m = (TH3F*)gROOT­>FindObject("hIntDistortionP");  if(m) m­>Delete();
  m = (TH3F*)gROOT­>FindObject("hIntDistortionZ");  if(m) m­>Delete();

  print_need_cout("map:%s.deleted via convoluted call.  sigh.\n",distortionFilename.c_str());
  for (int i=0;i<3;i++){
    //loop over which axis of the distortion to read
    for (int ax=0;ax<3;ax++){
      //loop over which plane to work in
      hIntDist[ax][i]->Delete();
    }
    hRDist[i]->Delete();
  }
  */
  std::cout << "wrote map and summary to " << filebase << std::endl;
  std::cout << "map: " << distortionFilename << std::endl;
  std::cout << "summary: " << summaryFilename << std::endl;
  // print_need_cout("wrote map and summary to %s and %s.\n",distortionFilename.c_str(),summaryFilename.c_str());
  return;
}

TVector3 AnnularFieldSim::swimTo(float zdest, const TVector3 &start, bool interpolate, bool useAnalytic)
{
  int defaultsteps = 100;
  int goodtostep = 0;
  if (useAnalytic)
  {
    return swimToInAnalyticSteps(zdest, start, defaultsteps, &goodtostep);
  }
  return swimToInSteps(zdest, start, defaultsteps, interpolate, &goodtostep);
}

TVector3 AnnularFieldSim::GetStepDistortion(float zdest, const TVector3 &start, bool interpolate, bool useAnalytic)
{
  // getting the distortion instead of the post-step position allows us to accumulate small deviations from the original position that might be lost in the large number

  // using second order langevin expansion from http://skipper.physics.sunysb.edu/~prakhar/tpc/Papers/ALICE-INT-2010-016.pdf
  // TVector3 (*field)[nr][ny][nz]=field_;
  int rt, pt, zt;  // these are filled by the checkbounds that follow, but are not used.
  BoundsCase zBound = GetZindexAndCheckBounds(start.Z(), &zt);
  if (GetRindexAndCheckBounds(start.Perp(), &rt) != InBounds || GetPhiIndexAndCheckBounds(FilterPhiPos(start.Phi()), &pt) != InBounds || (zBound != InBounds && zBound != OnHighEdge))
  {
    // print_need_cout("AnnularFieldSim::swimTo asked to swim particle from (xyz)=(%f,%f,%f) which is outside the ROI:\n", start.X(), start.Y(), start.Z());
    // print_need_cout(" -- %f <= r < %f \t%f <= phi < %f \t%f <= z < %f \n", rmin_roi * step.Perp(), rmax_roi * step.Perp(), phimin_roi * step.Phi(), phimax_roi * step.Phi(), zmin_roi * step.Z(), zmax_roi * step.Z());
    // print_need_cout("Returning original position.\n");
    return start;
  }

  // set the direction of the external fields.

  double zdist = zdest - start.Z();

  // short-circuit if there's no travel length:

  if (fabs(zdist) < ALMOST_ZERO * step.Z())
  {
    //   std::cout <<  boost::str(boost::format("Asked  particle from (%f,%f,%f) to z=%f, which is a distance of %fcm.  Returning zero.") %start.X() %start.Y() %start.Z() %zdest %zdist) << std::endl;
    return zero_vector;
  }

  TVector3 fieldInt;   // integral of E field along path
  TVector3 fieldIntB;  // integral of B field along path

  // note that using analytic takes priority over interpolating todo: clean this up to use a status rther than a pair of flags
  if (useAnalytic)
  {
    fieldInt = analyticFieldIntegral(zdest, start, Efield);
    fieldIntB = analyticFieldIntegral(zdest, start, Bfield);
  }
  else if (interpolate)
  {
    fieldInt = interpolatedFieldIntegral(zdest, start, Efield);
    fieldIntB = interpolatedFieldIntegral(zdest, start, Bfield);
  }
  else
  {
    fieldInt = fieldIntegral(zdest, start, Efield);
    fieldIntB = fieldIntegral(zdest, start, Bfield);
  }

  if (abs(fieldInt.Z() / zdist) < ALMOST_ZERO)
  {
    std::cout << "GetStepDistortion is attempting to swim with no drift field:" << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: (%2.4f,%2.4f,%2.4f) to z=%2.4f") % start.X() % start.Y() % start.Z() % zdest) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: fieldInt=(%E,%E,%E)") % fieldInt.X() % fieldInt.Y() % fieldInt.Z()) << std::endl;
    exit(1);
  }
  // float fieldz=field_[in3(x,y,0,fx,fy,fz)].Z()+E.Z();// *field[x][y][zi].Z();
  double EfieldZ = fieldInt.Z() / zdist;  // average field over the path.
  double BfieldZ = fieldIntB.Z() / zdist;
  // double fieldz=Enominal; // ideal field over path.

  // these values should be with real, not nominal field?
  // double mu=abs(vdrift/Enominal);//vdrift in [cm/s], field in [V/cm] hence mu in [cm^2/(V*s)];  should be a positive value.  drift velocity over field magnitude, not field direction.
  // double omegatau=-mu*Bnominal;//minus sign is for electron charge.
  double omegatau = omegatau_nominal;  // don't compute this every time!
  // or:  omegatau=-10*(10*B.Z()/Tesla)*(vdrift/(cm/us))/(fieldz/(V/cm)); //which is the same as my calculation up to a sign.
  // std::cout << "omegatau=" << omegatau << std::endl;

  double T1om = langevin_T1 * omegatau;
  double T2om2 = langevin_T2 * omegatau * langevin_T2 * omegatau;
  double c0 = 1 / (1 + T2om2);           //
  double c1 = T1om / (1 + T1om * T1om);  // Carlos gets this term wrong.  It should be linear on top, quadratic on bottom.
  double c2 = T2om2 / (1 + T2om2);

  TVector3 EintOverEz = 1 / EfieldZ * fieldInt;   // integral of E/E_z= integral of E / integral of E_z * delta_z
  TVector3 BintOverBz = 1 / BfieldZ * fieldIntB;  // should this (and the above?) be BfieldZ or Bnominal?

  // really this should be the integral of the ratio, not the ratio of the integrals.
  // and should be integrals over the B field, btu for now that's fixed and constant across the region, so not necessary
  // there's no reason to do this as r phi.  This is an equivalent result, since I handle everything in vectors.
  double deltaX = c0 * EintOverEz.X() + c1 * EintOverEz.Y() - c1 * BintOverBz.Y() + c2 * BintOverBz.X();
  double deltaY = c0 * EintOverEz.Y() - c1 * EintOverEz.X() + c1 * BintOverBz.X() + c2 * BintOverBz.Y();
  // strictly, for deltaZ we want to integrate v'(E)*(E-E0)dz and v''(E)*(E-E0)^2 dz, but over a short step the field is constant, and hence this can be a product of the integral and not an integral of the product:

  double vprime = (5000 * cm / s) / (100 * V / cm);  // hard-coded value for 50:50.  Eventually this needs to be part of the constructor, as do most of the repeated math terms here.
  double vdoubleprime = 0;                           // neglecting the v'' term for now.  Fair? It's pretty linear at our operating point, and it would require adding an additional term to the field integral.

  // note: as long as my step is very small, I am essentially just reading the field at a point and multiplying by the step size.
  // hence integral of P dz, where P is a function of the fields, int|P(E(x,y,z))dz=P(int|E(x,y,z)dz/deltaZ)*deltaZ
  // hence: , eg, int|E^2dz=(int|Edz)^2/deltaz
  double deltaZ = vprime / vdrift * (fieldInt.Z() - zdist * Enominal) + vdoubleprime / vdrift * (fieldInt.Z() - Enominal * zdist) * (fieldInt.Z() - Enominal * zdist) / (2 * zdist) - 0.5 / zdist * (EintOverEz.X() * EintOverEz.X() + EintOverEz.Y() * EintOverEz.Y()) + c1 / zdist * (EintOverEz.X() * BintOverBz.Y() - EintOverEz.Y() * BintOverBz.X()) + c2 / zdist * (EintOverEz.X() * BintOverBz.X() + EintOverEz.Y() * BintOverBz.Y()) + c2 / zdist * (BintOverBz.X() * BintOverBz.X() + BintOverBz.Y() * BintOverBz.Y());  // missing v'' term.

  if (false)
  {
    std::cout << boost::str(boost::format("GetStepDistortion:  (c0,c1,c2)=(%E,%E,%E)") % c0 % c1 % c2) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion:  EintOverEz==(%E,%E,%E)") % EintOverEz.X() % EintOverEz.Y() % EintOverEz.Z()) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion:  BintOverBz==(%E,%E,%E)") % BintOverBz.X() % BintOverBz.Y() % BintOverBz.Z()) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: (%2.4f,%2.4f,%2.4f) to z=%2.4f") % start.X() % start.Y() % start.Z() % zdest) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: fieldInt=(%E,%E,%E)") % fieldInt.X() % fieldInt.Y() % fieldInt.Z()) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: delta=(%E,%E,%E)") % deltaX % deltaY % deltaZ) << std::endl;
  }

  // if (abs(deltaZ / zdist) > 0.25)
  if (false)
  {
    std::cout << "GetStepDistortion produced a very large zdistortion!" << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: zdist=%2.4f, deltaZ=%2.4f, Delta/z=%2.4f") % zdist % deltaZ % (deltaZ / zdist)) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: average field Z:  Ez=%2.4fV/cm, Bz=%2.4fT") % (EfieldZ / (V / cm)) % (BfieldZ / Tesla)) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion:  (c0,c1,c2)=(%E,%E,%E)") % c0 % c1 % c2) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion:  EintOverEz==(%E,%E,%E)") % EintOverEz.X() % EintOverEz.Y() % EintOverEz.Z()) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion:  BintOverBz==(%E,%E,%E)") % BintOverBz.X() % BintOverBz.Y() % BintOverBz.Z()) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: (%2.4f,%2.4f,%2.4f) to z=%2.4f") % start.X() % start.Y() % start.Z() % zdest) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: EfieldInt=(%E,%E,%E)") % fieldInt.X() % fieldInt.Y() % fieldInt.Z()) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: BfieldInt=(%E,%E,%E)") % fieldIntB.X() % fieldIntB.Y() % fieldIntB.Z()) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: delta=(%E,%E,%E)") % deltaX % deltaY % deltaZ) << std::endl;
    // assert(false);
  }

  if (abs(deltaX) < 1E-20 && !(chargeCase == NoSpacecharge))
  {
    std::cout << boost::str(boost::format("GetStepDistortion produced a very small deltaX: %E") % deltaX) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion:  (c0,c1,c2)=(%E,%E,%E)") % c0 % c1 % c2) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion:  EintOverEz==(%E,%E,%E)") % EintOverEz.X() % EintOverEz.Y() % EintOverEz.Z()) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion:  BintOverBz==(%E,%E,%E)") % BintOverBz.X() % BintOverBz.Y() % BintOverBz.Z()) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: (%2.4f,%2.4f,%2.4f) to z=%2.4f") % start.X() % start.Y() % start.Z() % zdest) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: fieldInt=(%E,%E,%E)") % fieldInt.X() % fieldInt.Y() % fieldInt.Z()) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: delta=(%E,%E,%E)") % deltaX % deltaY % deltaZ) << std::endl;
    // assert(1==2);
  }

  // if (!(abs(deltaX) < 1E3))
  if (false)
  {
    std::cout << boost::str(boost::format("GetStepDistortion produced a very large deltaX: %E") % deltaX) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion:  (c0,c1,c2)=(%E,%E,%E)") % c0 % c1 % c2) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion:  EintOverEz==(%E,%E,%E)") % EintOverEz.X() % EintOverEz.Y() % EintOverEz.Z()) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion:  BintOverBz==(%E,%E,%E)") % BintOverBz.X() % BintOverBz.Y() % BintOverBz.Z()) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: (%2.4f,%2.4f,%2.4f) (rp)=(%2.4f,%2.4f) to z=%2.4f") % start.X() % start.Y() % start.Z() % start.Perp() % start.Phi() % zdest) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: fieldInt=(%E,%E,%E)") % fieldInt.X() % fieldInt.Y() % fieldInt.Z()) << std::endl;
    std::cout << boost::str(boost::format("GetStepDistortion: delta=(%E,%E,%E)") % deltaX % deltaY % deltaZ) << std::endl;
    exit(1);
  }

  // deltaZ=0;//temporary removal.

  TVector3 shift(deltaX, deltaY, deltaZ);
  if (debug_distortionScale.Mag() > 0)
  {  // debug code to scale the resulting distortions
    shift.RotateZ(-start.Phi());
    // TVector3 localScale=debug_distortionScale;
    // localScale.RotateZ(start.Phi());
    shift.SetXYZ(shift.X() * debug_distortionScale.X(), shift.Y() * debug_distortionScale.Y(), shift.Z() * debug_distortionScale.Z());
    shift.RotateZ(start.Phi());
  }

  return shift;
}

// putting all the getters here out of the way:
const std::string AnnularFieldSim::GetLookupString()
{
  if (lookupCase == LookupCase::Full3D)
  {
    return boost::str(boost::format("Full3D (%d x %d x %d) with (%d x %d x %d) roi") % nr % nphi % nz % nr_roi % nphi_roi % nz_roi);
  }

  if (lookupCase == LookupCase::PhiSlice)
  {
    return boost::str(boost::format("PhiSlice (%d x %d x %d) with (%d x 1 x %d) roi") % nr % nphi % nz % nr_roi % nz_roi);
  }

  return "broken";
}
const std::string AnnularFieldSim::GetGasString()
{
  return boost::str(boost::format("vdrift=%2.2fcm/us, Enom=%2.2fV/cm, Bnom=%2.2fT, omtau=%2.4E") % (vdrift / (cm / us)) % (Enominal / (V / cm)) % (Bnominal / Tesla) % omegatau_nominal);
}

const std::string AnnularFieldSim::GetFieldString()
{
  return boost::str(boost::format("%s, %s") % Efieldname % Bfieldname);
}

// NOLINTNEXTLINE(misc-no-recursion)
TVector3 AnnularFieldSim::GetFieldAt(const TVector3 &pos)
{
  // assume pos is in native units (see header)

  int r, p, z;

  if (GetRindexAndCheckBounds(pos.Perp(), &r) == BoundsCase::OutOfBounds)
  {
    return zero_vector;
  }
  if (GetPhiIndexAndCheckBounds(FilterPhiPos(pos.Phi()), &p) == BoundsCase::OutOfBounds)
  {
    return zero_vector;
  }
  if (GetZindexAndCheckBounds(pos.Z(), &z) == BoundsCase::OutOfBounds)
  {
    if (hasTwin)
    {
      return twin->GetFieldAt(pos);
    }
    return zero_vector;
  }
  return Efield->Get(r, p, z);
}

// NOLINTNEXTLINE(misc-no-recursion)
TVector3 AnnularFieldSim::GetBFieldAt(const TVector3 &pos)
{
  // assume pos is in native units (see header)

  int r, p, z;

  if (GetRindexAndCheckBounds(pos.Perp(), &r) == BoundsCase::OutOfBounds)
  {
    return zero_vector;
  }
  if (GetPhiIndexAndCheckBounds(FilterPhiPos(pos.Phi()), &p) == BoundsCase::OutOfBounds)
  {
    return zero_vector;
  }
  if (GetZindexAndCheckBounds(pos.Z(), &z) == BoundsCase::OutOfBounds)
  {
    if (hasTwin)
    {
      return twin->GetBFieldAt(pos);
    }
    return zero_vector;
  }
  return Bfield->Get(r, p, z);
}

// NOLINTNEXTLINE(misc-no-recursion)
float AnnularFieldSim::GetChargeAt(const TVector3 &pos)
{
  int z;
  BoundsCase zbound = GetZindexAndCheckBounds(pos.Z(), &z);  //==BoundsCase::OutOfBounds) return zero_vector;
  if (zbound == OutOfBounds)
  {
    if (hasTwin)
    {
      return twin->GetChargeAt(pos);
    }
    std::cout << boost::str(boost::format("Caution:  tried to read charge at zbin=%d!  No twin available to handle this") % z) << std::endl;
    return -999;
  }

  return q->GetChargeAtPosition(pos.Perp(), FilterPhiPos(pos.Phi()), pos.Z());  // because tvectors take position to be -phi to phi, we always have to filter.
  // actually, we should probably just yield to that assumption in more places to speed this up.
  /*
  //assume pos is in native units (see header)
  int r, p, z;

  //get the bounds, but we don't want to actually check the cases, because the charge can go outside the vector region.
  r = GetRindex(pos.Perp());
  p = GetPhiIndex(pos.Phi());

  BoundsCase zbound = GetZindexAndCheckBounds(pos.Z(), &z);  //==BoundsCase::OutOfBounds) return zero_vector;
  if (zbound == OutOfBounds && hasTwin) return twin->GetChargeAt(pos);
  return q->Get(r, p, z);
  */
}
