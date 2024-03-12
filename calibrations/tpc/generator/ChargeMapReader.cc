#include "ChargeMapReader.h"

#include "MultiArray.h"

#include <TAxis.h>
#include <TFile.h>
#include <TH2.h>
#include <TH3.h>

#include <boost/format.hpp>

#include <cassert>
#include <cmath>
#include <iostream>

#define DEBUG false

ChargeMapReader::ChargeMapReader()
  : ChargeMapReader(20, 20.0, 78.0, 20, 0, 2 * M_PI, 40, -105.5, 105.5)
{
  std::cout << "made a new ChargeMapReader with default values -- cascading to next constructor" << std::endl;
  return;
}

ChargeMapReader::ChargeMapReader(int _nr, float _rmin, float _rmax, int _nphi, float _phimin, float _phimax, int _nz, float _zmin, float _zmax)
{
  std::cout << boost::str(boost::format("made a new ChargeMapReader with defined values:\n %d %.2f %.2f\n %d %.2f %.2f\n %d %.2f %.2f") % _nr % _rmin % _rmax % _nphi % _phimin % _phimax % _nz % _zmin % _zmax) << std::endl;
  SetOutputParameters(_nr, _rmin, _rmax, _nphi, _phimin, _phimax, _nz, _zmin, _zmax);
  return;
}

bool ChargeMapReader::CanInterpolateAt(float r, float phi, float z)
{
  return CanInterpolateAt(phi, r, z, hChargeDensity);
}

// a convenient method to check whether it's safe to interpolate for a particular histogram.
bool ChargeMapReader::CanInterpolateAt(float x, float y, float z, TH3* h)
{
  if (h == nullptr)
  {
    return false;
  }
  float pos[3] = {x, y, z};
  // todo: is it worth keeping these values somewhere for ease of access?
  TAxis* ax[3] = {nullptr, nullptr, nullptr};
  ax[0] = h->GetXaxis();
  ax[1] = h->GetYaxis();
  ax[2] = h->GetZaxis();

  int nbins[3];
  for (int i = 0; i < 3; i++)
  {
    nbins[i] = ax[i]->GetNbins();  // number of bins, not counting under and overflow.
    if (nbins[i] == 1)
    {
      return false;  // if there's only one non over/under bin, then no place is safe to interpolate.
    }
  }

  //   0     1     2   ...   n-1    n    n+1
  // under|first|   ..|  ..  |  .. |last| over
  for (int i = 0; i < 3; i++)
  {  // check each axis:
    int axbin = ax[i]->FindBin(pos[i]);

    if (axbin < 1 || axbin > nbins[i])
    {
      return false;  // before the first bin, or after the last bin
    }
    if (axbin > 1 && axbin < nbins[i])
    {
      continue;  // if we're in a middle bin, we're fine on this axis.
    }

    // now we need to check if we're in the safe parts of the first and last bins:
    float low = ax[i]->GetBinLowEdge(axbin);
    float high = ax[i]->GetBinUpEdge(axbin);
    float binmid = 0.5 * (high + low);

    if (axbin == 1 && pos[i] <= binmid)
    {
      return false;  // we're in the first bin, but below the midpoint, so we would interpolate out of bounds
    }
    if (axbin == nbins[i] && pos[i] >= binmid)
    {
      return false;  // we're in the last bin, but above the midpoint, so we would interpolate out of bounds
    }
  }

  // if we get here, we are in the okay-range of all three axes.
  return true;
}

void ChargeMapReader::FillChargeHistogram(TH3* h)
{
  // fills the provided histogram with the data from the internal representation.

  //   0     1     2     ...   n-1
  // first|   ..|  ..  |  .. |last
  float dphi, dr, dz;  // bin widths in each dimension.  Got too confusing to make these an array.
  if (DEBUG)
  {
    std::cout << "filling chargehistogram" << std::endl;
  }

  dr = binWidth[0];
  dphi = binWidth[1];
  dz = binWidth[2];

  float phimid, zmid;  // midpoints at each step.
  int i[3];
  for (i[0] = 0; i[0] < nBins[0]; i[0]++)
  {  // r
    float rmid = lowerBound[0] + (i[0] + 0.5) * dr;
    for (i[1] = 0; i[1] < nBins[1]; i[1]++)
    {  // phi
      phimid = lowerBound[1] + (i[1] + 0.5) * dphi;
      for (i[2] = 0; i[2] < nBins[2]; i[2]++)
      {  // z
        zmid = lowerBound[2] + (i[2] + 0.5) * dz;
        h->Fill(phimid, rmid, zmid, charge->Get(i[0], i[1], i[2]));
      }  // z
    }    // phi
  }      // r
  return;
}

void ChargeMapReader::RegenerateCharge()
{
  // Builds the charge 3D array (internal representation) from the internal charge density map.
  // either the density map has changed, or the binning of the output has changed (hopefully not the latter, because that's a very unusual thing to change mid-run.
  // we want to rebuild the charge per bin of our output representation in any case.  Generally, we will interpolate from the charge density that we know we have, but we need to be careful not to ask to interpolate in regions where that is not allowed.
  if (DEBUG)
  {
    std::cout << boost::str(boost::format("regenerating charge array contents with axis scale=%1.2E and charge scale=%1.2E") % inputAxisScale % inputChargeScale) << std::endl;
  }
  if (hChargeDensity == nullptr)
  {
    // we don't have charge information, so set everything to zeroes
    std::cout << "no charge data found.  Setting all charges to 0.0" << std::endl;
    charge->SetAll(0);  // otherwise set the array to all zeroes.
    return;
  }

  //   0     1     2     ...   n-1
  // first|   ..|  ..  |  .. |last
  float dphi, dr, dz;  // internal rep bin widths in each dimension.  Got too confusing to make these an array.
  dr = binWidth[0];
  dphi = binWidth[1];
  dz = binWidth[2];
  float dzhist, drhist;
  dzhist = dz / inputAxisScale;
  drhist = dr / inputAxisScale;

  float phimid, zmid;  // position of the center of each fixed-width array bin, in the input histogram units
  // note that since we computed density using the hist units, we must use those units for the volume term again here.
  int i[3];
  for (i[0] = 0; i[0] < nBins[0]; i[0]++)
  {  // r
    float rmid = (lowerBound[0] + (i[0] + 0.5) * dr) / inputAxisScale;
    // float rlow = (lowerBound[0] + dr * i[0]) / inputAxisScale;
    float histBinVolume = dzhist * dphi * rmid * drhist;  // volume of our output bin in histogram units.
    // note that since we have equal bin widths, the volume term depends only on r.
    //  Q_bin=Density_ion_interp*bin_volume*(coulomb/ion)
    float scaleFactor = histBinVolume * inputChargeScale;  // hence the total scale factor is the volume term times the charge scale factor
    for (i[1] = 0; i[1] < nBins[1]; i[1]++)
    {  // phi
      phimid = lowerBound[1] + (i[1] + 0.5) * dphi;
      for (i[2] = 0; i[2] < nBins[2]; i[2]++)
      {  // z
        zmid = (lowerBound[2] + (i[2] + 0.5) * dz) / inputAxisScale;
        float q = 0;
        if (CanInterpolateAt(phimid, rmid, zmid, hChargeDensity))
        {  // interpolate if we can
          if (false)
          {  // deep debug statements.
            std::cout << boost::str(boost::format("function said we could interpolate at (r,phi,z)=(%.2f, %.2f,%.2f), bounds are:") % rmid % phimid % zmid);
            std::cout << boost::str(boost::format("  r: %.2f < %.2f < %.2f") % hChargeDensity->GetYaxis()->GetXmin() % rmid % (hChargeDensity->GetYaxis()->GetXmax()));
            std::cout << boost::str(boost::format("  p: %.2f < %.2f < %.2f") % hChargeDensity->GetXaxis()->GetXmin() % phimid % hChargeDensity->GetXaxis()->GetXmax());
            std::cout << boost::str(boost::format("  z: %.2f < %.2f < %.2f") % hChargeDensity->GetZaxis()->GetXmin() % zmid % hChargeDensity->GetZaxis()->GetXmax()) << std::endl;
          }
          q = scaleFactor * hChargeDensity->Interpolate(phimid, rmid, zmid);
        }
        else
        {  // otherwise, just take the central value and assume it's flat.  Better than a zero.
          q = scaleFactor * hChargeDensity->GetBinContent(hChargeDensity->FindBin(phimid, rmid, zmid));
        }
        if (false)
        {  // deep debug statements.
          int global = hSourceCharge->FindBin(phimid, rmid, zmid);
          if (CanInterpolateAt(phimid, rmid, zmid, hChargeDensity))
          {
            std::cout << boost::str(boost::format("density debug report (interp) (r,phi,z)=(%.2f, %.2f,%.2f), glob=%d, q_dens=%E") % rmid % phimid % zmid % global % q);
            std::cout << boost::str(boost::format(", density=%E, vol=%E") % (hChargeDensity->Interpolate(phimid, rmid, zmid)) % histBinVolume);
            std::cout << boost::str(boost::format(", q_bin=%E, q_bin_coul=%E") % (hSourceCharge->GetBinContent(hSourceCharge->FindBin(phimid, rmid, zmid))) % (hSourceCharge->GetBinContent(hSourceCharge->FindBin(phimid, rmid, zmid)) * inputChargeScale));
            std::cout << boost::str(boost::format(", q_interp=%E, q_bin_coul/vol=%E") % (hSourceCharge->Interpolate(phimid, rmid, zmid)) % (hSourceCharge->GetBinContent(hSourceCharge->FindBin(phimid, rmid, zmid)) / histBinVolume)) << std::endl;
          }
          else
          {
            std::cout << boost::str(boost::format("density debug report (getbin) (r,phi,z)=(%.2f, %.2f,%.2f), glob=%d, q_dens=%E") % rmid % phimid % zmid % global % q);
            std::cout << boost::str(boost::format(", density=%E, vol=%E") % hChargeDensity->GetBinContent(hChargeDensity->FindBin(phimid, rmid, zmid)) % histBinVolume);
            std::cout << boost::str(boost::format(", q_bin=%E, q_bin_ions=%E") % hSourceCharge->GetBinContent(hSourceCharge->FindBin(phimid, rmid, zmid)) % (hSourceCharge->GetBinContent(hSourceCharge->FindBin(phimid, rmid, zmid)) / inputChargeScale));
            std::cout << boost::str(boost::format(", q_bin_ions/vol=%E") % (hSourceCharge->GetBinContent(hSourceCharge->FindBin(phimid, rmid, zmid)) / scaleFactor)) << std::endl;
          }
        }
        charge->Set(i[0], i[1], i[2], q);

      }  // z
    }    // phi
  }      // r

  if (DEBUG)
  {
    std::cout << "done regenerating charge array contents" << std::endl;
  }

  return;
}

void ChargeMapReader::RegenerateDensity()
{
  // assume the input map has changed, so we need to rebuild our internal representation of the density.
  // this is done by cloning the SourceCharge histogram and dividing each bin in it by its volume
  if (DEBUG)
  {
    std::cout << "regenerating density histogram" << std::endl;
  }

  // if we have one already, delete it.
  if (hChargeDensity != nullptr)
  {
    if (DEBUG)
    {
      std::cout << "deleting old density histogram" << std::endl;
    }
    delete hChargeDensity;
  }
  if (hSourceCharge == nullptr)
  {
    // the source data doesn't exist, so we will fail if we try to clone
    std::cout << "no source charge data file is open, or the histogram was not found." << std::endl;
    return;
  }

  // clone this from the source histogram, which we assume is open.
  hChargeDensity = static_cast<TH3*>(hSourceCharge->Clone("hChargeDensity"));
  hChargeDensity->Reset();

  // then go through it, bin by bin, and replace each bin content with the corresponding density, so we can interpolate correctly.
  // TODO:  Does this mean we once again need 'guard' bins?  Gross.

  TAxis* ax[3] = {nullptr, nullptr, nullptr};
  ax[0] = hChargeDensity->GetXaxis();
  ax[1] = hChargeDensity->GetYaxis();
  ax[2] = hChargeDensity->GetZaxis();

  int nbins[3];
  for (int i = 0; i < 3; i++)
  {
    nbins[i] = ax[i]->GetNbins();  // number of bins, not counting under and overflow.
  }

  //   0     1     2   ...   n-1    n    n+1
  // under|first|   ..|  ..  |  .. |last| over
  int i[3];
  float low[3], high[3];
  float dr, dz;  // bin widths in each dimension.  Got too confusing to make these an array.
  // note that all of this is done in the native units of the source data, specifically, the volume element is in hist units, not our internal units.

  for (i[0] = 1; i[0] <= nbins[0]; i[0]++)
  {  // phi
    int a = 0;
    low[a] = ax[a]->GetBinLowEdge(i[a]);
    high[a] = ax[a]->GetBinUpEdge(i[a]);
    float dphi = high[a] - low[a];
    for (i[1] = 1; i[1] <= nbins[1]; i[1]++)
    {  // r
      a = 1;
      low[a] = ax[a]->GetBinLowEdge(i[a]);
      high[a] = ax[a]->GetBinUpEdge(i[a]);
      dr = high[a] - low[a];
      float rphiterm = dphi * (low[1] + 0.5 * dr) * dr;
      for (i[2] = 1; i[2] <= nbins[2]; i[2]++)
      {  // z
        a = 2;
        low[a] = ax[a]->GetBinLowEdge(i[a]);
        high[a] = ax[a]->GetBinUpEdge(i[a]);
        dz = high[a] - low[a];
        // float volume=dz*dphi*(low[1]+0.5*dr)*dr;
        float volume = dz * rphiterm;
        int globalBin = hSourceCharge->GetBin(i[0], i[1], i[2]);
        float q = hSourceCharge->GetBinContent(globalBin);
        hChargeDensity->SetBinContent(globalBin, q / volume);
        if (false)
        {  // deep debug statements.
          std::cout << boost::str(boost::format("iprz=(%d,%d,%d),glob=%d") % i[0] % i[1] % i[2] % globalBin);
          std::cout << boost::str(boost::format("edges=[%.2f,%.2f],[%.1f,%.1f],[%.1f,%f.1],") % low[0] % high[0] % low[1] % high[1] % low[2] % high[2]);
          std::cout << boost::str(boost::format("\tq=%E,vol=%E,dens=%E") % q % volume % hChargeDensity->GetBinContent(globalBin)) << std::endl;
        }
      }
    }
  }
  if (DEBUG)
  {
    std::cout << "done regenerating density histogram" << std::endl;
  }

  return;
}

bool ChargeMapReader::ReadSourceCharge(const std::string& filename, const std::string& histname, float axisScale, float contentScale)
{
  // load the charge-per-bin data from the specified file.
  inputAxisScale = axisScale;
  inputChargeScale = contentScale;
  TFile* inputFile = TFile::Open(filename.c_str(), "READ");
  inputFile->GetObject(histname.c_str(), hSourceCharge);
  if (hSourceCharge == nullptr)
  {
    return false;
  }
  RegenerateDensity();
  RegenerateCharge();
  inputFile->Close();
  // after this, the source histogram doesn't exist anymore.
  return true;
}

bool ChargeMapReader::ReadSourceCharge(TH3* sourceHist, float axisScale, float contentScale)
{
  if (DEBUG)
  {
    std::cout << "reading charge from " << sourceHist->GetName() << std::endl;
  }
  inputAxisScale = axisScale;
  inputChargeScale = contentScale;
  hSourceCharge = sourceHist;  // note that this means we don't own this histogram!
  if (hSourceCharge == nullptr)
  {
    return false;
  }
  RegenerateDensity();
  RegenerateCharge();

  if (DEBUG)
  {
    std::cout << "done reading charge from " << sourceHist->GetName() << std::endl;
  }

  return true;
}

bool ChargeMapReader::ReadSourceAdc(const std::string& adcfilename, const std::string& adchistname, const std::string& ibfgainfilename, const std::string& ibfgainhistname, float axisScale, float contentScale)
{
  // load the adc-per-bin data from the specified file.
  TFile* adcInputFile = TFile::Open(adcfilename.c_str(), "READ");
  TH3* hAdc = nullptr;
  adcInputFile->GetObject(adchistname.c_str(), hAdc);
  if (hSourceCharge == nullptr)
  {
    std::cout << "Source Charge hist or file " << adcfilename << " not found!" << std::endl;
    return false;
  }
  // load the conversion factor from adc sum to ibf
  TFile* ibfGainInputFile = TFile::Open(ibfgainfilename.c_str(), "READ");
  TH2* hIbfGain = nullptr;
  ibfGainInputFile->GetObject(ibfgainhistname.c_str(), hIbfGain);
  if (hIbfGain == nullptr)
  {
    std::cout << "IBF Gain hist or file " << ibfgainfilename << " not found!" << std::endl;
    return false;
  }

  // and then hand the task off to the next function
  bool ret = ReadSourceAdc(hAdc, hIbfGain, axisScale, contentScale);
  adcInputFile->Close();
  // after this, the source histograms don't exist anymore.
  return ret;
}

bool ChargeMapReader::ReadSourceAdc(TH3* adcHist, TH2* gainHist, float axisScale, float contentScale)
{
  if (DEBUG)
  {
    std::cout << "reading ADCs from " << adcHist->GetName() << ", ibfGain from "
              << gainHist->GetName() << std::endl;
  }
  inputAxisScale = axisScale;
  inputChargeScale = contentScale;
  hSourceCharge = adcHist;  // note that this means we don't own this histogram!  It may disappear.
  if (hSourceCharge == nullptr)
  {
    return false;
  }

  // scale the voxels of the sourcecharge by the gainhist data.

  TAxis* ax[3] = {nullptr, nullptr, nullptr};
  ax[0] = hSourceCharge->GetXaxis();
  ax[1] = hSourceCharge->GetYaxis();
  ax[2] = hSourceCharge->GetZaxis();

  int nbins[3];
  for (int i = 0; i < 3; i++)
  {
    nbins[i] = ax[i]->GetNbins();  // number of bins, not counting under and overflow.
  }

  //   0     1     2   ...   n-1    n    n+1
  // under|first|   ..|  ..  |  .. |last| over
  int i[3];

  for (i[0] = 1; i[0] <= nbins[0]; i[0]++)
  {  // phi
    for (i[1] = 1; i[1] <= nbins[1]; i[1]++)
    {  // r
      int globalBin2D = gainHist->GetBin(i[0], i[1]);
      float scalefactor = gainHist->GetBinContent(globalBin2D);
      for (i[2] = 1; i[2] <= nbins[2]; i[2]++)
      {  // z
        int globalBin = hSourceCharge->GetBin(i[0], i[1], i[2]);
        float q = hSourceCharge->GetBinContent(globalBin);
        hSourceCharge->SetBinContent(globalBin, q * scalefactor);
        if (false)
        {  // deep debug statements.
          std::cout << boost::str(boost::format("applying gain to adc input:  iprz=(%d,%d,%d), glob3d=%d, glob2d=%d, adc=%f, scale=%f") % i[0] % i[1] % i[2] % globalBin % globalBin2D % q % scalefactor) << std::endl;
        }
      }  // z
    }    // r
  }      // phi
  if (DEBUG)
  {
    std::cout << "done converting adc hist to ions" << std::endl;
  }

  RegenerateDensity();
  RegenerateCharge();

  if (DEBUG)
  {
    std::cout << "done converting ADCs from " << adcHist->GetName() << std::endl;
  }

  return true;
}

bool ChargeMapReader::SetOutputParameters(int _nr, float _rmin, float _rmax, int _nphi, float _phimin, float _phimax, int _nz, float _zmin, float _zmax)
{
  // change all the parameters of our output array and rebuild the array from scratch.
  if (!(_rmax > _rmin) || !(_phimax > _phimin) || !(_zmax > _zmin))
  {
    return false;  // the bounds are not well-ordered.
  }
  if (_nr < 1 || _nphi < 1 || _nz < 1)
  {
    return false;  // must be at least one bin wide.
  }

  nBins[0] = _nr;
  nBins[1] = _nphi;
  nBins[2] = _nz;
  lowerBound[0] = _rmin;
  lowerBound[1] = _phimin;
  lowerBound[2] = _zmin;
  upperBound[0] = _rmax;
  upperBound[1] = _phimax;
  upperBound[2] = _zmax;

  for (int i = 0; i < 3; i++)
  {
    binWidth[i] = (upperBound[i] - lowerBound[i]) / (1.0 * nBins[i]);
  }

  // if the array exists, delete it.
  if (charge != nullptr)
  {
    if (DEBUG)
    {
      std::cout << "charge array existed.  deleting" << std::endl;
    }
    delete charge;
    charge = nullptr;
  }
  if (DEBUG)
  {
    std::cout << "building new charge array" << std::endl;
  }
  if (DEBUG)
  {
    std::cout << "should have " << (nBins[0] * nBins[1] * nBins[2])
              << " elements" << std::endl;
  }

  charge = new MultiArray<float>(nBins[0], nBins[1], nBins[2]);

  if (hChargeDensity != nullptr)
  {
    if (DEBUG)
    {
      std::cout << "charge density data exists, regenerating charge" << std::endl;
    }
    RegenerateCharge();  // fill the array with the charge data if available
  }
  else
  {
    charge->SetAll(0);  // otherwise set the array to all zeroes.
  }
  if (DEBUG)
  {
    std::cout << "finished building array" << std::endl;
  }

  return true;
}

bool ChargeMapReader::SetOutputBounds(float _rmin, float _rmax, float _phimin, float _phimax, float _zmin, float _zmax)
{
  // change all the bounds of our output array and rebuild the array from scratch, leaving the original binning

  if (!(_rmax > _rmin) || !(_phimax > _phimin) || !(_zmax > _zmin))
  {
    return false;  // the bounds are not well-ordered.
  }

  lowerBound[0] = _rmin;
  lowerBound[1] = _phimin;
  lowerBound[2] = _zmin;
  upperBound[0] = _rmax;
  upperBound[1] = _phimax;
  upperBound[2] = _zmax;

  for (int i = 0; i < 3; i++)
  {
    binWidth[i] = (upperBound[i] - lowerBound[i]) / (1.0 * nBins[i]);
  }

  // if the array exists, delete it.
  if (charge != nullptr)
  {
    delete charge;
    charge = nullptr;
  }
  charge = new MultiArray<float>(nBins[0], nBins[1], nBins[2]);

  if (hChargeDensity != nullptr)
  {
    RegenerateCharge();  // fill the array with the charge data if available
  }
  else
  {
    charge->SetAll(0);  // otherwise set the array to all zeroes.
  }
  return true;

  return true;
}

bool ChargeMapReader::SetOutputBins(int _nr, int _nphi, int _nz)
{
  // change the number of bins of our output array and rebuild the array from scratch, leaving the bounds alone.
  if (_nr < 1 || _nphi < 1 || _nz < 1)
  {
    return false;  // must be at least one bin wide.
  }
  nBins[0] = _nr;
  nBins[1] = _nphi;
  nBins[2] = _nz;

  for (int i = 0; i < 3; i++)
  {
    binWidth[i] = (upperBound[i] - lowerBound[i]) / (1.0 * nBins[i]);
  }

  // if the array exists, delete it.
  if (charge != nullptr)
  {
    delete charge;
    charge = nullptr;
  }
  charge = new MultiArray<float>(nBins[0], nBins[1], nBins[2]);

  if (hChargeDensity != nullptr)
  {
    RegenerateCharge();  // fill the array with the charge data if available
  }
  else
  {
    charge->SetAll(0);  // otherwise set the array to all zeroes.
  }
  return true;
}

void ChargeMapReader::AddChargeInBin(int r, int phi, int z, float q)
{
  assert(r > 0 && r < nBins[0]);
  assert(phi > 0 && phi < nBins[1]);
  assert(z > 0 && z < nBins[2]);
  if (DEBUG)
  {
    std::cout << boost::str(boost::format("adding charge in array element %d %d %d to %.2E") % r % phi % z % q) << std::endl;
  }

  charge->Add(r, phi, z, q);
  return;
}

void ChargeMapReader::AddChargeAtPosition(float r, float phi, float z, float q)
{
  // bounds checking are handled by the binwise function, so no need to do so here:
  AddChargeInBin((r - lowerBound[0]) / binWidth[0], (phi - lowerBound[1]) / binWidth[1], (z - lowerBound[2]) / binWidth[2], q);
  return;
}

float ChargeMapReader::GetChargeInBin(int r, int phi, int z)
{
  if (!(r >= 0 && r < nBins[0]))
  {
    std::cout << "requested rbin " << r << ", but bounds are 0 to "
              << nBins[0] << ". Failing." << std::endl;
    assert(r >= 0 && r < nBins[0]);
  }
  if (!(phi >= 0 && phi < nBins[1]))
  {
    std::cout << "requested phibin " << phi << ", but bounds are 0 to "
              << nBins[1] << ". Failing." << std::endl;
    assert(phi >= 0 && phi < nBins[1]);
  }
  if (!(z >= 0 && z < nBins[2]))
  {
    std::cout << "requested rbin " << z << ", but bounds are 0 to "
              << nBins[2] << ". Failing." << std::endl;
    assert(z >= 0 && z < nBins[2]);
  }

  if (DEBUG)
  {
    std::cout << "getting charge in array element " << r << " "
              << phi << " " << z << std::endl;
  }

  return charge->Get(r, phi, z);
}

float ChargeMapReader::GetChargeAtPosition(float r, float phi, float z)
{
  // bounds checking are handled by the binwise function, so no need to do so here:
  return GetChargeInBin((r - lowerBound[0]) / binWidth[0], (phi - lowerBound[1]) / binWidth[1], (z - lowerBound[2]) / binWidth[2]);
}

void ChargeMapReader::SetChargeInBin(int r, int phi, int z, float q)
{
  if (!(r >= 0 && r < nBins[0]))
  {
    std::cout << "requested rbin " << r << ", but bounds are to " << nBins[0]
              << ". Failing." << std::endl;
    assert(r >= 0 && r < nBins[0]);
  }
  if (!(phi >= 0 && phi < nBins[1]))
  {
    std::cout << "requested phibin " << phi << ", but bounds are 0 to "
              << nBins[1] << ". Failing." << std::endl;
    assert(phi >= 0 && phi < nBins[1]);
  }
  if (!(z >= 0 && z < nBins[2]))
  {
    std::cout << "requested rbin " << z << ", but bounds are 0 to "
              << nBins[2] << ". Failing." << std::endl;
    assert(z >= 0 && z < nBins[2]);
  }
  if (DEBUG)
  {
    std::cout << boost::str(boost::format("setting charge in array element %d %d %d to %.2E") % r % phi % z % q) << std::endl;
  }

  charge->Set(r, phi, z, q);
  return;
}

void ChargeMapReader::SetChargeAtPosition(float r, float phi, float z, float q)
{
  // bounds checking are handled by the binwise function, so no need to do so here:
  SetChargeInBin((r - lowerBound[0]) / binWidth[0], (phi - lowerBound[1]) / binWidth[1], (z - lowerBound[2]) / binWidth[2], q);
  return;
}
