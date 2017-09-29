#include <vector>
#include <map>
#include <utility>
#include <iostream>

#include <TMath.h>

#include "TPCDataTypes.h"
#include "TPCConstants.h"
#include "TPCPadMap.h"

using namespace TPCDataTypes;
using namespace TPCConstants;

//==========
TPCPadMap::TPCPadMap()
{
  //std::cout << "Constructing TPCPadMap... " << std::endl;
  GeneratePadsV1();
  //std::cout << "Constructing TPCPadMap... [DONE]" << std::endl;
}
//==========
void TPCPadMap::DumpModule(Module_t module)
{
  Module_t mod = GetReducedModule(module);
  std::cout << "=== DUMPING MODULE " << int(mod) << " ===" << std::endl;
  std::cout << " PAD  | Row | Col | X,Y | Rad,Phi | dW,dH" << std::endl;
  for(int i=0; i!=GetNumberOfPads(mod); ++i ) {
    Pad_t pad = Pad_t(i);
    int row = GetRow(mod,pad);
    int col = GetCol(mod,pad);
    PairOfFloats_t xy = GetXY(mod, pad);
    PairOfFloats_t rp = GetRP(mod, pad);
    PairOfFloats_t dx = GetDXDY(mod, pad);
    float xx = xy.first;
    float yy = xy.second;
    float ra = rp.first;
    float ph = rp.second;
    float dw = dx.first;
    float dh = dx.second;
    std::cout << pad << " | " << row << " " << col << " | ";
    std::cout << xx << "," << yy << " | " << ra << "," << ph << " | ";
    std::cout << dw << "," << dh << std::endl;
  }
}
//==========
void TPCPadMap::GeneratePadsV1()
{
  for(int mod=0; mod!=kNModulesPerPlate; ++mod) {
    fNPadRows[mod] = kNPadRowsPerModule;
    for(int prow=0; prow!=kNPadRowsPerModule; ++prow) {
      fNPadCols[mod].push_back( kNColsPerMS[GetSection(mod)] );
    }
  }
  //ofstream ofile("db_padmap.txt");
  float rini[3], rfin[3];
  rini[0] = kGasInnerRadius + kModuleMargin + kModuleSpacing;
  rini[1] = rini[0] + kModuleSectionStep;
  rini[2] = rini[1] + kModuleSectionStep;
  for(int mod=0; mod!=3; ++mod) rfin[mod] = rini[mod] + kModuleDeltaRadius;
  // creating map to global coordinates
  for(int section=0; section!=kNSections; ++section) { // looping on sections
    int prows = kNPadRowsPerModule;
    int pcols = kNColsPerMS[section];
    float dr_prow = ( rfin[section] - rini[section] ) / prows;
    for(int sector=0; sector!=kNSectors; ++sector) { // looping on sectors
      int module = section*kNSectors + sector;
      for(int prow=0; prow!=prows; ++prow) { // looping on padrows
	float rad = rini[section] + dr_prow * prow + dr_prow/2.0;
	float dphi = TMath::TwoPi()/kNSectors;
	float pini = (kModuleSpacing/2.0 + kModuleMargin)/rad + sector*dphi;
	float dp_pcol = (dphi-(kModuleSpacing+kModuleMargin*3)/rad) / pcols;
	for(int pcol=0; pcol!=pcols; ++pcol) { // looping on padcols
	  float phi = pini + dp_pcol * pcol + dp_pcol/2.0;
	  float dyy = dr_prow/2.0;
	  float dxx = (rad-dr_prow/2.0)*dp_pcol/2.0;
	  float xxx = rad*TMath::Cos(phi);
	  float yyy = rad*TMath::Sin(phi);
	  Pad_t gidx = pcols*prow+pcol;
	  fRP[module].insert( std::make_pair( gidx, std::make_pair(rad,phi) ) );
	  fXY[module].insert( std::make_pair( gidx, std::make_pair(xxx,yyy) ) );
	  fDXDY[module].insert( std::make_pair( gidx, std::make_pair(dxx,dyy) ) );
	  //ofile << module << " " <<  gidx << " " << xxx << " " << yyy;
	  //ofile << " " << rad << " " << phi << " " << dxx << " " << dyy << endl;
	}
      }
    }
  }
}
//==========
ModuleRange_t TPCPadMap::FindModuleRP(float rad, float phi, float rms, float ele)
{
  // returns list of potential modules that are
  // contained in box defined by 5*sigma contour
  // to gaussian distribution centred at rad, phi
  // please make sure rms is positive define
  float window = 5.0*rms;
  ModuleRange_t mRange;
  float rmin = rad - window - kModuleStartRadius; // reduced radius
  float rmax = rmin + 2*window;
  if(rmin<0) rmin = 0;
  if(rmax<0) rmax = 0;
  int scSta = rmin/kModuleSectionStep; // contains the margin and spacing 
  int scEnd = rmax/kModuleSectionStep; // which is useful for border cases
  float pmin = phi - window/rad;
  float pmax = phi + window/rad;
  if(pmin<0) pmin += TMath::TwoPi();
  if(pmax>TMath::TwoPi()) pmax -= TMath::TwoPi();
  int stSta = pmin/kModuleSectorStep; // contains the margin and spacing
  int stEnd = pmax/kModuleSectorStep; // which is useful for border cases
  // worst case scenario, though very unlikely,
  // the cloud is on a corner and touches four modules
  // we are woking on this assumption first
  mRange.push_back( static_cast<Module_t>(scSta*kNSectors+stSta) );
  mRange.push_back( static_cast<Module_t>(scSta*kNSectors+stEnd) );
  mRange.push_back( static_cast<Module_t>(scEnd*kNSectors+stSta) );
  mRange.push_back( static_cast<Module_t>(scEnd*kNSectors+stEnd) );
  // now we remove duplicates.
  sort( mRange.begin(), mRange.end() );
  mRange.erase( unique( mRange.begin(), mRange.end() ), mRange.end() );
  return mRange;
}
//==========
PadQuotaRange_t TPCPadMap::GetPadQuotasRP( float ele, Module_t mod, float rad, float phi, float rms )
{
  // return quota for pads that are contained
  // in box degined by 5*sigma contour to gaussian
  // distribution centred at rad, phi
  // pdf = 1/(2*pi*sigma^2) Exp(-(x-mux)^2/(2*sigma^2)) Exp(-(y-muy)^2/(2*sigma^2))
  // pdf = 1/(2*pi*sigma^2) Exp( -(r^2 + R0^2 - 2rR0*Cos(dphi)) / (2*sigma^2) )
  // integral in r: Exp( - (R0*Sin(dphi))^2 /(2*sigma^2) ) *
  //                (-Erf( (Rmin-R0*Cos(dphi))/(sqrt2*sigma) ) + Erf( (Rmax-R0*Cos(dphi))/(sqrt2*sigma) ) /
  //                (2*sqrt2*pi*sigma)
  // for dphi small: f * Gauss( rphi, rphi0, sigma )  ,
  // where f = 0.5 ( Erf((Rmax-R0)/(sqrt2*sigma)) + Erf((R0-Rmin)/(sqrt2*sigma)) )
  float window = 5.0*rms;
  PadQuotaRange_t qRange;
  PadQuota_t q;
  float start = GetSection(mod)*kModuleSectionStep;
  float rmin = rad - window - kModuleStartRadius - start; // reduced radius
  float rmax = rmin + 2.0*window;
  if(rmin<0) rmin = 0;
  if(rmax<0) rmax = 0;
  if(rmin>kModuleDeltaRadius-1e-2) rmin = kModuleDeltaRadius-1e-2; // subtracting 100um
  if(rmax>kModuleDeltaRadius-1e-2) rmax = kModuleDeltaRadius-1e-2; // prevents overflow
  Pad_t padRowSta = TMath::FloorNint((rmin/kModuleDeltaRadius)*kNPadRowsPerModule);
  Pad_t padRowEnd = TMath::FloorNint((rmax/kModuleDeltaRadius)*kNPadRowsPerModule);
  //std::cout << "padRowSta " << padRowSta << " padRowEnd " << padRowEnd << std::endl;
  //---
  float phistart = GetSector(mod)*kModuleSectorStep + kModuleSpacing/rad/2.0 + kModuleMargin/rad; // not accurate but fast
  //std::cout << "phi " << phi << "|| phistart " << phistart << std::endl;
  //std::cout << "window/rad " << window/rad << std::endl;
  float pmin = phi - window/rad - phistart;
  float pmax = phi + window/rad - phistart;
  if(pmin<0) pmin = 0;
  if(pmax<0) pmax = 0;
  float invSigmaSqrt2 = kInverseOfSqrt2/rms;
  //std::cout << "pmin " << pmin << " pmax " << pmax << std::endl;
  float moduleDeltaPhi = kModuleSectorStep - (2*kModuleMargin + kModuleSpacing)/rad;
  if(pmin>moduleDeltaPhi) pmin = moduleDeltaPhi;
  if(pmax>moduleDeltaPhi) pmax = moduleDeltaPhi;
  float arcmu = rad*phi;
  for(Pad_t row = padRowSta; row!=padRowEnd+1; ++row) {
    //std::cout << "row " << row << " ncols=" << GetNumberOfCols(mod,row) << std::endl;
    Pad_t padColSta = TMath::FloorNint(pmin*GetNumberOfCols(mod,row)/moduleDeltaPhi);
    Pad_t padColEnd = TMath::FloorNint(pmax*GetNumberOfCols(mod,row)/moduleDeltaPhi+1);
    //std::cout << "padColSta " << padColSta << " padColEnd " << padColEnd << std::endl;
    PairOfFloats_t wh = GetDXDY(mod,0);
    float halfWidth = wh.first;
    if(0) {
      halfWidth *= 1.5;
    }
    //std::cout << "halfWidth " << halfWidth << std::endl;
    for(Pad_t col = padColSta; col!=padColEnd+1; ++col) {
      q.first = GetPad(mod,row,col);
      PairOfFloats_t rp = GetRP(mod,q.first);
      float arc1 = rad*rp.second - halfWidth;
      float arc2 = rad*rp.second + halfWidth;
      //std::cout << " arc1 arc2 arcmu " << arc1 << " " << arc2 << " " << arcmu << std::endl;
      //std::cout << " || (arc1-arcmu)/rms (arc2-arcmu)/rms " << (arc1-arcmu)/rms << " " << (arc2-arcmu)/rms << std::endl;
      float at1 = TMath::Erf((arc1-arcmu)*invSigmaSqrt2);
      float at2 = TMath::Erf((arc2-arcmu)*invSigmaSqrt2);
      float rad1 = rp.first-wh.second - rad;
      float rad2 = rp.first+wh.second - rad;
      float rat1 = TMath::Erf(rad1*invSigmaSqrt2);
      float rat2 = TMath::Erf(rad2*invSigmaSqrt2);
      float factorr = 0.5*(rat2-rat1);
      q.second = factorr*0.5*(at2-at1)*ele; // integral of pdf from (1) to (2)
      if(q.second<10) continue; // skip small signals
      qRange.push_back( q );
    }
  }
  return qRange;
}
//==========
TimeQuotaRange_t TPCPadMap::GetTimeQuotas(float ele, float mean, float rms0, float rms1)
{
  // Asummes a gaussian distribution with sigma0
  // for rise and sigma1 for tail
  // open a 5 sigma window and distributes
  // the electrons in bins of time.
  // Returns range of quotas.
  // pdf = 1/(2*pi*sigma^2) Exp(-(x-mu)^2/(2*sigma^2))
  // cdf = 1/2 * ( 1 + erf( (x-mu) / (sqrt2*sigma) ) )
  TimeQuotaRange_t tRange;
  TimeQuota_t tQ;
  float myrange0 = 5.0*rms0;
  float myrange1 = 5.0*rms1;
  float invSigma0Sqrt2 = kInverseOfSqrt2/rms0;
  float invSigma1Sqrt2 = kInverseOfSqrt2/rms1;
  float start = mean - myrange0;
  float end = mean + myrange1;
  //std::cout << "   TIMEQUOTA  mean " <<  mean << " [" << start << "," << end << "]" << std::endl;
  if(start<0) start = 0;
  if(end>(kTimeHalfOfLength-1e-2)) start = kTimeHalfOfLength-1e-2; // minus 10 ps
  BinTime_t stbin = static_cast<BinTime_t>(start*kInverseOfTimeBinWidth);
  BinTime_t enbin = static_cast<BinTime_t>(end*kInverseOfTimeBinWidth);
  //std::cout << "              bins [" <<  int(stbin) << "," << int(enbin) << "]" << std::endl;
  for(int i=stbin; i!=(enbin+1); ++i) { // opens time interval
    float x1 = i*kTimeBinWidth; // cm begin of time bin
    float x2 = x1 + kTimeBinWidth; // cm end of time bin
    float dx1 = x1-mean;
    float dx2 = x2-mean;
    float at1, at2;
    if(dx1<0) at1 = TMath::Erf(dx1*invSigma0Sqrt2);
    else at1 = TMath::Erf(dx1*invSigma1Sqrt2);
    if(dx2<0) at2 = TMath::Erf(dx2*invSigma0Sqrt2);
    else at2 = TMath::Erf(dx2*invSigma1Sqrt2);
    tQ.second = static_cast<Adc_t> (0.5*(at2-at1)*ele); // integral of pdf from (1) to (2)
    if(tQ.second<1) continue; // skip empty bins
    tQ.first = static_cast<BinTime_t> (i);
    tRange.push_back( tQ );
  }
  return tRange;
}
//==========
int TPCPadMap::GetNumberOfModules() {
  return kNModulesPerPlate*2;
}
int TPCPadMap::GetReducedModule(int mod)
{
  return mod%kNModulesPerPlate;
}
int TPCPadMap::GetNumberOfPads(Module_t mod)
{
  return int(fDXDY[GetReducedModule(mod)].size());
}
int TPCPadMap::GetNumberOfRows(Module_t mod)
{
  return fNPadRows[GetReducedModule(mod)];
}
int TPCPadMap::GetNumberOfCols(Module_t mod, Pad_t prow)
{
  return fNPadCols[GetReducedModule(mod)].at(prow);
}
int TPCPadMap::GetSection(Module_t mod)
{
  return GetReducedModule(mod)/int(kNSectors);
}
int TPCPadMap::GetSector(Module_t mod)
{
  return GetReducedModule(mod)%int(kNSectors);
}
Pad_t TPCPadMap::GetPad(Module_t mod, Pad_t prow, Pad_t pcol)
{
  return GetNumberOfCols(mod,prow)*prow+pcol;
}
Pad_t TPCPadMap::GetPad(Module_t mod, float rad, float phi)
{
  return 0;
}
PairOfFloats_t TPCPadMap::GetXY(Module_t mod, Pad_t pad)
{
  return fXY[GetReducedModule(mod)][pad];
}
PairOfFloats_t TPCPadMap::GetRP(Module_t mod, Pad_t pad)
{
  return fRP[GetReducedModule(mod)][pad];
}
PairOfFloats_t TPCPadMap::GetDXDY(Module_t mod, Pad_t pad)
{
  return fDXDY[GetReducedModule(mod)][pad];
}
PairOfFloats_t TPCPadMap::GetXY(Module_t mod, Pad_t prow, Pad_t pcol)
{
  return GetXY(mod,GetPad(mod,prow,pcol));
}
PairOfFloats_t TPCPadMap::GetRP(Module_t mod, Pad_t prow, Pad_t pcol)
{
  return GetRP(mod,GetPad(mod,prow,pcol));
}
PairOfFloats_t TPCPadMap::GetDXDY(Module_t mod, Pad_t prow, Pad_t pcol)
{
  return GetDXDY(mod,GetPad(mod,prow,pcol));
}
int TPCPadMap::GetRow(Module_t mod, Pad_t pad) {
  return static_cast<int>(pad/kNColsPerMS[GetSection(mod)]);
}
int TPCPadMap::GetCol(Module_t mod, Pad_t pad) {
  return static_cast<int>(pad%kNColsPerMS[GetSection(mod)]);
}
