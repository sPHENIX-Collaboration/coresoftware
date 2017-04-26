// TPC definitions
// holds the constants used in the TPC definition, simulation and reconstruction
// Author: Carlos Perez
#ifndef __TPCCONSTANTS_H__
#define __TPCCONSTANTS_H__

#include <TMath.h>
#include "TPCDataTypes.h"

using namespace TPCDataTypes;

namespace TPCConstants {
  const Float_t kTwoPi = TMath::TwoPi();
  const Float_t kInverseOfSqrt2 = 1.0/TMath::Sqrt(2);
  const Float_t kInverseOfSqrt12 = 1.0/TMath::Sqrt(12);
  const Int_t kNSectors = 12;
  const Int_t kNSections = 3;
  const Int_t kNModulesPerPlate = kNSectors*kNSections;
  const Int_t kNPadRowsPerModule = 16;
  const Int_t kNMaxPadCols = 2304;
  const Int_t kNColsPerMS[3] = { 96, 128, 192 }; // { 1152, 1536, 2304 } / 12;
  const Float_t kGasHalfOfLength = 105.5; // cm
  const Float_t kInverseOfDriftVelocity = 1000.0/6; // ns*cm-1 (6 cm/us)
  const Float_t kElectronsPerKeV = 30.0;
  const Float_t kInverseOfTimeBinWidth = 2.0/100; // ns-1 (10MSPS)
  const Float_t kTimeBinWidth = 1.0/kInverseOfTimeBinWidth; // ns
  const Float_t kcm2TimeBin = kInverseOfDriftVelocity*kInverseOfTimeBinWidth; // cm-1
  const Float_t kTimeBin2cm = 1.0/kcm2TimeBin; // cm
  const Float_t kTimeShapeRiseSquare = 40*40; // ns^s (120/3 ns)
  const Float_t kTimeShapeTailSquare = 100*100; // ns^2 (300/3 ns)
  const Float_t kTimeHalfOfLength = kGasHalfOfLength*kInverseOfDriftVelocity; // ns
  const Float_t kGasDiffusionTransverseSquare = 0.012*0.012; //cm^2
  const Float_t kGasDiffusionLongitudinalSquare = 0.012*0.012; //cm^2
  const Time_t kNumberOfTimeBins = static_cast<Time_t>(kGasHalfOfLength*kcm2TimeBin);
  const Float_t kGasInnerRadius = 21.0; // cm
  const Float_t kGasOuterRadius = 77.0; // cm
  const Float_t kModuleMargin = 0.2; // cm
  const Float_t kModuleSpacing = 0.1; // cm
  const Float_t kModuleStartRadius = kGasInnerRadius + kModuleMargin + kModuleSpacing; // cm
  const Float_t kModuleDeltaRadius = (kGasOuterRadius-kGasInnerRadius-2*kNSections*kModuleMargin-(kNSections+1)*kModuleSpacing)/kNSections; // cm
  const Float_t kPadRowStartRadius = kModuleStartRadius + kModuleMargin; // cm
  const Float_t kPadRowStep = kModuleDeltaRadius/kNPadRowsPerModule; // cm
  const Float_t kPadRowDeltaRadius = 0.9*kPadRowStep; // cm
  const Float_t kModuleSectionStep = kModuleDeltaRadius + 2*kModuleMargin + kModuleSpacing; // cm
  const Float_t kModuleSectorStep = kTwoPi/kNSectors; // rad
};

#endif
