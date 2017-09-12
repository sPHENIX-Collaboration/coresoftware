// holds the constants used in the TPC definition, simulation and reconstruction
// Author: Carlos Perez
#ifndef __TPCCONSTANTS_H__
#define __TPCCONSTANTS_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include "TPCDataTypes.h"

using namespace TPCDataTypes;

namespace TPCConstants {
  // math
  const float kPi = 4.0*std::atan(1.0);
  const float kTwoPi = 2.0*kPi;
  const float kInverseOfSqrt2 = 1.0/std::sqrt(2.0);
  const float kInverseOfSqrt12 = 1.0/std::sqrt(12.0);

  // geometry
  const float kGasHalfOfLength = 105.5; // cm
  const float kGasInnerRadius = 21.0; // cm
  const float kGasOuterRadius = 77.0; // cm

  // nominal fields
  const float kNominalEField = 400; // V/cm
  const float kNominalBField = 1.4; // Tesla

  // gas properties
  const float kInverseOfDriftVelocity = 1000.0/6; // ns*cm-1 (6 cm/us)
  const float kGasIonMobility = 4; // cm^2/(V.s)
  const float kGasDiffusionTransverseSquare = 0.012*0.012; //cm^2
  const float kGasDiffusionLongitudinalSquare = 0.012*0.012; //cm^2
  const float kElectronsPerKeV = 30.0;
  const int   kMaxElectronsPerHit = 3;

  // amplification
  const float kAmplificationNominal = 2000;
  const float kAmplificationSmearing = 0.01; //cm

  // padplane
  const int   kNSectors = 12;
  const int   kNSections = 3;
  const int   kNModulesPerPlate = kNSectors*kNSections;
  const int   kNPadRowsPerModule = 16;
  const int   kNMaxPadCols = 2304;
  const int   kNColsPerMS[3] = { 96, 128, 192 }; // { 1152, 1536, 2304 } / 12;
  const float kModuleMargin = 0.2; // cm
  const float kModuleSpacing = 0.1; // cm
  const float kModuleStartRadius = kGasInnerRadius + kModuleMargin + kModuleSpacing; // cm
  const float kModuleDeltaRadius = (kGasOuterRadius-kGasInnerRadius-2*kNSections*kModuleMargin-(kNSections+1)*kModuleSpacing)/kNSections; // cm
  const float kPadRowStartRadius = kModuleStartRadius + kModuleMargin; // cm
  const float kPadRowStep = kModuleDeltaRadius/kNPadRowsPerModule; // cm
  const float kPadRowDeltaRadius = 0.9*kPadRowStep; // cm
  const float kModuleSectionStep = kModuleDeltaRadius + 2*kModuleMargin + kModuleSpacing; // cm
  const float kModuleSectorStep = kTwoPi/kNSectors; // rad

  // electronics
  const float kInverseOfTimeBinWidth = 1.0/100; // ns-1 (10MSPS)
  const float kTimeBinWidth = 1.0/kInverseOfTimeBinWidth; // ns
  const float kcm2TimeBin = kInverseOfDriftVelocity*kInverseOfTimeBinWidth; // cm-1
  const float kTimeBin2cm = 1.0/kcm2TimeBin; // cm
  const float kTimeShapeRiseSquare = 60*60; // ns^s (60 ns)
  const float kTimeShapeTailSquare = 130*130; // ns^2 (130 ns)
  const float kTimeHalfOfLength = kGasHalfOfLength*kInverseOfDriftVelocity; // ns
  const BinTime_t kNumberOfTimeBins = static_cast<BinTime_t>(kGasHalfOfLength*kcm2TimeBin);

};

#endif
