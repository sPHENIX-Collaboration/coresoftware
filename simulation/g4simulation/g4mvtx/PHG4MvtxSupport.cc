#include "PHG4MvtxSupport.h"

#include "PHG4MvtxCable.h"
#include "PHG4MvtxDisplayAction.h"
#include "PHG4MvtxServiceStructure.h"

#include <g4main/PHG4Detector.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4Polycone.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4UnionSolid.hh>

#include <TString.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/format.hpp>
#pragma GCC diagnostic pop

#include <algorithm>  // for max
#include <cmath>      // for M_PI, atan, cos, sin
#include <cstddef>    // for NULL
#include <iostream>   // for operator<<, basic_...
#include <utility>    // for pair, make_pair
#include <regex>      // for regex

class G4VSolid;

namespace ServiceProperties
{
//  std::string materials[] = {"G4_Cu", "MVTX_CarbonFiber$", "G4_POLYETHYLENE"};
//  const int nMaterials = sizeof(materials) / sizeof(materials[0]);

  const double sEndWheelSNHolesZdist = 308.0 * mm;  // ALIITSUP0187,0177,0143
  const double sEndWStepHoleZpos = 4. * mm;
  const double sEndWStepHoleZdist = 4. * mm;
  const double sEndWheelNLen = 30. * mm;

  const double sCYSSFlgSsfFlgNsf = 606.59 * mm; // sf -> south face
  const double sCYSSFlgSsfCylsf = 181. * mm;
  const double sCYSSFlgSsfConesf = 7. * mm;
  const double sCYSSFlgSsfRibsf = 50. * mm;
  const double sPP2sfSBsf = 791.77 * mm;
  const double sPP2Len = 80 * mm;


  float ServiceEnd = -7.21 * cm ;
  float ServiceOffset = -16.0 * cm;
  float BarrelOffset = 18.679 * cm;
  float BarrelRadius = 10.33 * cm;     //Inner radious of service barrel
  float BarrelThickness = 0.436 * cm;  //Thickness in cm
  float BarrelLength = 1214.4 * mm;    //Length of cylinder in cm
  float BarrelCableStart = sEndWheelSNHolesZdist / 2 - ( sEndWStepHoleZpos + sEndWStepHoleZdist ) \
                         + sEndWheelNLen - sCYSSFlgSsfFlgNsf - BarrelLength + sPP2sfSBsf + sPP2Len / 2;
  float BarrelCableEnd = BarrelCableStart - ( sPP2sfSBsf + sPP2Len / 2 );
  float LayerThickness = 0.1 * cm;  //
  float CYSSConeThickness = 0.216 * cm;
  float CYSSRibThickness = 0.170 * cm;
  float cableRotate[3] = {10., 5., 5.};  //Rotate the cables to line up with the staves in deg
}  // namespace ServiceProperties

using namespace ServiceProperties;

//________________________________________________________________________________
PHG4MvtxSupport::PHG4MvtxSupport( PHG4MvtxDisplayAction *dispAct, bool overlapCheck )
  : m_DisplayAction(dispAct)
  , m_avSupport(nullptr)
  , m_avBarrelCable(nullptr)
{
  m_avLayerCable = { nullptr, nullptr, nullptr };
  m_overlapCheck = overlapCheck;
}

//________________________________________________________________________________
PHG4MvtxSupport::~PHG4MvtxSupport()
{
  if ( m_avSupport ) delete m_avSupport;
  if ( m_avBarrelCable ) delete m_avBarrelCable;
  for ( auto av : m_avLayerCable )
  {
    if ( av ) delete av;
  }
}

//________________________________________________________________________________
void PHG4MvtxSupport::CreateMvtxSupportMaterials()
{
  G4double density;
  G4int natoms;
  G4String symbol;
  G4String material = "MVTX_EW_Al$";
  if ( ! PHG4Detector::GetDetectorMaterial( material.c_str(), false ) )
  {
    auto MVTX_Al = new G4Material( material.c_str(), density = 2.7 * g / cm3,
                                   natoms = 1, kStateSolid );
    MVTX_Al->AddMaterial( PHG4Detector::GetDetectorMaterial( "G4_Al" ), 100 * perCent );
  }

  material = "MVTX_Epoxy$";
  if ( ! PHG4Detector::GetDetectorMaterial( material.c_str(), false ) )
  {
    auto MVTX_Epoxy = new G4Material( material.c_str(), density = 1.56 * g / cm3, natoms = 4 );
    MVTX_Epoxy->AddElement( PHG4Detector::GetDetectorElement( "H", true ), 32 );  // Hydrogen
    MVTX_Epoxy->AddElement( PHG4Detector::GetDetectorElement( "C", true ), 2 );   // Nitrogen
    MVTX_Epoxy->AddElement( PHG4Detector::GetDetectorElement( "N", true ), 4 );   // Oxygen
    MVTX_Epoxy->AddElement( PHG4Detector::GetDetectorElement( "O", true ), 15 );  // Carbon
  }

  material = "MVTX_CarbonFiber$";
  if ( ! PHG4Detector::GetDetectorMaterial( material.c_str(), false ) )
  {
    auto MVTX_CF = new G4Material( material.c_str(), density = 1.987 * g / cm3,
                                   natoms = 2 );
    MVTX_CF->AddMaterial( PHG4Detector::GetDetectorMaterial( "G4_C" ), 70 * perCent ); // Carbon (NX-80-240)
    MVTX_CF->AddMaterial( PHG4Detector::GetDetectorMaterial( "MVTX_Epoxy$" ), 30 * perCent ); // Epoxy (EX-1515)
  }
}

//________________________________________________________________________________
void PHG4MvtxSupport::CreateEndWheelsSideN( G4AssemblyVolume *& av )
{
  for ( unsigned int iLay = 0; iLay < PHG4MvtxDefs::kNLayers; iLay++ )
  {
    GetEndWheelSideN( iLay, av );
  }
}

//________________________________________________________________________________
void PHG4MvtxSupport::GetEndWheelSideN( const int lay, G4AssemblyVolume *&endWheel )
{
  //
  // Creates the single End Wheel on Side North
  // for a given layer of the MVTX detector
  // (Layer 0: MVTX-2-S-00001)
  // (Layer 1: MVTX-2-S-00016)
  // (Layer 2: MVTX-2-S-00023)
  //
  // Input:
  //         iLay : the layer number
  //         endWheel : the whole end wheel volume
  //                    where to place the current created wheel
  // Output:
  //
  // Return:
  //
  // Created:      10 Jan 2023 Yasser Corrales Morales
  //               (partially based on M. Sitta implementation in ITS2)
  //

  // The Basis Cone North Side is two halves,
  // but for sake of simplicity, so here they are made as a single cone.
  const double sEndWheelNRmax[3] = { 30.57 * mm, 38.59 * mm, 46.30 * mm };
  const double sEndWheelNRmin[3] = { 22.83 * mm, 30.67 * mm, 38.70 * mm };
  const double sEndWheelNWallR[3] = { 29.57 * mm, 37.59 * mm, 45.30 * mm };
  const double sEndWheelNThick = 2.8 * mm;

  const int sEndWNBaseNBigHoles = 6;
  const int sEndWNBaseNSmallHoles = 5;
  const double sEndWNBaseBigHoleD = 4 * mm;
  const double sEndWNBaseSmallHoleD = 1.6 * mm;
  const double sEndWNBaseHolesR[3] = { 27. * mm, 35. * mm, 42.7 * mm };
  const double sEndWNBaseHolesPhi = 15.; // Deg

  const int    sEndWNWallNHoles[3] = { 6, 8, 10 };
  const double sEndWNWallHoleD = 4.5 * mm;

  // The End Wheel Steps
  const double sEndWNStepXlen[3] = { 9.88 * mm, 11.2 * mm, 8.38 * mm };
  const double sEndWNStepYdispl[3] = { 26.521 * mm, 33.891 * mm, 41.64 * mm };
  const double sEndWNStepHolePhi[3] = { 30.0, 22.5, 18.0 }; // Deg
  const double sEndWNStepHoleTilt[3] = { 0., 0.295, 0.297 };  // Rad

  const double sEndWNStepZlen = sEndWheelNLen - sEndWheelNThick;

  // local varianbles
  double rmin, rmax, phimin, dphi;
  float xpos, ypos, zpos, rpos;
  float xlen, ylen;

  auto Ta = G4ThreeVector();
  auto Ra = G4RotationMatrix();

  const unsigned short nZplanes   = 4;
  const double zPlane[nZplanes] = { 0. * mm,
                                    sEndWheelNLen-sEndWheelNThick,
                                    sEndWheelNLen-sEndWheelNThick,
                                    sEndWheelNLen };

  const double rInner[nZplanes] = { sEndWheelNWallR[lay], sEndWheelNWallR[lay],
                                    sEndWheelNRmin[lay], sEndWheelNRmin[lay] };

  const double rOuter[nZplanes] = { sEndWheelNRmax[lay], sEndWheelNRmax[lay],
                                    sEndWheelNRmax[lay], sEndWheelNRmax[lay] };

  G4VSolid *endWNBasis = new G4Polycone( Form("endwnbasis%d", lay), 0., 2. * M_PI * rad,
                                         nZplanes, zPlane, rInner, rOuter );

  // The holes in the veritcal wall
  auto endwnwalhol = new G4Tubs( Form( "endwnwalhol%d", lay ), 0, sEndWNWallHoleD / 2,
                                 1.5 * ( sEndWheelNRmax[lay] - sEndWheelNWallR[lay] ),
                                 0, 2 * M_PI * rad );

  rmin = ( sEndWheelNRmax[lay] + sEndWheelNWallR[lay] ) / 2;
  zpos = sEndWStepHoleZpos;
  dphi = 180. / sEndWNWallNHoles[lay];
  phimin = dphi / 2;
  for ( int ihole = 0; ihole < 2 * sEndWNWallNHoles[lay]; ihole++ )
  {
    double phi = phimin + ihole * dphi;
    xpos = rmin * sin( phi * deg );
    ypos = rmin * cos( phi * deg );
    Ra.set( -phi * deg, 90 * deg, 0. );
    Ta.set( xpos, ypos, zpos );
    endWNBasis = new G4SubtractionSolid( Form( "endwnbasis-endwnwalhol%dl%d", ihole, lay ),
                                         endWNBasis, endwnwalhol, &Ra, Ta );
  }

  // The holes in the base
  auto endwnbasBhol = new G4Tubs( Form( "endwnbasBhol%d", lay ), 0, sEndWNBaseBigHoleD / 2,
                                  1.5 * sEndWheelNThick, 0, 2 * M_PI * rad );

  auto endwnbasShol = new G4Tubs( Form( "endwnbasShol%d", lay ), 0, sEndWNBaseSmallHoleD / 2,
                                  1.5 * sEndWheelNThick, 0, 2 * M_PI * rad );

  rmin = sEndWNBaseHolesR[lay];
  zpos = sEndWheelNLen - ( sEndWheelNThick / 2 );
  phimin = 0.;
  for ( int ihole = 0; ihole < ( sEndWNBaseNBigHoles + sEndWNBaseNSmallHoles ); ihole++ )
  {
    phimin += sEndWNBaseHolesPhi;
    xpos = rmin * cos( phimin * deg );
    ypos = rmin * sin( phimin * deg );
    Ra.set( 0., 0., 0. );
    Ta.set( xpos, ypos, zpos );
    auto endwnbashol =  ( ihole < 2 || ihole == 5 || ihole > 8 ) ? endwnbasShol : endwnbasBhol;
    endWNBasis = new G4SubtractionSolid( Form( "endwnbasis-endwnwalhol%dl%da", ihole, lay ),
                                         endWNBasis, endwnbashol, &Ra, Ta );
    Ta = G4ThreeVector( -xpos, -ypos, zpos );
    endWNBasis = new G4SubtractionSolid( Form( "endwnbasis-endwnwalhol%dl%db", ihole, lay ),
                                         endWNBasis, endwnbashol, &Ra, Ta );
  }

  // Now the Step as a Composite Shape (subtraction of a Pcon from a BBox)
  // (cutting volume should be slightly larger than desired region)
  rmin = sEndWheelNWallR[lay];
  xlen = sEndWNStepXlen[lay] - ( ( lay == 1 ) ? 0.01 * mm : 0. * mm );
  double xDispl = sqrt( rmin * rmin - sEndWNStepYdispl[lay] * sEndWNStepYdispl[lay] ) - xlen;
  ylen = sqrt( rmin * rmin - xDispl * xDispl ) - sEndWNStepYdispl[lay];

  auto stepBoxNSh = new G4Box( Form("stepBoxNSh%d", lay), xlen / 2, ylen / 2,
                               ( sEndWNStepZlen / 2 ) - 0.05 * mm );

  xpos = xDispl + stepBoxNSh->GetXHalfLength();
  ypos = sEndWNStepYdispl[lay] + stepBoxNSh->GetYHalfLength();
  rpos = sqrt( xpos * xpos + ypos * ypos );

  phimin = ( 90. * deg ) - acos( sEndWNStepYdispl[lay] / rmin ) * rad - 5 * deg;
  dphi = ( 90. * deg ) - asin( xDispl / rmin ) * rad - phimin + 5 * deg;
  rmax = rmin + 2 * stepBoxNSh->GetYHalfLength();

  auto stepTubNSh = new G4Tubs( Form( "stepTubNSh%d", lay ), rmin, rmax,
                                1.5 * sEndWNStepZlen / 2, phimin, dphi );
  Ra.set( 0., 0., 0. );
  Ta.set( -xpos, -ypos, 0 );
  auto stepNSh = new G4SubtractionSolid( Form( "stepNSh%d", lay ),
                                         stepBoxNSh, stepTubNSh, &Ra, Ta );

  auto matAl = PHG4Detector::GetDetectorMaterial( "MVTX_EW_Al$" );

  auto endWheelNvol = new G4LogicalVolume( endWNBasis, matAl, Form( "EndWheelNBasis%d", lay ) );
  m_DisplayAction->AddVolume( endWheelNvol, "red" );

  auto stepNShLogVol = new G4LogicalVolume( stepNSh , matAl, Form( "StepNL%d_LOGIC", lay ) );
  m_DisplayAction->AddVolume( stepNShLogVol, "red" );

  // Finally put everything in the mother volume
  zpos = ( sEndWheelSNHolesZdist / 2 )  - ( sEndWStepHoleZpos + sEndWStepHoleZdist );
  Ra.set( 0., 0., 0. );
  Ta.set( 0, 0, zpos );
  endWheel->AddPlacedVolume( endWheelNvol, Ta, &Ra );

  dphi = atan( ypos/xpos );
  double phi0 = PHG4MvtxDefs::mvtxdat[lay][PHG4MvtxDefs::kPhi0];
  const int numberOfStaves = PHG4MvtxDefs::mvtxdat[lay][PHG4MvtxDefs::kNStave];
  zpos += ( stepBoxNSh->GetZHalfLength() );
  for ( int j = 0; j < numberOfStaves; j++ )
  {
    double phi = phi0 * rad;
    phi += j * sEndWNStepHolePhi[lay] * deg;
    phi -= sEndWNStepHoleTilt[lay];
    xpos = rpos * cos( phi );
    ypos = rpos * sin( phi );
    Ra.set( 0., 0., dphi - phi );
    Ta.set( xpos, ypos, zpos );
    endWheel->AddPlacedVolume( stepNShLogVol, Ta, &Ra );
  }

  return;
}

//________________________________________________________________________________
void PHG4MvtxSupport::CreateEndWheelsSideS( G4AssemblyVolume *&av )
{
  for ( unsigned int iLay = 0; iLay < PHG4MvtxDefs::kNLayers; iLay++ )
  {
    GetEndWheelSideS( iLay, av );
  }
}

//________________________________________________________________________________
void PHG4MvtxSupport::GetEndWheelSideS( const int lay, G4AssemblyVolume *&endWheel )
{
  //
  // Creates the single End Wheel on Side South
  // for a given layer of the MVTX detector
  // (Layer 0: MVTX-2-S-00002)
  // (Layer 1: MVTX-2-S-00017)
  // (Layer 2: MVTX-2-S-00024)
  //
  // Input:
  //         iLay : the layer number
  //         endWheel : the whole end wheel volume
  //                    where to place the current created wheel
  // Output:
  //
  // Return:
  //
  // Created:      17 Jan 2023 Yasser Corrales Morales
  //               (partially based on M. Sitta implementation in ITS2)
  //

  // The Basis Cone South Side is two halves,
  // but for sake of simplicity, so here they are made as a single cone.
  const double sEWSExtSectRmax[3] = { 30.97 * mm, 38.99 * mm, 46.74 * mm };
  const double sEWSExtSectThick[3] = { 0.97 * mm, 1.20 * mm, 1.01 * mm };

  const double sEWSIntSectRmax[3] = { 29.77 * mm, 37.79 * mm, 45.54 * mm };
  const double sEWSIntSectRmin[3] = { 28.80 * mm, 36.80 * mm, 44.50 * mm };
  const double sEWSIntSectLongLen = 26 * mm;
  const double sEWSIntSectShortLen = 25 * mm;

  const double sEWSTotalLen = 42. * mm;

  const int    sEndWSWallNHoles[3] = { 6, 8, 10 };
  const double sEndWSWallHoleD = 4.5 * mm;

  const double sEndWSStepXlen[3] = { 10.80 * mm, 10.05 * mm, 9.45 * mm };
  const double sEndWSStepYdispl[3] = { 26.52 * mm, 33.89 * mm, 41.64 * mm };
  const double sEndWSStepHolePhi[3] = { 30.0, 22.5, 18.0 }; // Deg
  const double sEndWSStepHoleTilt[3] = {0., 0.295, 0.297};  // Rad

  const double sEndWSStepZlen = sEWSTotalLen - sEWSIntSectLongLen;
  const double sEndWheelNWallR = sEWSExtSectRmax[lay] - sEWSExtSectThick[lay] + 0.05 * mm;

  // local varianbles
  double rmin, rmax, phimin, dphi;
  float xpos, ypos, zpos, rpos;
  float xlen, ylen;

  auto Ta = G4ThreeVector();
  auto Ra = G4RotationMatrix();

  const unsigned short nZplanes = 6;
  const double zPlane[nZplanes] = { 0. * mm,
                                    sEWSTotalLen - sEWSIntSectLongLen,
                                    sEWSTotalLen - sEWSIntSectLongLen,
                                    sEWSTotalLen - sEWSIntSectShortLen,
                                    sEWSTotalLen - sEWSIntSectShortLen,
                                    sEWSTotalLen };

  const double rInner[nZplanes] = { sEndWheelNWallR,
                                    sEndWheelNWallR,
                                    sEWSIntSectRmin[lay],
                                    sEWSIntSectRmin[lay],
                                    sEWSIntSectRmin[lay],
                                    sEWSIntSectRmin[lay] };

  const double rOuter[nZplanes] = { sEWSExtSectRmax[lay],
                                    sEWSExtSectRmax[lay],
                                    sEWSExtSectRmax[lay],
                                    sEWSExtSectRmax[lay],
                                    sEWSIntSectRmax[lay],
                                    sEWSIntSectRmax[lay] };

  G4VSolid *endWSBasis = new G4Polycone( Form( "endwsbasis%d", lay ), 0., 2. * M_PI * rad,
                                         nZplanes, zPlane, rInner, rOuter);

  // The holes in the veritcal wall
  auto endwswalhol = new G4Tubs( Form( "endwswalhol%d", lay ), 0, sEndWSWallHoleD / 2,
                                 1.5 * sEWSExtSectThick[lay],
                                 0, 2 * M_PI * rad);

  rmin = sEWSExtSectRmax[lay] - sEWSExtSectThick[lay] / 2;
  zpos = sEndWStepHoleZpos;
  dphi = 180. / sEndWSWallNHoles[lay];
  phimin = dphi / 2;
  for ( int ihole = 0; ihole < 2 * sEndWSWallNHoles[lay]; ihole++ )
  {
    double phi = phimin + ihole * dphi;
    xpos = rmin * sin( phi * deg );
    ypos = rmin * cos( phi * deg );
    Ra.set( -phi * deg, 90 * deg, 0 );
    Ta.set( xpos, ypos, zpos );
    endWSBasis = new G4SubtractionSolid( Form( "endwsbasis-endwswalhol%dl%d", ihole, lay ),
                                         endWSBasis, endwswalhol, &Ra, Ta );
  }

  // Now the Step as a Composite Shape (subtraction of a Pcon from a BBox)
  // (cutting volume should be slightly larger than desired region)
  rmin = sEWSExtSectRmax[lay] - sEWSExtSectThick[lay];
  xlen = sEndWSStepXlen[lay];
  double xDispl = sqrt( rmin * rmin - sEndWSStepYdispl[lay] * sEndWSStepYdispl[lay] ) - xlen;
  ylen = sqrt( rmin * rmin - xDispl * xDispl ) - sEndWSStepYdispl[lay];

  auto stepBoxSSh = new G4Box( Form( "stepBoxSSh%d", lay ), xlen / 2, ylen / 2,
                               ( sEndWSStepZlen / 2 ) - 0.05 * mm );

  xpos = xDispl + stepBoxSSh->GetXHalfLength();
  ypos = sEndWSStepYdispl[lay] + stepBoxSSh->GetYHalfLength();
  rpos = sqrt( xpos * xpos + ypos * ypos );

  phimin = ( 90. * deg ) - acos( sEndWSStepYdispl[lay] / rmin ) * rad - 5 * deg;
  dphi = ( 90. * deg ) - asin( xDispl / rmin ) * rad - phimin + 5 * deg;
  rmax = rmin + 2 * stepBoxSSh->GetYHalfLength();

  auto stepTubSSh = new G4Tubs( Form( "stepTubSSh%d", lay ), rmin, rmax,
                                1.5 * sEndWSStepZlen / 2, phimin, dphi );
  Ra.set( 0., 0., 0. );
  Ta.set( -xpos, -ypos, 0. );
  auto stepSSh = new G4SubtractionSolid( Form( "stepSSh%d", lay ),
                                         stepBoxSSh, stepTubSSh, &Ra, Ta );

  auto matAl = PHG4Detector::GetDetectorMaterial( "MVTX_EW_Al$" );

  auto endWheelSvol = new G4LogicalVolume( endWSBasis, matAl, Form( "EndWheelSBasis%d", lay ) );
  m_DisplayAction->AddVolume( endWheelSvol, "red" );

  auto stepSShLogVol = new G4LogicalVolume( stepSSh , matAl, Form( "StepSL%d_LOGIC", lay ) );
  m_DisplayAction->AddVolume( stepSShLogVol, "red" );

  // Finally put everything in the mother volume
  zpos = ( sEndWheelSNHolesZdist / 2 )  - ( sEndWStepHoleZpos + sEndWStepHoleZdist );
  Ra.set( 0., M_PI * rad, 0. );
  Ta.set( 0, 0, -zpos );
  endWheel->AddPlacedVolume( endWheelSvol, Ta, &Ra );

  dphi = atan( ypos/xpos );
  double phi0 = PHG4MvtxDefs::mvtxdat[lay][PHG4MvtxDefs::kPhi0];
  const int numberOfStaves = PHG4MvtxDefs::mvtxdat[lay][PHG4MvtxDefs::kNStave];
  zpos += ( stepBoxSSh->GetZHalfLength() );
  for ( int j = 0; j < numberOfStaves; j++ )
  {
    double phi = phi0 * rad;
    phi += j * sEndWSStepHolePhi[lay] * deg;
    phi -= sEndWSStepHoleTilt[lay];
    xpos = rpos * cos( phi );
    ypos = rpos * sin( phi );
    Ra.set( 0., 0., dphi - phi );
    Ta.set( xpos, ypos, -zpos );
    endWheel->AddPlacedVolume( stepSShLogVol, Ta, &Ra );
  }

  return;
}

//________________________________________________________________________________
void PHG4MvtxSupport::CreateConeLayers( G4AssemblyVolume *& av )
{
  for ( unsigned int iLay = 0; iLay < PHG4MvtxDefs::kNLayers; iLay++ )
  {
    GetConeVolume( iLay, av );
  }
}

//________________________________________________________________________________
void PHG4MvtxSupport::GetConeVolume( int lay, G4AssemblyVolume *& av )
{
  //
  // Creates the single CF cone
  // for a given layer of the MVTX detector
  // (Layer 0: MVTX-2-S-00005)
  // (Layer 1: MVTX-2-S-00032)
  // (Layer 2: MVTX-2-S-00047)
  //
  // Input:
  //         iLay : the layer number
  //
  // Output:
  //
  // Return:
  //
  // Created:      20 Jan 2023 Yasser Corrales Morales
  //

  const double sEndWheelSExtSectLen = 17 * mm;

  const double sThickness[3] = { 1. * mm, 1.12 * mm, 1.12 * mm };
  const double sBigCylDmax[3] = { 103 * mm, 149 * mm, 195 * mm };
  const double sBigCylDmin = sBigCylDmax[lay] - 2 * sThickness[lay];
  const double sSmallCylDmin[3] = { 59.94 * mm, 75.98 * mm, 91.48 * mm };
  const double sTotalLen[3] = { 186.8 * mm, 179.74 * mm, 223 * mm };
  const double sSmallCylLen[3] = { 91.86 * mm, 89.38 * mm, 85.38 * mm };
  const double sConePlusSCLen[3] = { 165.79 * mm, 158.51 * mm, 152.06 * mm };

  const int nZplanes = 4;
  const double zPlane[nZplanes] = { 0 * mm,
                                    -sSmallCylLen[lay],
                                    -sConePlusSCLen[lay],
                                    -sTotalLen[lay] };

 const double rInner[nZplanes] = { sSmallCylDmin[lay] / 2,
                                   sSmallCylDmin[lay] / 2,
                                   sBigCylDmin / 2,
                                   sBigCylDmin / 2 };

 const double rOuter[nZplanes] = { sSmallCylDmin[lay] / 2 + sThickness[lay],
                                   sSmallCylDmin[lay] / 2 + sThickness[lay],
                                   sBigCylDmax[lay] / 2,
                                   sBigCylDmax[lay] / 2 };

  auto coneSolid = new G4Polycone( Form( "cfConeL%d", lay ), 0., 2. * M_PI * rad,
                                   nZplanes, zPlane, rInner, rOuter );

  auto matCF = PHG4Detector::GetDetectorMaterial( "MVTX_CarbonFiber$" );

  auto coneLogVol = new G4LogicalVolume( coneSolid, matCF,
                                         Form( "ConeL%d_LOGIC", lay ),
                                         nullptr, nullptr, nullptr );

  m_DisplayAction->AddVolume( coneLogVol, "MVTX_CarbonFiber$" );

  double zpos = sEndWheelSNHolesZdist / 2 - ( sEndWStepHoleZpos + sEndWStepHoleZdist ) \
                + sEndWheelSExtSectLen;
  G4ThreeVector Ta = G4ThreeVector( 0., 0., -zpos );
  G4RotationMatrix Ra;
  av->AddPlacedVolume( coneLogVol, Ta, &Ra );

  return;
}

//________________________________________________________________________________
void PHG4MvtxSupport::CreateCYSS( G4AssemblyVolume *& av )
{
  //
  // Creates the CYSS
  // (MVTX-s-S-00080)
  // (MVTX-2-S-00075)
  // (MVTX-s-S-00079)
  // (MVTX-s-S-00077)
  // (MVTX-s-S-00078)
  // (MVTX-s-S-00076)
  //
  // Input:
  //         av : assembly volume
  //
  // Output:
  //
  // Return:
  //
  //
  // Created:     20 Jan 2023 Yasser Corrales Morales
  //

  // Let's start creating the North CYSS Flange
  // We split the CYSS flange in 2 polycone,
  // Although it is a single piece.
  const double sFlgExtThick = 1 * mm;
  const double sFlgIntThick = 1.5 * mm;
  const double sFlgSringLen = 4.5 * mm;
  const double sFlgLringLen = 7.5 * mm;
  const double sFlgRmin = 23.6 * mm;
  const double sFlgRmax = 52.8 * mm;
  const double sFlgSringRmin = 46.4 * mm;
  const double sFlgSringRmax = 47.4 * mm;
  const double sFlgLringRmin = 50.6 * mm;
  const double sFlgLringRmax = 51.5 * mm;

  const int nZplanes = 4;
  const double zPlane_1[nZplanes] = { 0. * mm,
                                      sFlgIntThick,
                                      sFlgIntThick,
                                      sFlgSringLen };

  const double rInner_1[nZplanes] = { sFlgRmin,
                                      sFlgRmin,
                                      sFlgSringRmin,
                                      sFlgSringRmin };

  const double rOuter_1[nZplanes] = { sFlgLringRmin,
                                      sFlgLringRmin,
                                      sFlgSringRmax,
                                      sFlgSringRmax };

  const double zPlane_2[nZplanes] = { 0. * mm,
                                      sFlgExtThick,
                                      sFlgExtThick,
                                      sFlgLringLen };

  const double rInner_2[nZplanes] = { sFlgLringRmin,
                                      sFlgLringRmin,
                                      sFlgLringRmin,
                                      sFlgLringRmin };

  const double rOuter_2[nZplanes] = { sFlgRmax,
                                      sFlgRmax,
                                      sFlgLringRmax,
                                      sFlgLringRmax };

  auto flangeSolid_1 = new G4Polycone( "cyssFlangeNorth_1", 0., 2. * M_PI * rad,
                                       nZplanes, zPlane_1, rInner_1, rOuter_1 );

  auto flangeSolid_2 = new G4Polycone( "cyssFlangeNorth_2", 0., 2. * M_PI * rad,
                                       nZplanes, zPlane_2, rInner_2, rOuter_2 );

  auto matAl = PHG4Detector::GetDetectorMaterial( "MVTX_EW_Al$" );

  auto cyssFlangeLogVol_1 = new G4LogicalVolume( flangeSolid_1, matAl,
                                                 "cyssFlangeNorth_1_LOGIC",
                                                 nullptr, nullptr, nullptr );

  m_DisplayAction->AddVolume( cyssFlangeLogVol_1, "MVTX_EW_Al$" );

  auto cyssFlangeLogVol_2 = new G4LogicalVolume( flangeSolid_2, matAl,
                                                 "cyssFlangeNorth_2_LOGIC",
                                                 nullptr, nullptr, nullptr );

  m_DisplayAction->AddVolume( cyssFlangeLogVol_2, "MVTX_EW_Al$");

  double zpos = sEndWheelSNHolesZdist / 2 - ( sEndWStepHoleZpos + sEndWStepHoleZdist ) \
                + sEndWheelNLen + sFlgIntThick;
  G4ThreeVector Ta = G4ThreeVector( 0., 0., zpos );
  G4RotationMatrix Ra( 0, M_PI * rad, 0. );
  av->AddPlacedVolume( cyssFlangeLogVol_1, Ta, &Ra );
  av->AddPlacedVolume( cyssFlangeLogVol_2, Ta, &Ra );

  // Now we create the CYSS cylinder
  const double sCYSScylRmax = 52.82 * mm;
  const double sCYSScylRmin = 51.70 * mm;
  const double sCYSScylLen = 425.77 * mm;

  auto cyssCylSol = new G4Tubs( "CYSScyl_SOLID", sCYSScylRmin, sCYSScylRmax, sCYSScylLen / 2,
                                0., 2 * M_PI );

  auto matCF = PHG4Detector::GetDetectorMaterial( "MVTX_CarbonFiber$" );
  auto cyssCylLog = new G4LogicalVolume( cyssCylSol, matCF, "CYSScyl_LOGIC",
                                         nullptr, nullptr, nullptr );
  m_DisplayAction->AddVolume( cyssCylLog, "MVTX_CarbonFiber$" );

  zpos -= ( sFlgIntThick - sCYSScylLen / 2  + sCYSSFlgSsfFlgNsf - sCYSSFlgSsfCylsf );
  Ta.set( 0., 0., zpos );
  Ra.set( 0, 0., 0. );
  av->AddPlacedVolume( cyssCylLog, Ta, &Ra );

  // and the CYSS Cone
  const double sCYSSconFlgDmin = 106.03 * mm;
  const double sCYSSconIntDmin = 211.00 * mm;
  const double sCYSSconFlgThick = 1.12 * mm;
  const double sCYSSconFlgLen = 13.48 * mm;
  const double sCYSSconNoseLen = 24.39 * mm;
  const double sCYSSconThick = 2.16 * mm;
  const double sCYSSconSlopeLen = 95.0 * mm;
  const double sCYSSconLen = 200.28 * mm;

  const int nZplanesCone = 5;

  const double zPlanesCone[nZplanesCone] = { 0 * mm,
                                             sCYSSconFlgLen,
                                             sCYSSconNoseLen,
                                             sCYSSconSlopeLen,
                                             sCYSSconLen };

  const double rInnerCone[nZplanesCone] = { sCYSSconFlgDmin / 2,
                                            sCYSSconFlgDmin / 2,
                                            sCYSSconFlgDmin / 2,
                                            sCYSSconIntDmin / 2,
                                            sCYSSconIntDmin / 2 };

  const double rOuterCone[nZplanesCone] = { ( sCYSSconFlgDmin + sCYSSconFlgThick ) / 2,
                                            ( sCYSSconFlgDmin + sCYSSconThick ) / 2,
                                            ( sCYSSconFlgDmin + sCYSSconThick ) / 2,
                                            ( sCYSSconIntDmin + sCYSSconThick ) / 2,
                                            ( sCYSSconIntDmin + sCYSSconThick ) / 2 };

  auto cyssConeSol = new G4Polycone( "cyssConeSol", 0., 2. * M_PI * rad,
                                     nZplanesCone, zPlanesCone, rInnerCone, rOuterCone );

  auto cyssConeLog = new G4LogicalVolume( cyssConeSol, matCF, "CYSScone_LOGIC",
                                          nullptr, nullptr, nullptr );
  m_DisplayAction->AddVolume( cyssConeLog, "MVTX_CarbonFiber$" );

  zpos -= ( sCYSScylLen / 2 + sCYSSFlgSsfCylsf - sCYSSconLen - sCYSSFlgSsfConesf );
  Ta.set( 0., 0., zpos );
  Ra.set( 0, M_PI * rad, 0. );
  av->AddPlacedVolume( cyssConeLog, Ta, &Ra );

  // Create Rib
  const double sCYSSribRint = 97.62 * mm;
  const double sCYSSribRext = 99.32 * mm;
  const double sCYSSribRmax = 105.30 * mm;
  const double sCYSSribWidth = 55.23 * mm;
  const double sCYSSribStep1pos = ( 15.6 + 7.7 ) / 2 * mm;
  const double sCYSSribStep2pos = ( 47.6 + 39.7 ) / 2  * mm;
  const double sCYSSribThick = 1.7 * mm;

  const int nZplanesRib = 10;
  const double zPlanesRib[nZplanesRib] = { 0. * mm,
                                           sCYSSribStep1pos,
                                           sCYSSribStep1pos,
                                           sCYSSribStep1pos + sCYSSribThick,
                                           sCYSSribStep1pos + sCYSSribThick,
                                           sCYSSribStep2pos,
                                           sCYSSribStep2pos,
                                           sCYSSribStep2pos + sCYSSribThick,
                                           sCYSSribStep2pos + sCYSSribThick,
                                           sCYSSribWidth };

  const double rInnerRib[nZplanesRib] = { sCYSSribRint,
                                          sCYSSribRint,
                                          sCYSSribRint,
                                          sCYSSribRint,
                                          sCYSSribRmax - sCYSSribThick,
                                          sCYSSribRmax - sCYSSribThick,
                                          sCYSSribRint,
                                          sCYSSribRint,
                                          sCYSSribRint,
                                          sCYSSribRint };

  const double rOuterRib[nZplanesRib] = { sCYSSribRext,
                                          sCYSSribRext,
                                          sCYSSribRmax,
                                          sCYSSribRmax,
                                          sCYSSribRmax,
                                          sCYSSribRmax,
                                          sCYSSribRmax,
                                          sCYSSribRmax,
                                          sCYSSribRext,
                                          sCYSSribRext };

  auto cyssRibSol = new G4Polycone( "cyssRibSol", 0., 2. * M_PI * rad,
                                     nZplanesRib, zPlanesRib, rInnerRib, rOuterRib );

  auto cyssRibLog = new G4LogicalVolume( cyssRibSol, matCF, "CYSSrib_LOGIC",
                                          nullptr, nullptr, nullptr );
  m_DisplayAction->AddVolume( cyssRibLog, "MVTX_CarbonFiber$" );

  zpos -= ( sCYSSconLen + sCYSSFlgSsfConesf - sCYSSFlgSsfRibsf );
  Ta.set( 0., 0., zpos );
  Ra.set( 0, 0., 0. );
  av->AddPlacedVolume( cyssRibLog, Ta, &Ra );

  // Flange South
  const double sFlgSRmax = 107.7 * mm;
  const double sFlgSRmin = 95.0 * mm;
  const double sFlgSRstep = 105.3 * mm;
  const double sFlgSChamfeEnd = 9. * mm;
  const double sFlgSTotalWidth = 20. * mm;
  const double sFlgSRimWidth = 7. * mm;

  const int nZplanesFlgS = 4;
  const double zPlanesFlgS[nZplanesFlgS] = { 0. * mm,
                                             sFlgSRimWidth,
                                             sFlgSRimWidth,
                                             sFlgSTotalWidth };

  const double rInnerFlgS[nZplanesFlgS] = { sFlgSRmin,
                                            sFlgSRmin,
                                            sFlgSRmin,
                                            sFlgSRstep - sFlgSChamfeEnd };

  const double rOuterFlgS[nZplanesFlgS] = { sFlgSRmax,
                                            sFlgSRmax,
                                            sFlgSRstep,
                                            sFlgSRstep };

  auto cyssFlgSSol = new G4Polycone( "cyssFlgSSol", 0., 2. * M_PI * rad,
                                     nZplanesFlgS, zPlanesFlgS, rInnerFlgS, rOuterFlgS );

  auto cyssFlgSLog = new G4LogicalVolume( cyssFlgSSol, matAl, "CYSSFlgS_LOGIC",
                                          nullptr, nullptr, nullptr );
  m_DisplayAction->AddVolume( cyssFlgSLog, "MVTX_EW_Al$" );

  zpos -= ( sCYSSFlgSsfRibsf );
  Ta.set( 0., 0., zpos );
  Ra.set( 0, 0., 0. );
  av->AddPlacedVolume( cyssFlgSLog, Ta, &Ra );

  return;
}

//________________________________________________________________________________
void PHG4MvtxSupport::CreateServiceBarrel( G4AssemblyVolume *& av )
{
  //
  // Creates the ServiceBarrel
  // (MVTX-2-S-00181)
  // (MVTX-2-S-00081)
  //
  // Input:
  //         av : assembly volume
  //
  // Output:
  //
  // Return:
  //
  //
  // Created:     21 Jan 2023 Yasser Corrales Morales
  //

  // Local Variables
  double zpos;

  G4ThreeVector Ta;
  G4RotationMatrix Ra;

  // MVTX North SB Flange
  const double sSBFlgNRmax = 107.7 * mm;
  const double sSBFlgNRmin = 95.0 * mm;
  const double sSBFlgNRimRmax = 105.3 * mm;
  const double sSBFlgNRimRmin = 103.3 * mm;
  const double sSBFlgNIntThick = 4.92 * mm;
  const double sSBFlgNExtThick = 7.00 * mm;
  const double sSBFlgNTotThick = 16.0 * mm;

  const int nZplanesFlgN = 5;
  const double zPlanesFlgN[nZplanesFlgN] = { 0. * mm,
                                            sSBFlgNIntThick,
                                            sSBFlgNExtThick,
                                            sSBFlgNExtThick,
                                            sSBFlgNTotThick };

  const double rInnerFlgN[nZplanesFlgN] = { sSBFlgNRmin,
                                            sSBFlgNRmin,
                                            sSBFlgNRimRmin,
                                            sSBFlgNRimRmin,
                                            sSBFlgNRimRmin };

  const double rOuterFlgN[nZplanesFlgN] = { sSBFlgNRmax,
                                            sSBFlgNRmax,
                                            sSBFlgNRmax,
                                            sSBFlgNRimRmax,
                                            sSBFlgNRimRmax };

  auto sbFlgNSol = new G4Polycone( "sbFlgNSol", 0., 2. * M_PI * rad,
                                     nZplanesFlgN, zPlanesFlgN, rInnerFlgN, rOuterFlgN );

  auto matAl = PHG4Detector::GetDetectorMaterial( "MVTX_EW_Al$" );

  auto sbFlgNLog = new G4LogicalVolume( sbFlgNSol, matAl, "sbFlgN_LOGIC",
                                          nullptr, nullptr, nullptr );
  m_DisplayAction->AddVolume( sbFlgNLog, "MVTX_CarbonFiber$" );

  zpos = sEndWheelSNHolesZdist / 2 - ( sEndWStepHoleZpos + sEndWStepHoleZdist ) \
         + sEndWheelNLen - sCYSSFlgSsfFlgNsf;
  Ta.set( 0., 0., zpos );
  Ra.set( 0, M_PI * rad, 0. );
  av->AddPlacedVolume( sbFlgNLog, Ta, &Ra );

  // SB Cylinder (To confirm)
  const double sSBcylRmin = 105.5 * mm;
  const double sSBcylRmax = 107.66 * mm;
  const double sSBcylLen = 1197.0 * mm;

  auto sbCylSol = new G4Tubs( "SBcyl_SOLID", sSBcylRmin, sSBcylRmax,
                              sSBcylLen / 2, 0., 2 * M_PI );

  auto matCF = PHG4Detector::GetDetectorMaterial( "MVTX_CarbonFiber$" );
  auto sbCylLog = new G4LogicalVolume( sbCylSol, matCF, "SBcyl_LOGIC",
                                         nullptr, nullptr, nullptr );
  m_DisplayAction->AddVolume( sbCylLog, "MVTX_CarbonFiber$" );

  zpos -= ( sSBcylLen / 2 + sSBFlgNExtThick );
  Ta.set( 0., 0., zpos );
  Ra.set( 0, 0., 0. );
  av->AddPlacedVolume( sbCylLog, Ta, &Ra );

  return;
}

//________________________________________________________________________________
void PHG4MvtxSupport::CreateCable( PHG4MvtxCable *object, G4AssemblyVolume &assemblyVolume )
{
  std::string cableMaterials[2] = {object->get_coreMaterial(), "G4_POLYETHYLENE"};

  float dX = object->get_xNorth() - object->get_xSouth();
  float dY = object->get_yNorth() - object->get_ySouth();
  float dZ = object->get_zNorth() - object->get_zSouth();

  float rotY = dZ != 0. ? std::atan(dX / dZ) : 0.;
  float rotZ = dX != 0. ? std::atan(dY / dX) : 0.;

  float setX = (object->get_xSouth() + object->get_xNorth()) / 2;
  float setY = (object->get_ySouth() + object->get_yNorth()) / 2;
  float setZ = (object->get_zSouth() + object->get_zNorth()) / 2;

  float length = std::sqrt(dX * dX + dY * dY + dZ * dZ);
  float IR[2] = {0, object->get_coreRadius()};
  float OR[2] = {object->get_coreRadius(), object->get_sheathRadius()};

  G4RotationMatrix rot;
  rot.rotateY(rotY);
  rot.rotateZ(rotZ);
  G4ThreeVector place;
  place.setX(setX);
  place.setY(setY);
  place.setZ(setZ);
  G4Transform3D transform(rot, place);
  // we need just one of these but have multiple calls to this method
  static G4UserLimits *g4userLimits = new G4UserLimits(0.01);

  for (int i = 0; i < 2; ++i)
  {
    G4Material *trackerMaterial = PHG4Detector::GetDetectorMaterial(cableMaterials[i]);

    G4VSolid *cylinderSolid = new G4Tubs(G4String(object->get_name() + "_SOLID"),
                                         IR[i], OR[i], (length / 2.), 0, 2 * M_PI);

    G4LogicalVolume *cylinderLogic = new G4LogicalVolume(cylinderSolid, trackerMaterial,
                                                         G4String(object->get_name() + "_LOGIC"), nullptr, nullptr, g4userLimits);

    if (i == 0)
    {
      m_DisplayAction->AddVolume(cylinderLogic, cableMaterials[i]);
    }
    else
    {
      m_DisplayAction->AddVolume(cylinderLogic, object->get_color());
    }

    assemblyVolume.AddPlacedVolume(cylinderLogic, transform);
  }
}

//________________________________________________________________________________
void PHG4MvtxSupport::CreateCableBundle(G4AssemblyVolume &assemblyVolume, const std::string &superName,
                                        bool enableSignal, bool enableCooling, bool enablePower,
                                        float x1, float x2, float y1, float y2, float z1, float z2)  //, float theta)
{
  //Set up basic MVTX cable bundle (24 Samtec cables, 1 power cable, 2 cooling cables)
  float samtecCoreRadius = 0.01275 * cm;
  float samtecSheathRadius = 0.05 * cm;
  float coolingStaveCoreRadius = 0.056 * cm;
  float coolingStaveSheathRadius = 0.1 * cm;
  float coolingCoreRadius = 0.125 * cm;
  float coolingSheathRadius = 0.2 * cm;  //?
  float powerLargeCoreRadius = 0.069 * cm;
  float powerLargeSheathRadius = 0.158 * cm;
  float powerMediumCoreRadius = 0.033 * cm;
  float powerMediumSheathRadius = 0.082 * cm;
  float powerSmallCoreRadius = 0.028 * cm;
  float powerSmallSheathRadius = 0.0573 * cm;  //?

  float globalShiftX = 0.;
  float globalShiftY = -0.0984 * cm;
  float samtecShiftX = -6 * samtecSheathRadius + globalShiftX;
  float samtecShiftY = 1 * samtecSheathRadius + globalShiftY;
  float coolingShiftX = -3 * coolingSheathRadius + globalShiftX;
  float coolingShiftY = -1 * coolingSheathRadius + globalShiftY;
  float powerShiftX = 3.5 * powerLargeSheathRadius + globalShiftX;
  float powerShiftY = 6.1 * powerLargeSheathRadius + globalShiftY;

  //Samtec cables (we use 24 as there are 12 twinax)
  if (enableSignal)
  {
    unsigned int nSamtecWires = 24;
    unsigned int nRow = 8;
    unsigned int nCol = nSamtecWires / nRow;
    for (unsigned int iRow = 0; iRow < nRow; ++iRow)
    {
      for (unsigned int iCol = 0; iCol < nCol; ++iCol)
      {
        float deltaX = samtecShiftX + ((iCol + 1) * (samtecSheathRadius * 2.6));
        float deltaY = samtecShiftY - ((iRow + 1) * (samtecSheathRadius * 2.1));
        PHG4MvtxCable *cable = new PHG4MvtxCable( boost::str(boost::format( "%s_samtec_%d_%d" ) \
                                                  % superName.c_str() % iRow % iCol),
                                                  "G4_Cu", samtecCoreRadius, samtecSheathRadius,
                                                  x1 + deltaX, x2 + deltaX, y1 + deltaY,
                                                  y2 + deltaY, z1, z2, "blue" );
        CreateCable(cable, assemblyVolume);
        delete cable;
      }
    }
  }

  //Cooling Cables
  if (enableCooling)
  {
    unsigned int nCool = 2;
    std::string cooling_color[2] = {"red", "white"};
    //std::regex_constants::icase - TO IGNORE CASE.
    auto rx = std::regex{ "MVTX_L([0-2])", std::regex_constants::icase };
    bool smallCooling = std::regex_search( superName, rx );
    PHG4MvtxCable *cable = nullptr;
    for ( unsigned int iCool = 0; iCool < nCool; ++iCool )
    {
      float coreRadius =  smallCooling ? coolingStaveCoreRadius : coolingCoreRadius;
      float sheathRadius = smallCooling ? coolingStaveSheathRadius : coolingSheathRadius;
      if ( ! smallCooling )
      {
        float deltaX = coolingShiftX + ( ( iCool + 1 ) * ( sheathRadius * 2 ) );
        float deltaY = coolingShiftY + ( sheathRadius * 2 );
        cable = new PHG4MvtxCable( boost::str( boost::format( "%s_cooling_%d" ) \
                                                              % superName.c_str() % iCool ),
                                                  "G4_WATER", coreRadius, sheathRadius,
                                                  x1 + deltaX, x2 + deltaX, y1 + deltaY,
                                                  y2 + deltaY, z1, z2, cooling_color[iCool] );
      } else {
        float deltaX = coolingShiftX + ( sheathRadius * 2 );
        float deltaY = coolingShiftY + ( ( iCool + 1 ) * ( sheathRadius * 2 ) );
        cable = new PHG4MvtxCable( boost::str( boost::format( "%s_cooling_%d" ) \
                                                              % superName.c_str() % iCool ),
                                                  "G4_WATER", coreRadius, sheathRadius,
                                                  x1 + deltaX, x2 + deltaX, y1 + deltaY,
                                                  y2 + deltaY, z1, z2, cooling_color[iCool] );
      }
      CreateCable( cable, assemblyVolume );
      delete cable;
    }
  }

  //Power Cables
  if (enablePower)
  {
    using PowerCableParameters = std::pair<std::pair<std::string, std::string>, std::pair<float, float>>;
    std::vector<PowerCableParameters> powerCables;
    std::vector<std::vector<float>> powerCableColors;

    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_digiReturn") % superName.c_str()), "Large"), std::make_pair((-2.5 * powerLargeSheathRadius) + powerShiftX, (-2.5 * powerLargeSheathRadius) + powerShiftY));
    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_digiSupply") % superName.c_str()), "Large"), std::make_pair((-4.5 * powerLargeSheathRadius) + powerShiftX, (-1.5 * powerLargeSheathRadius) + powerShiftY));
    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_anaReturn") % superName.c_str()), "Medium"), std::make_pair((-4 * powerLargeSheathRadius) + powerShiftX, (-5.75 * powerMediumSheathRadius) + powerShiftY));
    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_anaSupply") % superName.c_str()), "Medium"), std::make_pair((-5.1 * powerLargeSheathRadius) + powerShiftX, (-5.75 * powerMediumSheathRadius) + powerShiftY));
    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_digiSense") % superName.c_str()), "Small"), std::make_pair((-10 * powerSmallSheathRadius) + powerShiftX, (-1 * powerSmallSheathRadius) + powerShiftY));
    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_anaSense") % superName.c_str()), "Small"), std::make_pair((-8 * powerSmallSheathRadius) + powerShiftX, (-2 * powerSmallSheathRadius) + powerShiftY));
    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_bias") % superName.c_str()), "Small"), std::make_pair((-6 * powerSmallSheathRadius) + powerShiftX, (-3 * powerSmallSheathRadius) + powerShiftY));
    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_ground") % superName.c_str()), "Small"), std::make_pair((-4 * powerSmallSheathRadius) + powerShiftX, (-4 * powerSmallSheathRadius) + powerShiftY));

    for (PowerCableParameters &powerCable : powerCables)
    {
      float coreRad, sheathRad;
      std::string cableColor;
      std::string cableType = powerCable.first.second;
      std::string cableName = powerCable.first.first;
      if (cableType == "Small")
      {
        coreRad = powerSmallCoreRadius;
        sheathRad = powerSmallSheathRadius;
      }
      else if (cableType == "Medium")
      {
        coreRad = powerMediumCoreRadius;
        sheathRad = powerMediumSheathRadius;
      }
      else
      {
        coreRad = powerLargeCoreRadius;
        sheathRad = powerLargeSheathRadius;
      }

      if (cableName == boost::str(boost::format("%s_digiReturn") % superName.c_str())) cableColor = "black";
      if (cableName == boost::str(boost::format("%s_digiSupply") % superName.c_str())) cableColor = "red";
      if (cableName == boost::str(boost::format("%s_anaReturn") % superName.c_str())) cableColor = "black";
      if (cableName == boost::str(boost::format("%s_anaSupply") % superName.c_str())) cableColor = "red";
      if (cableName == boost::str(boost::format("%s_digiSense") % superName.c_str())) cableColor = "white";
      if (cableName == boost::str(boost::format("%s_anaSense") % superName.c_str())) cableColor = "green";
      if (cableName == boost::str(boost::format("%s_bias") % superName.c_str())) cableColor = "white";
      if (cableName == boost::str(boost::format("%s_ground") % superName.c_str())) cableColor = "green";

      PHG4MvtxCable *cable = new PHG4MvtxCable(powerCable.first.first, "G4_Cu", coreRad, sheathRad,
                                               (x1 + powerCable.second.first), (x2 + powerCable.second.first),
                                               (y1 + powerCable.second.second), (y2 + powerCable.second.second), z1, z2, cableColor);
      CreateCable(cable, assemblyVolume);
      delete cable;
    }
  }
}

//________________________________________________________________________________
G4AssemblyVolume *PHG4MvtxSupport::buildBarrelCable()
{
  G4AssemblyVolume *av = new G4AssemblyVolume();

  CreateCableBundle( *av, "barrelCable", true, true, false, 0, 0, 0, 0,
                     BarrelCableEnd, BarrelCableStart );
  CreateCableBundle( *av, "barrelCable", false, false, true, 0, 0, 0, 0,
                     BarrelCableEnd,
                     - sEndWheelSNHolesZdist / 2 + ( sEndWStepHoleZpos + sEndWStepHoleZdist ) - 40 * cm );
  return av;
}

//________________________________________________________________________________
G4AssemblyVolume *PHG4MvtxSupport::buildLayerCables( const int &lay )
{
  G4AssemblyVolume *av = new G4AssemblyVolume();
//  float rInner[3] = { 59.94 / 2 * mm, 75.98 / 2 * mm, 91.48 / 2 * mm };
  float rOuter[3] = { ( 103 - 2. ) / 2 * mm, ( 149 - 2.24 ) / 2 * mm, ( 195 - 2.24 ) / 2 * mm };
  float zConeLen[3] = { 186.8 * mm, 179.74 * mm, 223 * mm };
  float zMax = - sEndWheelSNHolesZdist / 2 + ( sEndWStepHoleZpos + sEndWStepHoleZdist )
               - 17 * mm - zConeLen[lay];
//  float zTransition2[3] = { -9.186 * cm, -8.938 *cm,  -8.538 * cm };
  CreateCableBundle( *av, Form( "MVTX_L%dCable", lay ), true, true, false,
                     rOuter[lay] - 3 * mm, rOuter[lay] - 5 * mm, 0, 0,
                     BarrelCableStart + 1 * mm, zMax );

//  float zCoolStart = - sEndWheelSNHolesZdist / 2 + ( sEndWStepHoleZpos + sEndWStepHoleZdist )
//                     - 17 * mm - zConeLen[lay];
//  CreateCableBundle( *av, Form( "MVTX_L%dCool_0", lay ), false, true, false, rInner[lay], rInner[lay], 0, 0,
//                     BarrelCableStart,  );
//  CreateCableBundle(*av, Form( "MVTX_L%dCable_1", lay ), true, false, false, rOuter[lay], rInner[lay], 0, 0, zTransition1[lay] + 0.1, zTransition2[lay]);
//  CreateCableBundle(*av, Form( "MVTX_L%dCable_2", lay ), true, false, false, rInner[lay], rInner[lay], 0, 0, zTransition2[lay] + 0.1, zMax);

  return av;
}

//________________________________________________________________________________
void PHG4MvtxSupport::ConstructMvtxSupport( G4LogicalVolume *&lv )
{
  CreateMvtxSupportMaterials();
  m_avSupport = new G4AssemblyVolume();

  CreateEndWheelsSideN( m_avSupport );
  CreateEndWheelsSideS( m_avSupport );
  CreateConeLayers( m_avSupport );
  CreateCYSS( m_avSupport );
  CreateServiceBarrel( m_avSupport );

  G4RotationMatrix Ra;
  G4ThreeVector Ta = G4ThreeVector();
  G4Transform3D Tr( Ra, Ta );
  m_avSupport->MakeImprint( lv, Tr, 0, m_overlapCheck );

  unsigned int nStaves[PHG4MvtxDefs::kNLayers];
  unsigned int totStaves = 0;
  for ( unsigned int i = 0; i < PHG4MvtxDefs::kNLayers; ++i )
  {
    nStaves[i] = (int) PHG4MvtxDefs::mvtxdat[i][PHG4MvtxDefs::kNStave];
    totStaves += nStaves[i];
  }

  m_avBarrelCable = buildBarrelCable();
  G4ThreeVector placeBarrelCable;
  for (unsigned int i = 0; i < totStaves; ++i)
  {
    float phi = (2.0 * M_PI / totStaves) * i;
    placeBarrelCable.setX((BarrelRadius - 1 * cm) * std::cos(phi));
    placeBarrelCable.setY((BarrelRadius - 1 * cm) * std::sin(phi));
    G4RotationMatrix rotBarrelCable;
    rotBarrelCable.rotateZ(phi + (-90. * deg));
    G4Transform3D transformBarrelCable(rotBarrelCable, placeBarrelCable);
    m_avBarrelCable->MakeImprint(lv, transformBarrelCable, 0, m_overlapCheck);
  }

  for (unsigned int iLayer = 0; iLayer < PHG4MvtxDefs::kNLayers; ++iLayer)
  {
    m_avLayerCable[iLayer] = buildLayerCables( iLayer );
    for (unsigned int iStave = 0; iStave < nStaves[iLayer]; ++iStave)
    {
      G4RotationMatrix rotCable;
      G4ThreeVector placeCable;
      float phi = ( 2.0 * M_PI / nStaves[iLayer] ) * iStave;
      placeCable.setX( std::cos(phi) );
      placeCable.setY( std::sin(phi) );
      rotCable.rotateZ(phi + ( ( 90. + cableRotate[iLayer] ) * deg ) );
      G4Transform3D transformCable(rotCable, placeCable);
      m_avLayerCable[iLayer]->MakeImprint(lv, transformCable, 0, m_overlapCheck);
    }
  }
}
