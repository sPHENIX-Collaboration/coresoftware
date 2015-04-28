// Implementation of an eta parameterization in a single phi slice for a Tubs.
// Shamelessly stolen from the G4ParameterisationTubs implementation (based on
// the Z subclass.  See $G4INSTALL/source/volumes/divisions)
//

#include "PHG4ParameterisationTubsEta.h"

#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4ReflectedSolid.hh>
#include <Geant4/G4Tubs.hh>

#include <iomanip>

//--------------------------------------------------------------------------
PHG4ParameterisationTubsEta::PHG4ParameterisationTubsEta( EAxis axis, G4int nDiv,
							  G4double width, G4double offset,
							  G4VSolid* msolid, DivisionType divType ) :
  G4VParameterisationTubs( axis, nDiv, width, offset, msolid, divType )
{ 
  CheckParametersValidity();
  SetType( "DivisionTubsZ" );

  //  G4Tubs* msol = (G4Tubs*)(fmotherSolid);
//   if( divType == DivWIDTH )
//   {
//     fnDiv = CalculateNDiv( 2*msol->GetZHalfLength(), width, offset );
//   }
//   else if( divType == DivNDIV )
//   {
//     fwidth = CalculateWidth( 2*msol->GetZHalfLength(), nDiv, offset );
//   }

  G4Tubs* msol = (G4Tubs*)(fmotherSolid);
  double radius = msol->GetInnerRadius();

  // TODO: allow for a placement offset in Z
  double _centerZ = 0.0;

  double minZ = -msol->GetZHalfLength() ;
  double maxZ = -minZ;
  double minR = sqrt(radius*radius + minZ*minZ);
  double minEta = -std::log((minR-minZ)/radius);
  double maxR = sqrt(radius*radius + maxZ*maxZ);
  double maxEta = -std::log((maxR-maxZ)/radius);

  double totalEta = maxEta - minEta;
  double dEta = totalEta / GetNoDiv();
  double zmin = _centerZ + radius * std::sinh(minEta);
  for(int i=0; i<GetNoDiv(); i++)
    {
      // Compute the edges of this eta bin
      double etaMin = minEta + dEta * i;
      double etaMax = etaMin + dEta;

      // Compute the corresponding Z positions of the edges
      //double zmin = _centerZ + radius * std::sinh(etaMin);
      double zmax = _centerZ + radius * std::sinh(etaMax);

      // Z positions is halfway between the edges
      double zpos = (zmin+zmax)/2.0;
      double zhalf = (zmax-zmin)/2.0;
      _zpos.push_back(zpos);
      _zhalf.push_back(zhalf);

//       std::cout << i << ": " << etaMin << " " << etaMax << " " 
// 		<< zmin << " " << zmax << " " << zpos << " +/- " << zhalf << std::endl;

      // Start the next iteration with zmax as the new value for zmin
      zmin = zmax;
    }

//   G4cout << "G4ParameterisationTubsEta: Compare internal to passed args:" << G4endl
// 	 << " # divisions " << fnDiv << " = " << nDiv << G4endl
// 	 << " Offset " << foffset << " = " << offset << G4endl
// 	 << " Width " << fwidth << " = " << width << G4endl;
//   std::cout << "Min/Max Z = " << _zpos.front()-_zhalf.front() << " / " << _zpos.back()+_zhalf.back() << std::endl;
//   std::cout << "Min/Max Eta = " << minEta << " / " << maxEta << std::endl;
}

//--------------------------------------------------------------------------
PHG4ParameterisationTubsEta::~PHG4ParameterisationTubsEta()
{
}

//------------------------------------------------------------------------
G4double PHG4ParameterisationTubsEta::GetMaxParameter() const
{
  G4Tubs* msol = (G4Tubs*)(fmotherSolid);
  return 2*msol->GetZHalfLength();
}

//--------------------------------------------------------------------------
void
PHG4ParameterisationTubsEta::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const
{
  //----- set translation: along Z axis
  //G4Tubs* motherTubs = (G4Tubs*)(fmotherSolid);
//   G4double posi = - motherTubs->GetZHalfLength() + OffsetZ() 
//                   + fwidth/2 + copyNo*fwidth;
  double posi = _zpos.at(copyNo);
  G4ThreeVector origin(0.,0.,posi); 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix: unit

#if 0
  {
    G4cout << " PHG4ParameterisationTubsEta::ComputeTransformation()" << G4endl
           << " Position: " << posi << " - copyNo: " << copyNo << G4endl
           << " foffset " << foffset/deg << " - fwidth " << fwidth/deg
           << G4endl;
  }
#endif 

  ChangeRotMatrix( physVol );

#if 0
  {
    G4cout << std::setprecision(8) << " PHG4ParameterisationTubsEta " << copyNo
           << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl; 
  }
#endif
}

//--------------------------------------------------------------------------
void
PHG4ParameterisationTubsEta::ComputeDimensions( G4Tubs& tubs, const G4int copyNo,
						const G4VPhysicalVolume* ) const
{
  G4Tubs* msol = (G4Tubs*)(fmotherSolid);

  G4double pRMin = msol->GetInnerRadius();
  G4double pRMax = msol->GetOuterRadius();
  //  G4double pDz = msol->GetZHalfLength() / GetNoDiv();
  //G4double pDz = fwidth/2.;
  G4double pDz = _zhalf.at(copyNo);
  G4double pSPhi = msol->GetStartPhiAngle();
  G4double pDPhi = msol->GetDeltaPhiAngle();

  tubs.SetInnerRadius( pRMin );
  tubs.SetOuterRadius( pRMax );
  tubs.SetZHalfLength( pDz );
  tubs.SetStartPhiAngle( pSPhi );
  tubs.SetDeltaPhiAngle( pDPhi );

#if 0
  {
    G4cout << " PHG4ParameterisationTubsEta::ComputeDimensions()" << G4endl
           << " pDz: " << pDz << G4endl;
    tubs.DumpInfo();
  }
#endif
}
