
#include "PHG4TrackInformation.h"

#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ios.hh>

G4Allocator<PHG4TrackInformation> aTrackInformationAllocator;

PHG4TrackInformation::PHG4TrackInformation()
{
  originalTrackID = 0;
  particleDefinition = 0;
  originalPosition = G4ThreeVector(0.,0.,0.);
  originalMomentum = G4ThreeVector(0.,0.,0.);
  originalEnergy = 0.;
  originalTime = 0.;
  trackingStatus = 1;
  sourceTrackID = -1;
  sourceTrackID = -1;
  sourceDefinition = 0;
  sourcePosition = G4ThreeVector(0.,0.,0.);
  sourceMomentum = G4ThreeVector(0.,0.,0.);
  sourceEnergy = 0.;
  sourceTime = 0.;
  _birthPlace = Undefined;
}

PHG4TrackInformation::PHG4TrackInformation(const G4Track* aTrack)
{
  originalTrackID = aTrack->GetTrackID();
  particleDefinition = aTrack->GetDefinition();
  originalPosition = aTrack->GetPosition();
  originalMomentum = aTrack->GetMomentum();
  originalEnergy = aTrack->GetTotalEnergy();
  originalTime = aTrack->GetGlobalTime();
  trackingStatus = 1;
  sourceTrackID = -1;
  sourceDefinition = 0;
  sourcePosition = G4ThreeVector(0.,0.,0.);
  sourceMomentum = G4ThreeVector(0.,0.,0.);
  sourceEnergy = 0.;
  sourceTime = 0.;
  _birthPlace = Undefined;
}

PHG4TrackInformation::PHG4TrackInformation(const PHG4TrackInformation* aTrackInfo)
{
  originalTrackID = aTrackInfo->originalTrackID;
  particleDefinition = aTrackInfo->particleDefinition;
  originalPosition = aTrackInfo->originalPosition;
  originalMomentum = aTrackInfo->originalMomentum;
  originalEnergy = aTrackInfo->originalEnergy;
  originalTime = aTrackInfo->originalTime;
  trackingStatus = aTrackInfo->trackingStatus;
  sourceTrackID = aTrackInfo->sourceTrackID;
  sourceDefinition = aTrackInfo->sourceDefinition;
  sourcePosition = aTrackInfo->sourcePosition;
  sourceMomentum = aTrackInfo->sourceMomentum;
  sourceEnergy = aTrackInfo->sourceEnergy;
  sourceTime = aTrackInfo->sourceTime;
  _birthPlace = aTrackInfo->_birthPlace;
}

PHG4TrackInformation::~PHG4TrackInformation()
{ ; }

PHG4TrackInformation& PHG4TrackInformation::operator =(const PHG4TrackInformation& aTrackInfo)
{
  originalTrackID = aTrackInfo.originalTrackID;
  particleDefinition = aTrackInfo.particleDefinition;
  originalPosition = aTrackInfo.originalPosition;
  originalMomentum = aTrackInfo.originalMomentum;
  originalEnergy = aTrackInfo.originalEnergy;
  originalTime = aTrackInfo.originalTime;
  trackingStatus = aTrackInfo.trackingStatus;
  sourceTrackID = aTrackInfo.sourceTrackID;
  sourceDefinition = aTrackInfo.sourceDefinition;
  sourcePosition = aTrackInfo.sourcePosition;
  sourceMomentum = aTrackInfo.sourceMomentum;
  sourceEnergy = aTrackInfo.sourceEnergy;
  sourceTime = aTrackInfo.sourceTime;
  _birthPlace = aTrackInfo._birthPlace;

  return *this;
}

void PHG4TrackInformation::SetSourceTrackInformation(const G4Track* aTrack)
{
  sourceTrackID = aTrack->GetTrackID();
  sourceDefinition = aTrack->GetDefinition();
  sourcePosition = aTrack->GetPosition();
  sourceMomentum = aTrack->GetMomentum();
  sourceEnergy = aTrack->GetTotalEnergy();
  sourceTime = aTrack->GetGlobalTime();
}

void PHG4TrackInformation::Print() const
{
  G4cout 
    << "Source track ID " << sourceTrackID << " (" << sourceDefinition->GetParticleName() << ","
    << sourceEnergy/GeV << "[GeV]) at " << sourcePosition << G4endl;
  G4cout
    << "Original primary track ID " << originalTrackID << " (" << particleDefinition->GetParticleName() << ","
    << originalEnergy/GeV << "[GeV])" << G4endl;
  G4cout << "Birthplace " << _birthPlace << G4endl;
}

