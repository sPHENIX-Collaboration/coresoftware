#ifndef __PHG4TrackUserInfoV1_H__
#define __PHG4TrackUserInfoV1_H__

#include <Geant4/G4VUserTrackInformation.hh>

// Made this with "V1" in the name in case we ever want to inherit from
// it with other versions...

// Use the UserTrackInformation to attach a flag telling the framework 
// to save the track in the truth output.  Other uses might include keeping
// track of the 

class PHG4TrackUserInfoV1 : public G4VUserTrackInformation
{
public:
  PHG4TrackUserInfoV1() : G4VUserTrackInformation("TrackUserInfoV1"),
			  trackidoffset(0), usertrackid(0), wanted(0), keep(0) {}
  virtual ~PHG4TrackUserInfoV1() {}
  void Print() const 
  {
    G4cout << "PHG4TrackUserInfoV1: " << std::endl;
    G4cout << "   TrackIdOffset = " << trackidoffset << std::endl;
    G4cout << "   UserTrackId = " << usertrackid << std::endl;
    G4cout << "   Wanted = " << wanted << std::endl;
    G4cout << "   Keep = " << keep << std::endl;
  }
  void SetTrackIdOffset (const int val) { trackidoffset = val; }
  //int GetTrackIdOffset() const {return trackidoffset;}
  void SetUserTrackId(const int val) {usertrackid = val;}
  int GetUserTrackId() const {return usertrackid;}
  void SetWanted(const int val) {wanted = val;}
  int GetWanted() const {return wanted;}
  void SetKeep(const int val) {keep = val;}
  int GetKeep() const {return keep;}

private:
  int trackidoffset;
  int usertrackid;
  int wanted;
  int keep;
};

// Utility function to wrap up the operations involved with adding user info
// to the track.
class G4Track;

namespace PHG4TrackUserInfo 
{
  void SetTrackIdOffset(G4Track* track, const int trkidoffset);
  void SetUserTrackId(G4Track* track, const int usertrackid);
  void SetWanted(G4Track* track, const int wanted);
  void SetKeep(G4Track* track, const int keep);
};

#endif
