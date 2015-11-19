#include <PHG4TrackUserInfoV1.h>

#include <Geant4/G4Track.hh>
#include <boost/lexical_cast.hpp>

namespace PHG4TrackUserInfo
{
  void SetTrackIdOffset(G4Track* track, const int trkidoffset)
  {
    if ( G4VUserTrackInformation* p = track->GetUserInformation() )
      {
	// User info exists, test it for something valid
	if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
	  {
	    pp->SetTrackIdOffset(trkidoffset);
	  }
	else
	  {
	    std::cout << "Unknown UserTrackInformation stored in track number "
		      << boost::lexical_cast<std::string>(track->GetTrackID())
		      << std::endl;
	  }
      }
    else
      {
	// User info does not exist, add it.
	PHG4TrackUserInfoV1* pp = new PHG4TrackUserInfoV1();
	pp->SetTrackIdOffset(trkidoffset);
	track->SetUserInformation(pp);
      }
  }
  void SetUserTrackId(G4Track* track, const int usertrackid)
  {
    if ( G4VUserTrackInformation* p = track->GetUserInformation() )
      {
	// User info exists, test it for something valid
	if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
	  {
	    pp->SetUserTrackId(usertrackid);
	  }
	else
	  {
	    std::cout << "Unknown UserTrackInformation stored in track number "
		      << boost::lexical_cast<std::string>(track->GetTrackID())
		      << std::endl;
	  }
      }
    else
      {
	// User info does not exist, add it.
	PHG4TrackUserInfoV1* pp = new PHG4TrackUserInfoV1();
	pp->SetUserTrackId(usertrackid);
	track->SetUserInformation(pp);
      }
  }
  void SetWanted(G4Track* track, const int trkid)
  {
    if ( G4VUserTrackInformation* p = track->GetUserInformation() )
      {
	// User info exists, test it for something valid
	if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
	  {
	    pp->SetWanted(trkid);
	  }
	else
	  {
	    std::cout << "Unknown UserTrackInformation stored in track number "
		 << boost::lexical_cast<std::string>(track->GetTrackID())
		 << std::endl;
	  }
      }
    else
      {
	// User info does not exist, add it.
	PHG4TrackUserInfoV1* pp = new PHG4TrackUserInfoV1();
	pp->SetWanted(trkid);
	track->SetUserInformation(pp);
      }

  }
  void SetKeep(G4Track* track, const int trkid)
  {
    if ( G4VUserTrackInformation* p = track->GetUserInformation() )
      {
	// User info exists, test it for something valid
	if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
	  {
	    pp->SetKeep(trkid);
	  }
	else
	  {
	    std::cout << "Unknown UserTrackInformation stored in track number "
		      << boost::lexical_cast<std::string>(track->GetTrackID())
		      << std::endl;
	  }
    }
  else
    {
      // User info does not exist, add it.
      PHG4TrackUserInfoV1* pp = new PHG4TrackUserInfoV1();
      pp->SetKeep(trkid);
      track->SetUserInformation(pp);
    }

}
};
