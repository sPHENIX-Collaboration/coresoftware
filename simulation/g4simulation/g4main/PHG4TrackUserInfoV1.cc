#include "PHG4TrackUserInfoV1.h"

#include <Geant4/G4Track.hh>
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation

#include <boost/lexical_cast.hpp>

#include <iostream>  // for operator<<, basic_ostream, endl, cout
#include <string>    // for string, operator<<

namespace PHG4TrackUserInfo
{
  void SetUserTrackId(G4Track* track, const int usertrackid)
  {
    if (G4VUserTrackInformation* p = track->GetUserInformation())
    {
      // User info exists, test it for something valid
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
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

  void SetUserParentId(G4Track* track, const int userparentid)
  {
    if (G4VUserTrackInformation* p = track->GetUserInformation())
    {
      // User info exists, test it for something valid
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
      {
        pp->SetUserParentId(userparentid);
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
      pp->SetUserParentId(userparentid);
      track->SetUserInformation(pp);
    }
  }

  void SetUserPrimaryId(G4Track* track, const int userprimaryid)
  {
    if (G4VUserTrackInformation* p = track->GetUserInformation())
    {
      // User info exists, test it for something valid
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
      {
        pp->SetUserPrimaryId(userprimaryid);
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
      pp->SetUserPrimaryId(userprimaryid);
      track->SetUserInformation(pp);
    }
  }

  void SetWanted(G4Track* track, const int trkid)
  {
    if (G4VUserTrackInformation* p = track->GetUserInformation())
    {
      // User info exists, test it for something valid
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
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
    if (G4VUserTrackInformation* p = track->GetUserInformation())
    {
      // User info exists, test it for something valid
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
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

  void SetShower(G4Track* track, PHG4Shower* shower)
  {
    if (G4VUserTrackInformation* p = track->GetUserInformation())
    {
      // User info exists, test it for something valid
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
      {
        pp->SetShower(shower);
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
      pp->SetShower(shower);
      track->SetUserInformation(pp);
    }
  }

}  // namespace PHG4TrackUserInfo
