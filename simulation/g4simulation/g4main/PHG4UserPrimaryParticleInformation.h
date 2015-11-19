#ifndef PHG4UserPrimaryParticleInformation_H__
#define PHG4UserPrimaryParticleInformation_H__

#include <Geant4/G4VUserPrimaryParticleInformation.hh>
#include <iostream>

class PHG4UserPrimaryParticleInformation : public G4VUserPrimaryParticleInformation
{
public:
  PHG4UserPrimaryParticleInformation(const int emb) : embed(emb),
						      usertrackid(0) {}
  void Print() const 
  {
    std::cout << "Embedding = " << embed << std::endl;
    std::cout << "User Track ID = " << usertrackid << std::endl;
  }
  int get_embed() const {return embed;}

  void set_user_track_id(int val) {usertrackid = val;}
  int get_user_track_id() const {return usertrackid;}
  
private:
  int embed;
  int usertrackid;
};

#endif
