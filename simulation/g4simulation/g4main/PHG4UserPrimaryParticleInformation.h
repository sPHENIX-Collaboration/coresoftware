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
    std::cout << "User Vertex ID = " << uservtxid << std::endl;
  }
  int get_embed() const {return embed;}

  void set_user_track_id(int val) {usertrackid = val;}
  int get_user_track_id() const {return usertrackid;}

  void set_user_vtx_id(int val) {uservtxid = val;}
  int get_user_vtx_id() const {return uservtxid;}
  
private:
  int embed;
  int usertrackid;
  int uservtxid;
};

#endif
