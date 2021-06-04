// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4USERPRIMARYPARTICLEINFORMATION_H
#define G4MAIN_PHG4USERPRIMARYPARTICLEINFORMATION_H

#include <Geant4/G4VUserPrimaryParticleInformation.hh>
#include <iostream>

class PHG4UserPrimaryParticleInformation : public G4VUserPrimaryParticleInformation
{
public:
  PHG4UserPrimaryParticleInformation(const int emb) : 
    embed(emb),
    usertrackid(0),
    uservtxid(0),
    barcode(-1)   {}

  void Print() const override 
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

  void set_user_barcode(int bcd) {barcode = bcd;}
  int get_user_barcode() const {return barcode;}
  
private:
  int embed;
  int usertrackid;
  int uservtxid;
  int barcode; 
};

#endif
