#ifndef PHG4UserPrimaryParticleInformation_H__
#define PHG4UserPrimaryParticleInformation_H__

#include <Geant4/G4VUserPrimaryParticleInformation.hh>
#include <iostream>

class PHG4UserPrimaryParticleInformation : public G4VUserPrimaryParticleInformation
{
public:
  PHG4UserPrimaryParticleInformation(const int emb) : embed(emb) {}
  void Print() const 
  {
    std::cout << "Embedding = " << embed << std::endl;
  }
  int get_embed() const {return embed;}

private:
  int embed;
};

#endif
