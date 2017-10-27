#ifndef PHRANDOMSEED_H
#define PHRANDOMSEED_H

#ifndef __CINT__
#include <random>
#endif

//! standard way to get a random seed:
//! `unsigned int seed = PHRandomSeed();`
//! It return fix seed sequence if recoConsts RANDOMSEED is set.
//! otherwise it return a random seed from /dev/urandom
class PHRandomSeed
{
 public:
  PHRandomSeed() {}
  virtual ~PHRandomSeed() {}
  //! conversion operator for `unsigned int seed = PHRandomSeed();`
  operator unsigned int() const
  {
    return GetSeed();
  }

  //! get a seed
  static unsigned int GetSeed();

 protected:
  static void InitSeed();

#ifndef __CINT__
  static bool fInitialized;
  static std::mt19937 fRandomGenerator;
  static std::uniform_int_distribution<unsigned int> fDistribution;
#endif
};

#endif
