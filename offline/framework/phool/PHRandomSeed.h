#ifndef PHOOL_PHRANDOMSEED_H
#define PHOOL_PHRANDOMSEED_H

//! standard way to get a random seed:
//! `unsigned int seed = PHRandomSeed();`
//! It return fix seed sequence if recoConsts RANDOMSEED is set.
//! If values are preloaded via PHRandomSeed::LoadSeed, they are returned in loaded order
//! otherwise it return a random seed from std::random_device rdev
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
  static void LoadSeed(const unsigned int iseed);
  static void Verbosity(const int iverb);
  static int Verbosity() {return verbose;};

 protected:
  static void InitSeed();
  static bool fFixed;
  static bool fInitialized;
  static int verbose;
};

#endif
