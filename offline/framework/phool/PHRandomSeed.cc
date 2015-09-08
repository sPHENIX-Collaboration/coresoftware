#include "PHRandomSeed.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

unsigned int PHRandomSeed()
{
  //  that's when we switch to c++11
  //  std::random_device rdev;
  //  uint32_t random_seed = rdev();
  unsigned int random_seed;
  // use /dev/urandom which unlike /dev/random is non blocking 
  // even if entropy is too low. Not for use in cryptographic apps but
  // good enough for a random seed
  FILE *fp = fopen("/dev/urandom", "r");
  if (!fp)
    {
      cout << "could not open /dev/urandom for seed" << endl;
      exit(1);
    }
  size_t readbytes = fread(&random_seed, sizeof(random_seed), 1,fp);
  fclose(fp);
  if (!readbytes)
    {
      cout << "PHRandomSeed: reading /dev/urandom failed" << endl;
      exit(1);
    }
  return random_seed;
}
