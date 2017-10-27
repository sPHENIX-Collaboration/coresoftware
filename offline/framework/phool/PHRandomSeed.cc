#include "PHRandomSeed.h"
#include "recoConsts.h"

#include <gsl/gsl_const.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <random>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

bool PHRandomSeed::fInitialized(false);
mt19937 PHRandomSeed::fRandomGenerator;
uniform_int_distribution<unsigned int> PHRandomSeed::fDistribution;

unsigned int PHRandomSeed::GetSeed()
{
  if (!fInitialized)
    InitSeed();

  assert(fInitialized);
  return fDistribution(fRandomGenerator);
}

void PHRandomSeed::InitSeed()
{
  recoConsts *rc = recoConsts::instance();
  if (rc->FlagExist("RANDOMSEED"))
  {
    // fixed init seed
    const int seed = rc->get_IntFlag("RANDOMSEED");
    fRandomGenerator.seed(seed);
    fInitialized = true;

    cout << "PHRandomSeed::InitSeed - initialized with fixed random seed " << seed << " from recoConsts.RANDOMSEED. "<<endl
        <<  "                         To reproduce this random number sequences, please add this line to Fun4All macro: "<<endl
        <<  "                         recoConsts::instance()->set_IntFlag(\"RANDOMSEED\", " << seed << ");";
  }
  else
  {
    // random init seed

    //  that's when we switch to c++11
    //  std::random_device rdev;
    //  uint32_t random_seed = rdev();
    unsigned int random_seed(0);
    // use /dev/urandom which unlike /dev/random is non blocking
    // even if entropy is too low. Not for use in cryptographic apps but
    // good enough for a random seed
    FILE *fp = fopen("/dev/urandom", "r");
    if (!fp)
    {
      cout << "could not open /dev/urandom for seed" << endl;
      exit(1);
    }
    size_t readbytes = fread(&random_seed, sizeof(random_seed), 1, fp);
    fclose(fp);
    if (!readbytes)
    {
      cout << "PHRandomSeed: reading /dev/urandom failed" << endl;
      exit(1);
    }

    fRandomGenerator.seed(random_seed);
    fInitialized = true;

    cout << "PHRandomSeed::InitSeed - initialized with random seed " << random_seed << " from /dev/urandom. "<<endl
         << "                         To reproduce this random number sequences, please add this line to Fun4All macro"<<endl
         << "                         recoConsts::instance()->set_IntFlag(\"RANDOMSEED\", " << random_seed << ");";
  }
}
