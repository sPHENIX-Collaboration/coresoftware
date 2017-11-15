#include "PHRandomSeed.h"
#include "recoConsts.h"

#include <cstdlib>
#include <iostream>
#include <queue>
#include <random>

using namespace std;

static queue<unsigned int> seedqueue;
static std::mt19937 fRandomGenerator;
static std::uniform_int_distribution<unsigned int> fDistribution;

bool PHRandomSeed::fInitialized(false);
bool PHRandomSeed::fFixed(false);

unsigned int PHRandomSeed::GetSeed()
{
  unsigned int iseed;
  if (!seedqueue.empty())
  {
    iseed = seedqueue.front();
    seedqueue.pop();
  }
  else
  {
    if (!fInitialized)
    {
      InitSeed();
    }
    if (fFixed)
    {
      iseed = fDistribution(fRandomGenerator);
    }
    else
    {
      std::random_device rdev;
      iseed = rdev();
    }
  }
  cout << "PHRandomSeed::GetSeed() seed: " << iseed << endl;
  return iseed;
}

void PHRandomSeed::InitSeed()
{
  recoConsts *rc = recoConsts::instance();
  if (rc->FlagExist("RANDOMSEED"))
  {
    // fixed init seed
    const unsigned int seed = rc->get_IntFlag("RANDOMSEED");
    cout << "PHRandomSeed: using fixed seed " << seed << endl;
    fRandomGenerator.seed(seed);
    fFixed = true;
    fInitialized = true;
  }
}

void PHRandomSeed::LoadSeed(unsigned int iseed)
{
  seedqueue.push(iseed);
}
