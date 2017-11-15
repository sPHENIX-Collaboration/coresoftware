#include "PHRandomSeed.h"
#include "recoConsts.h"

#include <random>
#include <cstdlib>
#include <iostream>
#include <queue>

using namespace std;

static queue<unsigned int> seedqueue;
static std::mt19937 fRandomGenerator;
static std::uniform_int_distribution<unsigned int> fDistribution;

bool PHRandomSeed::fInitialized(false);
bool PHRandomSeed::fFixed(false);

unsigned int PHRandomSeed::GetSeed()
{
  if (! seedqueue.empty())
  {
    unsigned int iseed = seedqueue.front();
    seedqueue.pop();
    return iseed;
  }
  if (!fInitialized)
  {
    InitSeed();
  }
  if (fFixed)
  {
  return fDistribution(fRandomGenerator);
  }
     std::random_device rdev;
      uint32_t random_seed = rdev();
      return random_seed;
}

void PHRandomSeed::InitSeed()
{
  recoConsts *rc = recoConsts::instance();
  if (rc->FlagExist("RANDOMSEED"))
  {
    // fixed init seed
    const int seed = rc->get_IntFlag("RANDOMSEED");
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
