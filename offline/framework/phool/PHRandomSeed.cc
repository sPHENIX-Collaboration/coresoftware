#include "PHRandomSeed.h"
#include "recoConsts.h"

#include <iostream>
#include <memory>  // for allocator
#include <queue>
#include <random>

static std::queue<unsigned int> seedqueue;
static std::mt19937 fRandomGenerator;
static std::uniform_int_distribution<unsigned int> fDistribution;

bool PHRandomSeed::fInitialized(false);
bool PHRandomSeed::fFixed(false);
int PHRandomSeed::verbose(1);

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
  if (verbose)
  {
    std::cout << "PHRandomSeed::GetSeed() seed: " << iseed << std::endl;
  }
  return iseed;
}

void PHRandomSeed::InitSeed()
{
  recoConsts *rc = recoConsts::instance();
  if (rc->FlagExist("RANDOMSEED"))
  {
    // fixed init seed
    const unsigned int seed = rc->get_IntFlag("RANDOMSEED");
    std::cout << "PHRandomSeed: using fixed seed " << seed << std::endl;
    fRandomGenerator.seed(seed);
    fFixed = true;
    fInitialized = true;
  }
}

void PHRandomSeed::LoadSeed(const unsigned int iseed)
{
  seedqueue.push(iseed);
}

void PHRandomSeed::Verbosity(const int iverb)
{
  verbose = iverb;
}
