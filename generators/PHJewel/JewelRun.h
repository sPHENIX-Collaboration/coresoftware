// JewelRun.h -- the one way PHJewel starts JEWEL:
//
//   1. JewelCard::validate   strict reading of the user's card pair
//   2. writeCanonical        <stem>.card.dat / <stem>.medium.dat next to the log
//   3. JewelCore::load/init  core library for the medium, canonical card, seed
//   4. cross-check           what JEWEL read back against what we validated
//
// Any failure throws std::runtime_error carrying every message, so a bad card
// reports all of its problems at once.
#ifndef PHJEWEL_JEWELRUN_H
#define PHJEWEL_JEWELRUN_H

#include "JewelCard.h"

#include <iostream>
#include <string>

class JewelCore;

struct JewelRunConfig
{
  std::string card;          // user's main card
  std::string medium;        // simple | brick | vac
  std::string stem;          // output stem: <stem>.log, <stem>.card.dat, ...
  int seed = 0;              // PYTHIA MRPY(1)
  std::string coreLib;       // "" -> libjewel_<medium>.so
  std::string fortranHepMC;  // "" or a path (validation)
  std::string evinfo;        // "" or a path
  bool allowTableGeneration = false;
  long isospinBlock = 1000;
};

// Returns the validation result (canonical cards, hash, notes).
JewelCard::Result startJewel(JewelCore &core, const JewelRunConfig &c, std::ostream &os = std::cout);

#endif
