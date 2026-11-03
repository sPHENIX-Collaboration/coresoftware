// JewelCard.h -- strict reading of JEWEL 2.6.0 parameter cards for PHJewel.
//
// JEWEL's own reader (init() and MEDINIT()) never fails: an unknown key, a
// value it cannot parse, a missing medium file or a forgotten key all leave
// the built-in default in place, and the defaults are LHC Pb+Pb (SQRTS 2760,
// MASS 208, TI 0.59, SIGMANN 6.2).  A RHIC card with one typo therefore runs
// to completion as a different physics sample.  This reader is the guard:
//
//   * every key must be known, appear once, and parse completely;
//   * the keys that define the sample have no default and must be present;
//   * keys PHJewel owns (seed, event count, output files) are noted, ignored;
//   * values that are wrong inside a detector simulation are refused
//     (WRITEDUMMIES T, DOSUBTRACTION T, HADRO F, PROCESS EEJJ, ...);
//   * main card and medium card are cross-checked (NF, MASS, recoils);
//   * the grammar accepted is the subset that standalone JEWEL reads the same
//     way, so one card file drives both and comparisons stay honest.
//
// The result is a canonical pair of cards -- every parameter explicit,
// defaults written out, absolute paths -- which is what JEWEL is actually
// given, kept next to the log, and hashed into RUN/PHGenIntegral.
//
// Only PHJewel uses it; the header is not installed.
#ifndef PHJEWEL_JEWELCARD_H
#define PHJEWEL_JEWELCARD_H

#include <map>
#include <set>
#include <string>
#include <vector>

namespace JewelCard
{
  // ---- the full, cross-checked configuration --------------------------------
  struct Options
  {
    std::string medium;                 // simple | brick | vac
    bool allowTableGeneration = false;  // let JEWEL integrate missing tables
  };

  struct Result
  {
    bool ok = false;
    std::vector<std::string> errors;
    std::vector<std::string> notes;

    std::string mainCanon;    // canonical main card text
    std::string mediumCanon;  // canonical medium card text ("" for vac)
    std::string hash;         // 16 hex digits over both canonical cards

    // the values PHJewel needs, and checks JEWEL's own reading against
    std::map<std::string, std::string> main;             // key -> canonical value
    std::map<std::string, std::string> medium;           // key -> canonical value
    std::set<std::string> explicitMain, explicitMedium;  // keys the user wrote
    std::string mediumPath;                              // absolute; "" for vac
  };

  // check the card pair; Result::ok false means r.errors says why
  Result validate(const std::string &cardPath, const Options &opt);

  // Write the canonical pair next to the log: <stem>.card.dat (+ .medium.dat).
  // Returns the path of the main card JEWEL is to read, "" (and err) on failure.
  std::string writeCanonical(const Result &r, const std::string &stem, const std::string &source,
                             const std::string &medium, std::string &err);

  // absolute path of p relative to base ("." = the working directory), with "." and ".." collapsed
  std::string absPath(const std::string &p, const std::string &base);

  // canonical values back as numbers
  double R(const std::map<std::string, std::string> &m, const char *k);
  long I(const std::map<std::string, std::string> &m, const char *k);
  bool L(const std::map<std::string, std::string> &m, const char *k);
}  // namespace JewelCard

#endif
