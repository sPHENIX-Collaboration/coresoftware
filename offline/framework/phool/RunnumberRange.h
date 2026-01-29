#ifndef PHOOL_RUNNUMBERRANGE_H
#define PHOOL_RUNNUMBERRANGE_H

/**
 * Defines run-number range constants and special run markers used to identify physics data-taking periods.
 *
 * Each constant names the first or last run number (or a special marker) for a given data-taking period.
 *
 * @var RUN2PP_FIRST First Run 2 proton-proton physics run.
 * @var RUN2PP_LAST  Last Run 2 proton-proton physics run.
 * @var RUN2AUAU_FIRST First Run 2 Au+Au (heavy-ion) physics run.
 * @var RUN2AUAU_LAST  Last Run 2 Au+Au (heavy-ion) physics run.
 * @var RUN3_TPCFW_CLOCK_CHANGE Run 3 marker for the TPC Forward clock change.
 * @var RUN3AUAU_FIRST First Run 3 Au+Au (heavy-ion) physics run.
 * @var RUN3AUAU_LAST  Last Run 3 Au+Au (heavy-ion) physics run.
 * @var RUN3PP_FIRST First Run 3 proton-proton (beam) physics run.
 * @var RUN3PP_LAST  Last Run 3 proton-proton physics run.
 * @var RUN3OO_FIRST Temporary placeholder for the first Run 3 OO run (to be updated once OO starts).
 * @var RUN3OO_LAST  Temporary upper bound for Run 3 OO runs.
 */
namespace RunnumberRange
{
  static const int RUN2PP_FIRST = 47286;
  static const int RUN2PP_LAST = 53880;
  static const int RUN2AUAU_FIRST = 54128;
  static const int RUN2AUAU_LAST = 54974;
  static const int RUN3_TPCFW_CLOCK_CHANGE = 58667;
  static const int RUN3AUAU_FIRST = 66457;
  static const int RUN3AUAU_LAST = 78954;
  static const int RUN3PP_FIRST = 79146; // first beam data
  static const int RUN3PP_LAST = 81668;
  static const int RUN3OO_FIRST = 82300; // TEMP (to be updated once OO starts)
  static const int RUN3OO_LAST = 200000;
}

#endif