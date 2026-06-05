#ifndef PHOOL_RUNNUMBERRANGE_H
#define PHOOL_RUNNUMBERRANGE_H

/**
 * Defines run-number range constants and special run markers used to identify physics data-taking periods.
 *
 * Each constant names the first or last run number (or a special marker) for a given data-taking period.
 *
 * @var RUN2PP_FIRST First Run 2 proton-proton physics run passing >=5m, >=100k evts.
 * @var RUN2PP_LAST  Last Run 2 proton-proton physics run.
 * @var RUN2AUAU_FIRST First Run 2 Au+Au (heavy-ion) physics run.
 * @var RUN2AUAU_LAST  Last Run 2 Au+Au (heavy-ion) physics run.
 * @var RUN3_TPCFW_CLOCK_CHANGE Run 3 marker for the TPC Forward clock change.
 * @var RUN3AUAU_FIRST First Run 3 Au+Au (heavy-ion) physics run.
 * @var RUN3AUAU_LAST  Last Run 3 Au+Au (heavy-ion) physics run.
 * @var RUN3PP_FIRST First Run 3 proton-proton (beam) physics run.
 * @var RUN3PP_LAST  Last Run 3 proton-proton physics run.
 * @var RUN3OO_FIRST First Run 3 O+O physics run.
 * @var RUN3OO_LAST  Last Run 3 O+O physics run.
 */
namespace RunnumberRange
{
  constexpr int RUN2PP_FIRST = 47287;
  constexpr int RUN2PP_LAST = 53880;
  constexpr int RUN2AUAU_FIRST = 54128;
  constexpr int RUN2AUAU_LAST = 54974;
  constexpr int RUN3_TPCFW_CLOCK_CHANGE = 58667;
  constexpr int RUN3AUAU_FIRST = 66457;
  constexpr int RUN3AUAU_LAST = 78954;
  constexpr int RUN3PP_FIRST = 79146; // first beam data
  constexpr int RUN3PP_LAST = 81668;
  constexpr int RUN3OO_FIRST = 82388; // after trigger settled down (run 82374 excluded);
  constexpr int RUN3OO_LAST = 82703;
}

#endif
