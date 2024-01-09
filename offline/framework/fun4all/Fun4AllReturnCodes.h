// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLRETURNCODES_H
#define FUN4ALL_FUN4ALLRETURNCODES_H

// In general negative return codes signal some fatal condition where continuing is just a waste of cpu

// SubsysReco module return codes:
// ABORTRUN: signals that processing should be aborted for this run, the process will exit with a non zero exit
// ABORTEVENT: abort reconstruction of the current event, reset everything and process the next event
// EVENT_OK: generic good return
// DISCARDEVENT: tell a connected output manager not to save this event but continue processing
// DONOTREGISTERSUBSYSTEM: during registration, module indicates it doesn't want to be registered

// Synchronization
// SYNC_OK: all good
// SYNC_FAIL: synchronization failed (leads to a search for matching events by the sync manager)
// SYNC_NOOBJECT: no synchronization object in this input manager, take current event without check
// RESET_NODE_TREE: internal use - under certain conditions during the synchronization process Fun4All needs to reset the node tree to clear the prvious event


namespace Fun4AllReturnCodes
{
  enum
  {
    ABORTRUN = -2,
    ABORTEVENT = -1,
    EVENT_OK = 0,
    DISCARDEVENT = 1,
    SYNC_OK = 0,
    SYNC_FAIL = -1,
    SYNC_NOOBJECT = 1,
    DONOTREGISTERSUBSYSTEM = -3,
    RESET_NODE_TREE = 1
  };
}

#endif /* FUN4ALL_FUN4ALLRETURNCODES_H */
