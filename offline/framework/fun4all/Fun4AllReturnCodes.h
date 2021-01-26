// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLRETURNCODES_H
#define FUN4ALL_FUN4ALLRETURNCODES_H

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

#endif /* __FUN4ALLRETURNCODES_H__ */
