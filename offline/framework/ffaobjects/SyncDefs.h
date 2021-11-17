// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_SYNCDEFS_H
#define FFAOBJECTS_SYNCDEFS_H

#include <string>

namespace syncdefs
{
  // __attribute__((unused)) prevents compiler from flagging variable as unused
  static const int NUM_SYNC_VARS __attribute__((unused)) = 4;
  static const char *SYNCVARS[] __attribute__((unused)) = {"eventcounter", "eventnumber", "runnumber", "segmentnumber"};
  static const std::string SYNCNODENAME __attribute__((unused)) = "Sync";
}  // namespace syncdefs

#endif
