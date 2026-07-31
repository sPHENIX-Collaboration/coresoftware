// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_SYNCDEFS_H
#define FFAOBJECTS_SYNCDEFS_H

#include <string>

namespace syncdefs
{
  // __attribute__((unused)) prevents compiler from flagging variable as unused
  static const int NUM_SYNC_VARS __attribute__((unused)) = 4;
// Sync object in old split TTree
//  static const char *SYNCVARS[] __attribute__((unused)) = {"eventcounter", "eventnumber", "runnumber", "segmentnumber"};
  static const char *SYNCVARS[] __attribute__((unused)) = {"DST#Sync"}; // TBranch name of sync object in unsplit TTree
  static const std::string SYNCNODENAME __attribute__((unused)) = "Sync";
}  // namespace syncdefs

#endif
