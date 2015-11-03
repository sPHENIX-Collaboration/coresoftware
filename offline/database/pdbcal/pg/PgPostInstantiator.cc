#include "PgPostApplication.hh"
#include "PgPostBankManager.hh"
#include "RunToTimePg.hh"

namespace {


    int PgPostApp  = PgPostApplication::Register();

    int PgPostBank = PgPostBankManager::Register();

    int PgPostrtt  = RunToTimePg::Register();
} 

//_xd __xd;


