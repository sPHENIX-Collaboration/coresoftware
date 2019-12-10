#include "PgPostApplication.h"
#include "PgPostBankManager.h"
#include "RunToTimePg.h"

namespace
{
int PgPostApp = PgPostApplication::Register();

int PgPostBank = PgPostBankManager::Register();

int PgPostrtt = RunToTimePg::Register();
}  // namespace

//_xd __xd;
