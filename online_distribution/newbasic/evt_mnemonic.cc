// ---------------------------------------------------------------------
#include "EventTypes.h"

const char *get_evt_mnemonic(const int id)
{

  int mid = id & ( CORRUPTEVENTMASK ^ 0xffff);
  switch (mid)
  {
     case (DATAEVENT):
     case (DATA2EVENT):
     case (DATA3EVENT):
        return "Data Event";

     case (SPILLONEVENT):
        return "Spill On Event";

     case (BEGRUNEVENT):
        return "Begin Run Event";

     case (RUNINFOEVENT):
        return "Run Info Event";

     case (ENDRUNEVENT):
        return "End Run Event";

     case (13):
        return "Movie Event";

     case (SCALEREVENT):
        return "Scaler Event";

     case (REJECTEDEVENT):
        return "Rejected Event";

     case (SPILLOFFEVENT):
        return "Spill Off Event";

     default:
        return "Unknown Event";
  }


}
