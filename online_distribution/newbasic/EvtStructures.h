#ifndef __EVT_STRUCTURES_H
#define __EVT_STRUCTURES_H

#include "phenixTypes.h"

typedef struct evt_data
 {
   unsigned int evt_length;
   int evt_type;
   int evt_sequence;
   int run_number;
   int date;
   int time;
   unsigned int reserved[2];
   PHDWORD data[99999];
} *evtdata_ptr;


#endif
