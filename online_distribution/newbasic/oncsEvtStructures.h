#ifndef __ONCSEVT_STRUCTURES_H
#define __ONCSEVT_STRUCTURES_H

#include "oncsEvtConstants.h"

typedef struct oncsevt_data
 {
   unsigned int evt_length;
   int evt_type;
   int evt_sequence;
   int run_number;
   int date;
   int time;
   int reserved[2];
   int data[99999];
} *oncsevtdata_ptr;


#endif
