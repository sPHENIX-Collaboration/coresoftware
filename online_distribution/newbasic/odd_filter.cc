#include "dpipe_filter.h"
#include <iostream>


class odd_filter : public DpipeFilter {

public:


  int  select( Event *e) 
  {
    if ( e->getEvtSequence() &1)  return 1;
    return 0;
  };

 const char * idString() const 
  {
    return "odd events dpipe filter";
  };

};

odd_filter of;

