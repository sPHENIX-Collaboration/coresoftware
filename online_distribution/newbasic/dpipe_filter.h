#ifndef __DPIPE_FILTER_H__
#define __DPIPE_FILTER_H__

#include "Event.h"
#include "msg_control.h"

class DpipeFilter;

void dpipe_register(DpipeFilter * );
void dpipe_unregister(DpipeFilter *);

/** This is the pure virtual parent class for any dpipe filter

Upon loading, it registers itself with dpipe.

*/


class DpipeFilter
{

 public:
  DpipeFilter()
    {
      dpipe_register(this);
    }
  
  virtual ~DpipeFilter()
    {
      dpipe_unregister(this);
    }

  virtual int select (Event *e) = 0;

  virtual const char * idString() const =0; // {return "generic dpipe filter";}; 

  
 protected:
  
  
};



#endif /* __DPIPE_FILTER_H__ */
