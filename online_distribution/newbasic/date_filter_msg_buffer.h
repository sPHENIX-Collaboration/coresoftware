#ifndef __DATE_FILTER_MSG_BUFFER_CC__
#define __DATE_FILTER_MSG_BUFFER_CC__



#include "filter_msg_buffer.h"

#define ON 1
#define OFF 0

/** This is the "date\_filter" msg\_buffer class which allows you to
filter messages based on their profile. Itacts much like the
filter\_msg\_buffer, but it also prepends a date tag to messages which have
a profile.
*/

class date_filter_msg_buffer : public   filter_msg_buffer {
  

public:

  /** The msglen parameter specifies the initial length of the 
      message string which is kept internally. If you exceed the length, 
      it is automautically extended.
  */
  date_filter_msg_buffer (const int msglen=256);

  /** This constructor defines a custom-sized matrix of 
      type/source/severities.  */
  date_filter_msg_buffer (const int type_max, const int source_max, 
		     const int sev_max, const int msglen=256);

  /// the virtual destructor
  virtual ~date_filter_msg_buffer();


  virtual int sync ();

};

#endif /* __DATE_FILTER_MSG_BUFFER_CC__ */



