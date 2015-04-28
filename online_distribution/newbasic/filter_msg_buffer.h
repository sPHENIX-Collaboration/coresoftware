#ifndef __FILTER_MSG_BUFFER_H__
#define __FILTER_MSG_BUFFER_H__



#include "msg_buffer.h"

#define ON 1
#define OFF 0

/** This is the "filter" msg\_buffer class which allows you to filter
messages based on their profile. Its default behavior is to let all messages pass. 
You can use the member functions to tailor the kind of messages filtered and passed on.

*/

class filter_msg_buffer : public   msg_buffer {

#ifndef __CINT__
protected:
  STREAMBUF * original_streambuf;
  int ***state;
  int msg_type_max;
  int msg_source_max;
  int msg_sev_max;
#endif


  //  int state[MSG_TYPE_MAX][MSG_SOURCE_MAX][MSG_SEV_MAX];

public:

  /** The msglen parameter specifies the initial length of the 
      message string which is kept internally. If you exceed the length, 
      it is automautically extended.
  */
  filter_msg_buffer (const int msglen=256);

  /** This constructor defines a custom-sized matrix of 
      type/source/severities.  */
  filter_msg_buffer (const int type_max, const int source_max, 
		     const int sev_max, const int msglen=256);

  /// the virtual destructor
  virtual ~filter_msg_buffer();

  /// the sync function overrides the streambuf's sync function

  // mlp -- the pubsync is what's needed for the new
  // iostream libraries - sync will no longer be a public method. 
  // for now we leave it as it was.


  virtual int sync ();


  /** the set function defines (ON or OFF) which messages are passed and
      which ones are filtered. This gives very fine-grained control. There 
      are other functions which give you more global control.
  */


  virtual int set (const int msg_type, 
	  const int msg_source,
	  const int msg_severity,
	  const int value = OFF);

  /// Globally set all messages below a certain severity threshold to ON or OFF 
  virtual int set_severity_below_threshold (const int threshold, 
					    const int value =OFF); 

  /// Globally switch ON or OFF all messages of a certain type
  virtual int set_type (const int type, 
			const int value =OFF); 

  /// Globally switch ON or OFF all messages from a certain source
  virtual int set_source (const int source, 
			const int value =OFF); 

  /// Globally switch all messages off
  virtual int all_off();

  /// Globally switch all messages on
  virtual int all_on();

};


#endif /* __FILTER_MSG_BUFFER_H__ */
