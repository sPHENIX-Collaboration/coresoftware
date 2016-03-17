#ifndef __MSG_CONTROL_H__
#define __MSG_CONTROL_H__

#include "msg_profile.h"
#include "event_io.h"
#include <string.h>

/** 
Objects of the msg\_control class maintain the profiling information 
for your messages. They add a few tokens to the output stream which can be
parsed by the msg\_buffer's routines, which so learns about the profile. 
Use it, for example, like this:
\begin{verbatim}
  msg_control *mvd_warning = new msg_control(MSG_TYPE_CODEDEBUG, 
					     MSG_SOURCE_MVD, 
					     MSG_SEV_WARNING);
  COUT <<  *mvd_info    << " this is a info message from MVD" << endl;
\end{verbatim}

The mvd\_info's "<<" operator adds now the profiling information to the message,
which is then parsed again and stripped off by our custom streambuf.

If we do not have a custom streambuf, the msg\_control objects maintain silence.
They have to be activated by one of our custom streambufs first.

* @VERSION: 2.0
* ejd add string for component source of message
*/

class msg_control
{
   friend class msg_buffer;

  friend OSTREAM& operator<< (OSTREAM& , msg_control &);

protected:

  int msg_type;
  int msg_source;
  int msg_severity;
  int storedseverity;
  char * msg_sourcecomponent;


  int activate();
  int deactivate();

public:


  msg_control(const int mtype = MSG_TYPE_DEFAULT
	      , const int source = MSG_SOURCE_DEFAULT
	      , const int severity = MSG_SEV_DEFAULT
	      , const char *sourcecomponent = "ONLINE");

  virtual ~msg_control();

  virtual void set_severity(const int severity); 
  virtual int get_severity() const { return msg_severity;};
  virtual void reset_severity() {msg_severity = storedseverity; };
  // set the message source id
  virtual void set_source(const int source) {msg_source = source; };
  virtual int  get_source() const { return msg_source; };
  virtual void set_sourcecomponent(const char * msgsourcecomponent = "ONLINE");
  virtual const char * get_sourcecomponent(){ return msg_sourcecomponent; };

  static int xx_active;

};

#endif /* __MSG_CONTROL_H__ */
