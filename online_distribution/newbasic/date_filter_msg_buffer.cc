
#include "date_filter_msg_buffer.h"

// in the constructor we allocate the specified amount of bytes for
// the output array.

date_filter_msg_buffer::date_filter_msg_buffer(const int msglen)
  : filter_msg_buffer(msglen)
{

}

date_filter_msg_buffer::date_filter_msg_buffer(
				      const int type_max, 
				      const int source_max, 
				      const int sev_max, 
				      const int msglen)
  : filter_msg_buffer(type_max, source_max,sev_max,msglen)
{}



// in the destructor, we restore the original streambuf for COUT.
date_filter_msg_buffer::~date_filter_msg_buffer()
{}

// this is the function which is called as a result of the << std::endl 
// and flush() operations

int date_filter_msg_buffer::sync ()
{ 

  msgProfile mp;
  int length, ix;


  char *x = format(&length, &mp);
  int ret = 0;

  if ( state[mp.type][mp.source][mp.severity] )
    {
      if (mp.type ) streambuf_add_date ( original_streambuf);
#if defined(LVL2_WINNT) || defined(STREAMBUF_NEW_IOSTREAM)
      for (ix =0; ix < length; ix++)
	original_streambuf->sputc( x[ix] );
      ret = original_streambuf->pubsync();
#else
      for (ix =0; ix < length; ix++)
	original_streambuf->overflow( x[ix] );
        ret = original_streambuf->sync();
#endif
    }
  delete [] x;

  pos = 0;
  
  return ret;
}







