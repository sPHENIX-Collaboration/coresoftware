
#include "filter_msg_buffer.h"

// in the constructor we allocate the specified amount of bytes for
// the output array.

filter_msg_buffer::filter_msg_buffer(const int msglen)
  : msg_buffer(msglen)
{
  // in the constructor we replace COUT's streambuf with this one.
  // we remember the old one because we may want to replace it.

#ifdef RDBUF_ACCEPTS_STREAMBUF
  original_streambuf = COUT.rdbuf ( (STREAMBUF *) this);
#else
  original_streambuf = COUT.rdbuf ();
	COUT = (STREAMBUF *) this;
#endif

  int i,j,k;
  msg_type_max = MSG_TYPE_MAX;
  msg_source_max = MSG_SOURCE_MAX;
  msg_sev_max = MSG_SEV_MAX;

  // allocate new memory for the matrix
  state = new int **[msg_type_max];
  for (i=0; i< msg_type_max; i++)
    {
      state[i] = new int* [msg_source_max];
      for ( j =0; j<msg_source_max; j++) state[i][j] = new int [msg_sev_max];
    }

  for (i=0; i < msg_type_max; i++)
    for (j=0; j < msg_source_max ; j++)
      for (k=0; k < msg_sev_max ; k++)
	state[i][j][k] = ON;
}

filter_msg_buffer::filter_msg_buffer(const int type_max, const int source_max, 
		     const int sev_max, const int msglen)
  : msg_buffer(msglen)
{
  // in the constructor we replace COUT's streambuf with this one.
  // we remember the old one because we may want to replace it.
#ifdef RDBUF_ACCEPTS_STREAMBUF
  original_streambuf = COUT.rdbuf ( (STREAMBUF *) this);
#else
  original_streambuf = COUT.rdbuf ();
	COUT = (STREAMBUF *) this;
#endif

  int i,j,k;
  msg_type_max = type_max;
  msg_source_max = source_max;
  msg_sev_max = sev_max;

  // allocate new memory for the matrix
  state = new int **[msg_type_max];
  for (i=0; i< msg_type_max; i++)
    {
      state[i] = new int* [msg_source_max];
      for ( j =0; j<msg_source_max; j++) state[i][j] = new int [msg_sev_max];
    }

  for (i=0; i < msg_type_max; i++)
    for (j=0; j < msg_source_max ; j++)
      for (k=0; k < msg_sev_max ; k++)
	state[i][j][k] = ON;
}



// in the destructor, we restore the original streambuf for COUT.
filter_msg_buffer::~filter_msg_buffer()
{
#ifdef RDBUF_ACCEPTS_STREAMBUF
  COUT.rdbuf ( original_streambuf);
#else
  COUT= original_streambuf;
#endif

  // free memory from the matrix
  int i,j;
  for (i=0; i< msg_type_max; i++)
    {
      for ( j =0; j<msg_source_max; j++) delete [] state[i][j];
      delete [] state[i] ;
    }
  delete [] state;

}

// this is the function which is called as a result of the << std::endl 
// and flush() operations
int filter_msg_buffer::sync ()
{ 

  msgProfile mp;
  int length, ix;


  char *x = format(&length, &mp);
  int ret = 0;

  if ( state[mp.type][mp.source][mp.severity] )
    {
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

int filter_msg_buffer::set (const int type, 
	  const int source,
	  const int severity,
	  const int value)
{

  if ( type < 0 ||  type >= msg_type_max) return -1;  
  if ( source < 0 ||  source>= msg_source_max) return -2;
  if ( severity < 0 ||  severity>= msg_sev_max) return -3;



  state[type][source][severity] = value;
  return 0;
}

int filter_msg_buffer::set_severity_below_threshold(const int threshold,
						    const int value)
{
  if ( threshold < 0 ||  threshold >= msg_sev_max) return -1;

  int i,j,k;

	// leave type = 0 (UNSPECIFIED) untouched
  for (i=1; i < msg_type_max; i++)
    for (j=0; j < msg_source_max ; j++)
      for (k=0; k < threshold ; k++)
				state[i][j][k] = value;
  return 0;
}

int filter_msg_buffer::set_type (const int type,
				 const int value)
{
  if ( type < 0 ||  type >= msg_type_max) return -1;

  int j,k;
    for (j=0; j < msg_source_max ; j++)
      for (k=0; k < msg_sev_max ; k++)
	state[type][j][k] = value;
  return 0;
}

int filter_msg_buffer::set_source (const int source,
				 const int value)
{
  if ( source< 0 ||  source >= msg_type_max) return -1;

  int i,k;
  for (i=0; i < msg_type_max; i++)
      for (k=0; k < msg_sev_max ; k++)
	state [i] [source] [k] = value;
  return 0;
}


int filter_msg_buffer::all_off()
{
  int i,j,k;
  for (i=1; i < msg_type_max; i++)
    for (j=0; j < msg_source_max ; j++)
      for (k=0; k < msg_sev_max ; k++)
				state[i][j][k] = OFF;
  return 0;
}

int filter_msg_buffer::all_on()
{
  int i,j,k;
  for (i=0; i < msg_type_max; i++)
    for (j=0; j < msg_source_max ; j++)
      for (k=0; k < msg_sev_max ; k++)
	state[i][j][k] = ON;
  return 0;
}
