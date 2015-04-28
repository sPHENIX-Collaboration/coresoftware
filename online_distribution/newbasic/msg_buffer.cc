
#include "msg_buffer.h"
#include <strnstr.h>

// in the constructor we allocate the specified amount of bytes for
// the output array.

msg_buffer::msg_buffer(const int msglen)
{
  maximum_position = msglen;
  oBuffer = new char[msglen];
  for (int i=0; i < msglen; i++)
    {
      oBuffer[i] = '\0';
    }

  pos =0; 

  m = new msg_control(0,0,0);
  m->activate();

}

// in the destructor, we deallocte the memory again
msg_buffer::~msg_buffer()
{
  delete [] oBuffer;
  if (m) 
    {
      m->deactivate();
      delete m;
    }
}

// this is the function which is called as a result of the << std::endl 
// and flush() operations
int msg_buffer::sync ()
{ 
  pos = 0;
  return 0;
}



// and this routine 
int msg_buffer::overflow (int ch)
{ 

  if (pos >= maximum_position) // oops, our buffer is too small...
    {
      int oldmaximum_position=maximum_position;
      maximum_position += MSG_EXTENSION_AMOUNT;  // we extend it
      char * n = new char[maximum_position];
      for (int ix =0; ix < oldmaximum_position; ix++) n[ix] = oBuffer[ix];
      for(int ix=oldmaximum_position;ix<maximum_position;ix++)n[ix]='\0';
      delete [] oBuffer;      // get rid of too small array
      oBuffer = n;            // and put new array in its place
    }

  oBuffer[pos++] = ch;

  return 0;
}

// and this routine 
char* msg_buffer::format(int * length,  msgProfile *mp)
{

  char *c, *s, *n,*src;
  int remainingLength = pos;
  int srcmsglen = 0;

  if ( ! (c = strnstr(oBuffer,remainingLength,"<|",2) ) )
    {
      mp->type      = MSG_TYPE_DEFAULT;
      mp->source    = MSG_SOURCE_DEFAULT;
      mp->severity  = MSG_SEV_DEFAULT;
      mp->reserved1 = 0;
      mp->reserved2 = 0;
      mp->reserved3 = 0;
      //     mp->sourcecomponent = new char[strlen("ONLINE")+1 ];
      //      strcpy(mp->sourcecomponent,"ONLINE");
      strcpy(mp->sourcecomponent,"ONLINE");
      *length = pos;
      c = oBuffer;
      n = new char[pos+1]; 
      s = n;
      for (int i=0; i<*length; i++) *s++ = *c++;
      n[pos] = '\0';
      return n;

    }
  
  remainingLength = pos - (int) ( c-oBuffer);
  // extract message type
  if ( ( s =  strnstr(c,remainingLength,"|A" , 2) ) )
    {
      sscanf (s+2, "%d", &mp->type);
    }
  else mp->type = MSG_TYPE_DEFAULT;

  // extract message source
  if ( ( s =  strnstr(c,remainingLength,"|B" , 2) ) )
    {
      sscanf (s+2, "%d", &mp->source);
    }
  else mp->source = MSG_SOURCE_DEFAULT;

  // extract message severity
  if ( ( s =  strnstr(c,remainingLength,"|C" , 2) ) )
    {
      sscanf (s+2, "%d", &mp->severity);
    }
  else mp->severity = MSG_SEV_DEFAULT;

  // extract message reserved1
  if ( ( s =  strnstr(c,remainingLength,"|D" , 2) ) )
    {
      sscanf (s+2, "%d", &mp->reserved1);
    }
  else mp->reserved1 = 0;

  // extract message reserved2
  if ( ( s =  strnstr(c,remainingLength,"|E" , 2) ) )
    {
      sscanf (s+2, "%d", &mp->reserved2);
    }
  else mp->reserved2 = 0;

  // extract message reserved3
  if ( ( s =  strnstr(c,remainingLength,"|F" , 2) ) )
    {
      sscanf (s+2, "%d", &mp->reserved3);
    }
  else mp->reserved3 = 0;

  // extract message sourcecomponent
  if ( ( s =  strnstr(c,remainingLength,"|G" , 2) ) )
    {
      src = strnstr(c,remainingLength,"|>",2);
      srcmsglen = (src - (s +2) );
      strncpy(mp->sourcecomponent,s +2, srcmsglen);
      mp->sourcecomponent[srcmsglen] = '\0';
      //      sscanf (s+2, "%s", &mp->sourcecomponent);
    }
  else 
    {
      // mp->sourcecomponent = new char[strlen("ONLINE")+1 ];
      strcpy(mp->sourcecomponent,"ONLINE");
    }

  remainingLength = pos;
  if ( ( c =  strnstr(oBuffer,remainingLength,"|>" , 2) ) )  
    {
      c+=2;
      *length = ( pos - (int)(c - oBuffer) ) ;
    }
  else
    {
      c = oBuffer;
      *length =  pos  ;
    }
  n = new char[pos+1]; 
  s = n;
  for (int i=0; i<*length; i++) *s++ = *c++;
  n[pos] = '\0';
  return n;
}

void streambuf_add_date (STREAMBUF * sb)
{

  
  struct timeb tp;
  ftime(&tp);
  char *timestr = ctime(&tp.time);
  int length = strlen(timestr);
  
  //  char timestr[128];
  //sprintf(timestr, "%d", time(0));
  //int length = strlen(timestr);

  char *cc = timestr;

#if defined(LVL2_WINNT) || defined(STREAMBUF_NEW_IOSTREAM)

  for ( int ix =0; ix < length-1; ix++) 
    sb->sputc( *cc++ );
  
  sb->sputc( ':' );
  sb->sputc( ' ' );

#else
  for ( int ix =0; ix < length-1; ix++) 
    sb->overflow( *cc++ );
  
  sb->overflow( ':' );
  sb->overflow( ' ' );
#endif
}

int printf (const char *format, ...)
{
  va_list ap;
  va_start (ap,format);
  char x[1000];
  vsprintf (x,format,ap);

  if (x[strlen(x)-1] == '\n') x[strlen(x)-1] = '\0';
  COUT << x << ENDL;
  return 0;

}

