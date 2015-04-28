//
// fileEventIterator   mlp 4/19/1997
//
// this iterator reads events froma data file. 


#include "fileEventiterator.h"
#include "oncsEventiterator.h"
#include <stddef.h>
#include <string.h>
//#include <sys/types.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "Cframe.h"
#include "framePackets.h"

// there are two similar constructors, one with just the
// filename, the other with an additional status value
// which is non-zero on return if anything goes wrong. 

//#ifndef LVL2_WINNT
#include <lzobuffer.h>
//#endif


fileEventiterator::~fileEventiterator()
{
  if ( legacy) delete legacy;

  if (fd) close (fd);
  if (thefilename != NULL) delete [] thefilename;
  if (bp != NULL ) delete [] bp;
  if (bptr != NULL ) delete bptr;
}  


fileEventiterator::fileEventiterator(const char *filename)
{
  legacy = 0;
  open_file ( filename);
}  

fileEventiterator::fileEventiterator(const char *filename, int &status)
{
  legacy = 0;
  status =  open_file ( filename);
}


int fileEventiterator::open_file(const char *filename)
{
  fd  = open (filename, O_RDONLY | O_LARGEFILE);
  bptr = 0;
  bp = 0;
  allocatedsize = 0;
  thefilename = NULL;
  events_so_far = 0;
  verbosity=0;
  if (fd > 0) 
    {
      PHDWORD cp[2048];
      int xc = read ( fd, cp, BUFFERBLOCKSIZE);
      if ( ! validFrameHdr( &cp[12]) )
	{
	  std::cout << "This doesn't look right" << std::endl;
	  close (fd);
	  int s;
	  legacy = new oncsEventiterator(filename, s);
	  return s;
	}

      lseek (fd, 0, SEEK_SET);

      thefilename = new char[strlen(filename)+1];
      strcpy (thefilename, filename);
      last_read_status = 0;
      current_index = 0;
      return 0;
    }
  else
    {
      last_read_status = 1;
    }
  return 1;

}



void fileEventiterator::identify (OSTREAM &os) const
{ 
  if ( legacy) 
    {
      legacy->identify(os);
      return;
    }
  os << "fileEventiterator reading from " << thefilename << std::endl;

};


const char * fileEventiterator::getCurrentFileName() const
{ 
  if ( legacy) 
    {
      return legacy->getCurrentFileName();
    }

  static char namestr[512];
  if ( thefilename == NULL)
    {
      return " ";
    }
  else
    {
      strcpy (namestr, thefilename);
      return namestr;
    }
};




char *  fileEventiterator::getIdTag () const
{ 
  if ( legacy) 
    {
      return legacy->getIdTag();
    }

  //  sprintf (idline, " -- fileEventiterator reading from %s", thefilename);
  return "fileEventiterator";
};



// and, finally, the only non-constructor member function to
// retrieve events from the iterator.

Event * fileEventiterator::getNextEvent()
{
  if ( legacy) 
    {
      return legacy->getNextEvent();
    }
  Event *evt = 0;

  // if we had a read error before, we just return
  if (last_read_status) return NULL;

  // see if we have a buffer to read
  if (bptr == 0) 
    {
      if ( (last_read_status = read_next_buffer()) !=0 )
	{
	  return NULL;
	}
    }

  while (last_read_status == 0)
    {
      if (bptr) evt =  bptr->getEvent();
      if (evt) 
	{
	  events_so_far++;
	  return evt;
	}
      last_read_status = read_next_buffer();
    }

  return NULL;

}

// -----------------------------------------------------
// this is a private function to read the next buffer
// if needed. 

int fileEventiterator::read_next_buffer()
{
  unsigned int ip = 8192;
 
  buffer_size = 0;
  
  if (bptr) 
    {
      if (verbosity >0)
	{
	  int ecount =  bp[2] & 0xffff;
	  int atpid = ( bp[2] >> 16) & 0xffff ;
	  std::cout << "Length: " <<  bp[0]
		    << " Atp id: " << atpid
		    << " Events in header: " <<  ecount
		    << " Events counted: " << events_so_far;
	  
	  if ( ecount != events_so_far )
	    {
	      std::cout << " ****";
	    }
	  std::cout << std::endl;
	}
      delete bptr;
      bptr = 0;
    }
  events_so_far = 0;
	
  // set the pointer to char to the destination buffer
  char *cp = (char *) initialbuffer;

  int xc;

  while (buffer_size == 0 )
    {  
      // read the first record
      xc = read ( fd, cp, BUFFERBLOCKSIZE);
      
      // error of EoF?
      if ( xc < BUFFERBLOCKSIZE  ) 
	{
	  //      COUT << "ferror" << std::endl;
	  return -1;
	}


      // get the length into a dedicated variable
      if (initialbuffer[1] == BUFFERMARKER || initialbuffer[1]== GZBUFFERMARKER ||  initialbuffer[1]== LZO1XBUFFERMARKER) 
	{
	  buffer_size = initialbuffer[0];
	}
      else
	{
	  unsigned int  marker = buffer::u4swap(initialbuffer[1]);
	  if (marker == BUFFERMARKER || marker == GZBUFFERMARKER || marker ==  LZO1XBUFFERMARKER )
	    {
	      buffer_size = buffer::u4swap(initialbuffer[0]);
	    }
	}
    }


  int i;
  if (bp) 
    {
      if  (buffer_size > allocatedsize*4)
	{

	  delete [] bp;
	  i = (buffer_size +BUFFERBLOCKSIZE-1) /BUFFERBLOCKSIZE;
	  allocatedsize = i * BUFFERBLOCKSIZE/4;
	  bp = new PHDWORD[allocatedsize];
	  //  std::cout << __FILE__ << "  " << __LINE__ << " new bp pointer is " << bp << "  length value "  << bp[-1]<< std::endl;
	}
    }
  else
    {
      i = (buffer_size +BUFFERBLOCKSIZE-1) /BUFFERBLOCKSIZE;
      allocatedsize = i * BUFFERBLOCKSIZE/4;
      bp = new PHDWORD[allocatedsize];

    }
  for (i = 0; i<BUFFERBLOCKSIZE/4; i++ ) bp[i] = initialbuffer[i];

  cp = (char *) bp;

  // and update the destination buffer pointer
  cp += BUFFERBLOCKSIZE;

  PHDWORD read_so_far =  BUFFERBLOCKSIZE;

  int errorinread=0;

  // now we read records until the whole buffer is read 
  while ( ip < buffer_size)
    {
      //      COUT << "ip is " << ip << std::endl;
      // read the next record
      xc = read ( fd, cp, BUFFERBLOCKSIZE);
      if ( xc < BUFFERBLOCKSIZE ) 
	{
	  COUT << "error in buffer, salvaging" << std::endl;
	  bp[0] = read_so_far; 
	  errorinread =1;
	  break;
	}

      // update the pointer and byte count
      cp += BUFFERBLOCKSIZE;
      ip += BUFFERBLOCKSIZE;
      read_so_far += BUFFERBLOCKSIZE;
    }

  // and initialize the current_index to be the first event
//#ifndef WIN32
  if ( ( initialbuffer[1]== GZBUFFERMARKER || 
       buffer::u4swap(initialbuffer[1])== GZBUFFERMARKER ||
       initialbuffer[1]== LZO1XBUFFERMARKER || 
       buffer::u4swap(initialbuffer[1])== LZO1XBUFFERMARKER )
       && errorinread  )
    {
      bptr = 0;
      return -3;
    }
  else if ( initialbuffer[1]== GZBUFFERMARKER || buffer::u4swap(initialbuffer[1])== GZBUFFERMARKER )
    {
      bptr = new gzbuffer(bp, allocatedsize );
    }

  else if ( initialbuffer[1]== LZO1XBUFFERMARKER || buffer::u4swap(initialbuffer[1])== LZO1XBUFFERMARKER )
    {
      bptr = new lzobuffer ( bp, allocatedsize );
    }

  else
    {

//#endif
      bptr = new buffer ( bp, allocatedsize );
//#ifndef WIN32
    }
//#endif
  return 0;
}

