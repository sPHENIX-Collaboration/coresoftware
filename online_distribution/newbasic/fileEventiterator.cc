//
// fileEventIterator   mlp 4/19/1997
//
// this iterator reads events froma data file. 


#include <stddef.h>
#include <string.h>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

// there are two similar constructors, one with just the
// filename, the other with an additional status value
// which is non-zero on return if anything goes wrong. 
#include "fileEventiterator.h"

//#ifndef LVL2_WINNT
#include "lzobuffer.h"
//#endif


fileEventiterator::~fileEventiterator()
{
     if (fd) close (fd);
     if (thefilename != NULL) delete [] thefilename;
     if (bp != NULL ) delete [] bp;
     if (bptr != NULL ) delete bptr;
}  


fileEventiterator::fileEventiterator(const char *filename)
{
  open_file ( filename);
}  

fileEventiterator::fileEventiterator(const char *filename, int &status)
{
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
  _defunct = 0;

  if (fd > 0) 
    {
      thefilename = new char[strlen(filename)+1];
      strcpy (thefilename, filename);
      last_read_status = 0;
      current_index = 0;
      return 0;
    }
  else
    {
      last_read_status = 1;
      _defunct = 1;
    }
  return 1;

}



void fileEventiterator::identify (OSTREAM &os) const
{ 
  os << "fileEventiterator reading from " << thefilename;
  if ( _defunct ) os << " *** defunct";
  os<< std::endl;

};


const char * fileEventiterator::getCurrentFileName() const
{ 
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




const char *  fileEventiterator::getIdTag () const
{ 
  //  sprintf (idline, " -- fileEventiterator reading from %s", thefilename);
  return "fileEventiterator";
};



// and, finally, the only non-constructor member function to
// retrieve events from the iterator.

Event * fileEventiterator::getNextEvent()
{
  if ( _defunct ) return 0;
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
  PHDWORD initialbuffer[BUFFERBLOCKSIZE/4];

  unsigned int ip = 8192;
 
  buffer_size = 0;
  
  if (bptr) 
    {
      delete bptr;
      bptr = 0;
    }
  events_so_far = 0;
	
  // we are now reading the first block (8192 bytes) into
  // initialbuffer. The buffer is at least that long. We
  // look up the buffersize etc from that, and 
  // start filling the actual buffer. We do not read 
  // the data into the eventual destination since we may 
  // have to delete it to get more room (we might find that
  // the buffer is larger than what we have allocated).


  // set the pointer to char to the destination buffer
  char *cp = (char *) initialbuffer;

  int xc;


  // this while loop implements the skipping of 8k records until
  // we find a valid buffer marker. (We usually find it right away). 
 
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


      // get the buffer length into a dedicated variable
      if (initialbuffer[1] == BUFFERMARKER || initialbuffer[1]== GZBUFFERMARKER 
	  ||  initialbuffer[1]== LZO1XBUFFERMARKER || initialbuffer[1]== ONCSBUFFERMARKER) 
	{
	  buffer_size = initialbuffer[0];
	}
      else
	{
	  unsigned int  marker = buffer::u4swap(initialbuffer[1]);
	  if (marker == BUFFERMARKER || marker == GZBUFFERMARKER || marker ==  LZO1XBUFFERMARKER || marker == ONCSBUFFERMARKER)
	    {
	      buffer_size = buffer::u4swap(initialbuffer[0]);
	    }
	}
    }


  int i;
  if (bp) 
    {
      // this tests if we have enough space in the existing buffer
      if  (buffer_size > allocatedsize*4)   // no, we delete and make a bigger one
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

  //  for (i = 0; i<BUFFERBLOCKSIZE/4; i++ ) bp[i] = initialbuffer[i];
  memcpy ( bp, initialbuffer, BUFFERBLOCKSIZE);
 
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

  if ( ( initialbuffer[1]== GZBUFFERMARKER || 
       buffer::u4swap(initialbuffer[1])== GZBUFFERMARKER ||
       initialbuffer[1]== LZO1XBUFFERMARKER || 
       buffer::u4swap(initialbuffer[1])== LZO1XBUFFERMARKER )
       && errorinread  )
    {
      bptr = 0;
      return -3;
    }

  return buffer::makeBuffer( bp, allocatedsize, &bptr);


}

