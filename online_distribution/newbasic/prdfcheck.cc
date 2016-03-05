 

#include "buffer.h"
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <csignal>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>



void exitmsg()
{
  std::cout << "** usage: prdfcheck [-v] file name" << std::endl;
  exit(0);
}


int 
main(int argc, char *argv[])
{

  unsigned int buffer[8*8192];

  int fd;

  fd = open(argv[1], O_RDONLY | O_LARGEFILE);

  int needs_swap = 0;
  int length;
  int bufseq;
  int ip;

  int total_read = 0;

  int xc;

  xc = read ( fd, (char *)buffer, 8192);
  //  total_read++;
  while ( xc == 8192 )
    {
 
      ip = 8192;


      if ( buffer[1] == BUFFERMARKER || 
	   buffer::i4swap(buffer[1]) == BUFFERMARKER || 
	   buffer[1] == (int) GZBUFFERMARKER || 
	   buffer::i4swap(buffer[1]) == (int) GZBUFFERMARKER ||
	   buffer[1] == (int) LZO1XBUFFERMARKER || 
	   buffer::i4swap(buffer[1]) == (int) LZO1XBUFFERMARKER )
	{


	  if ( buffer::i4swap(buffer[1]) == BUFFERMARKER || 
	       buffer::i4swap(buffer[1]) == (int) GZBUFFERMARKER ||
	       buffer::i4swap(buffer[1]) == (int) LZO1XBUFFERMARKER )
	    {
	      needs_swap = 1;
	    }

	  if (needs_swap)
	    {
	      length = buffer::i4swap(buffer[0]);
	      bufseq = buffer::i4swap(buffer[2]);
	    }
	  else
	    {
	      length = buffer[0];
	      bufseq = buffer[2];
	    }

	  if ( needs_swap ) 
	    {
	      std::cout << "buffer at record " << std::setw(4) << total_read 
	       << " length = " << std::setw(7) <<  buffer::i4swap(buffer[0]) 
			<<  "  " << 	(buffer::i4swap(buffer[0]) + 8192)/8192	
			<< " marker = " << std::hex <<  buffer::i4swap(buffer[1]) << std::dec << "  ";
	    }
	  else
	    {
	      std::cout << "buffer at record " << std::setw(4) << total_read 
	       << " length = " << std::setw(7) << buffer[0] 
			<< "  " << 	( buffer[0] + 8192)/8192	
			<< " marker = " << std::hex << buffer[1] << std::dec << "  ";
	    }

	  if ( buffer[1] == BUFFERMARKER || 
	       buffer::i4swap(buffer[1]) == BUFFERMARKER ) std::cout << "Uncomp. Marker" << std::endl;

	  else if ( buffer[1] == (int) GZBUFFERMARKER || 
		    buffer::i4swap(buffer[1]) == (int) GZBUFFERMARKER ) std::cout << "GZIP Marker" << std::endl;
	  else if ( buffer[1] == (int) LZO1XBUFFERMARKER || 
		    buffer::i4swap(buffer[1]) == (int) LZO1XBUFFERMARKER ) 
	    {
	      std::cout << "LZO Marker ";
	      std::cout << " Or.length: " << buffer[3];
	      float ratio = 100.*buffer[0]/buffer[3];
	      std::cout << "  " << ratio << "%";

	      int e = buffer[2] & 0xffff;
	      int atp = (buffer[2] >> 16) & 0xffff;
	      if ( atp) 
		{
		  std::cout << " events: " << e << " from ATP " << atp << std::endl;
		}
	      else
		{
		  std::cout << std::endl;
		}
	    }


	  while ( ip < length)
	    {
	      xc = read ( fd, (char *)buffer, 8192);
	      if ( xc < 8192 ) 
		{
		  std::cout << "end or error in read loop at rec " << total_read << std::endl;
		  exit(1);
		}
	      total_read++;
	      ip+= 8192;
	    }
	  xc = read ( fd, (char *)buffer, 8192);
	  if ( xc < 8192 ) 
	    {
	      std::cout << "legitimate end or error at rec " << total_read << std::endl;
	      exit(1);
	    }
	  total_read++;

	}

      else
	{
	  if (needs_swap)
	    {
	      std::cout << "found a non-buffer start..."<< total_read 
			<< " length = " << buffer::i4swap(buffer[0]) << " marker = " << buffer::i4swap(buffer[1]) ;
	    }
	  else
	    {
	      std::cout << "found a non-buffer start..."<< total_read 
			<< " length = " << buffer[0] << " marker = " << buffer[1] ;
	    }

	  int skipped = 0;
	  while (  buffer[1] != BUFFERMARKER &&
	   buffer::i4swap(buffer[1]) != BUFFERMARKER &&
	   buffer[1] != (int) GZBUFFERMARKER && 
	   buffer::i4swap(buffer[1]) != (int) GZBUFFERMARKER )
	    {
	      xc = read ( fd, (char *)buffer, 8192);
	      if ( xc < 8192 ) 
		{
		  std::cout << "end or error in salavge loop at rec " << total_read << std::endl;
		  exit(1);
		}
	      total_read++;
	      skipped++;
	    }
	  std::cout << " Skipped " << skipped << std::endl;
	}
    }
  return 0;

}

