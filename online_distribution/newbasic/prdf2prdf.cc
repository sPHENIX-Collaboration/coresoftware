
/*
**  readFrame.C
**
**    Reads a frame, accesses packets and fills structures 
**    like the drift chamber STAF table dDchDCM
**
*/
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <stdio.h>

#include "oBuffer.h"
#include "EventTypes.h"

#include "phenixOnline.h"

#include "Cframe.h"

#define MAXSIZE 8192
#define MAXBUFFERSIZE 256*8192


int main (int argc, char** argv) 
{
  int runnumber;

  char *filename;
  if (argc>=4) filename=argv[1];
  else 
    {
      COUT << "usage: " << argv[0] << " DATAFILE outfile run-number [frames to combine]"<< std::endl; 
      return 1; 
    }

  int fd =  open(filename, O_RDONLY | O_LARGEFILE);
  if (fd < 0  ) 
    {
      COUT << "ERROR: could not open " << filename << std::endl; 
      return 1;
    }
  unlink (argv[2]);
  int outfile = open(argv[2], O_WRONLY | O_CREAT | O_EXCL | O_LARGEFILE , 
		  S_IRWXU | S_IROTH | S_IRGRP ); 
  if (outfile < 0 ) 
    {
      COUT << "ERROR: could not open " << argv[2] << std::endl; 
      return 1;
    }

  // read the run number
  sscanf(argv[3], "%d", &runnumber);
  COUT << "will write run number " << runnumber << std::endl; 

   // read the run number
  int frames_to_combine = 1;
  if (argc==5)
    {
      sscanf(argv[4], "%d", &frames_to_combine);
      COUT << "will combine  " << frames_to_combine << " frames for each event" <<std::endl; 
    }
 static PHDWORD databuffer[MAXBUFFERSIZE];
  oBuffer *ob = new oBuffer(outfile, databuffer, MAXBUFFERSIZE, runnumber);

  // add the begin-run event
  ob->nextEvent(100, BEGRUNEVENT);

  PHDWORD frame_ptr[MAXSIZE];

  int eventnumber = 0;
  // read the first woird of the frame to see how long it is
  int nframe = 0;

  int nread;
  nread = read(fd, frame_ptr,4*6);    /* read frame from file into memory */ 

  int framelen;
  if (! checkFrameEndianism(frame_ptr) ) framelen = singleDwordByteSwap(frame_ptr[0]
);
  else framelen = frame_ptr[0];

  // now read the rest...
  nread = read(fd, &frame_ptr[6], 4*(framelen-6));
  if (! checkFrameEndianism(frame_ptr) ) byteSwapFrame(frame_ptr);

  ob->nextEvent(MAXBUFFERSIZE-100, DATAEVENT);
  eventnumber++;
  int reallength = 0;

  while ( nread > 0 )
    {
      // and add the one and only frame as a whole)
      ob->addFrame(frame_ptr);
      reallength += frame_ptr[0];
      nframe++;
		  
      // * now start to read in the next frame
		  
      nread = read(fd, frame_ptr,4*6);    /* read frame from file into memory */ 
      
      if (! checkFrameEndianism(frame_ptr) ) framelen = singleDwordByteSwap(frame_ptr[0]);
      else framelen = frame_ptr[0];

      // now read the rest...
      nread = read (fd, &frame_ptr[6], 4*(framelen-6));
      if (! checkFrameEndianism(frame_ptr) ) byteSwapFrame(frame_ptr);
      if ( nframe >=frames_to_combine)
	{
	  ob->nextEvent(reallength +  100, DATAEVENT);
	  eventnumber++;
	  reallength = 0;
	  nframe = 0;
	}
    }

  // add the end-run event
  ob->nextEvent(100, ENDRUNEVENT);

  COUT << eventnumber << " data events written " 
       << std::endl;


  // delete the oBuffer to wrap things up
  delete ob;

  return 0;
}


