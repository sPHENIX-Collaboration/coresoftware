/* 
** fakeFrame.C
** 
** Author: $Author: dave $  
**   Date: $Date: 2003/03/01 20:14:42 $ 
** 
** $Log: fakeFrame.C,v $
** Revision 1.2  2003/03/01 20:14:42  dave
** Put ifdef DEBUG around Debug_Output call instead of trying to define the call itself as a macro
**
** Revision 1.1.1.1  2000/07/21 01:51:12  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.6  1999/09/29 22:16:07  steinber
** mods to bring afs to nevis1 version
**
** Revision 1.4  1998/12/11 22:02:01  markacs
** (stephen markacs) adding log into cvs tags
** 
*/

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include "phenixOnline.h"
#include "fakeFrame.h"
#include "Cframe.h"
#include "CframeV1.h"
#include "packetRoutines.h"
#include "frameRoutines.h"
#include "framePublic.h"
#include "formatIO.h"

#define BUFFER_SIZE 65536

int fakeFrame(
              PHDWORD* frame_ptr,
              int bufferSize,
              UINT sourceId,
              ALIGNBLK alignBlk,
              int number_of_packets,
              int packet_ids[], 
              int packet_lengths[],
              int historySize,
              int errorSize
              )
{
  int remainingSpace;
  
  int total_frame_size;
  int maxFrameLen;
  PHDWORD* write_ptr;
  
  //  int total_data_words;
  PHDWORD data_arr[BUFFER_SIZE];
  UINT addr_arr[BUFFER_SIZE];
  int packetLength;
  int numBytes;
  
  int err;
  int i,j;
  int historyBlockSize;
  
  // Set registers and clear data buffer
  
  maxFrameLen = BUFFER_SIZE;
  dwordClear(&data_arr,BUFFER_SIZE);
  
  if ( (err = makeFrameHdr(frame_ptr,maxFrameLen,0x1,0x1, sourceId)) != 0 ){
    printf("FakeFrame: makeFrameHdr failed\n");
    return(-1);
  }
  
  setAlignBlock(frame_ptr,(PHDWORD*) &alignBlk,2);
  
  write_ptr = findFrameEnd(frame_ptr)+1;
  
  for (i = 0; i<number_of_packets;i++)
  {
    
    remainingSpace = bufferSize - (int) (write_ptr - frame_ptr);
    //printf("Remaining space = %d\n",remainingSpace);
    makeUnstructPacket(write_ptr,
      remainingSpace,
      packet_ids[i], //packetId
      sizeof(PHDWORD), //wordSize
      0  //hitformat
      );

//    numBytes = packet_lengths[i] * sizeof(PHDWORD);

    numBytes = packet_lengths[i] ; // number of "words", actually
    
    for (j=0; j<packet_lengths[i];j++){
//      data_arr[j]=packet_ids[i];
      data_arr[j]=j;
    }
    
    packetLength = storePacketHits(
      write_ptr,
      remainingSpace,
      addr_arr,
      (BYTE*) &data_arr,
      numBytes,
      0
      );
    extendFrameDataNopad(frame_ptr,bufferSize,packetLength);    
#ifdef DEBUG
    Debug_Output("Packet %2d : Length = %d Total=%d\n",
	     i,packetLength,getFrameDataLength(frame_ptr));
#endif
    write_ptr += packetLength;
  }
  
  //  printf("1: Frame data length: %d\n",getFrameDataLength(frame_ptr));
  remainingSpace = bufferSize - (int) (write_ptr - frame_ptr);
  historyBlockSize = storeFrameHistory(frame_ptr,
    remainingSpace,
    (PHDWORD*) &data_arr,
    historySize);
  //  printf("2: Frame data length: %d\n",getFrameDataLength(frame_ptr));
  write_ptr += historyBlockSize;
  
  /*
  printf("Frame lengths: Header,Data,History = %d, %d, %d\n",
	 getFrameHdrLength(frame_ptr),
   getFrameDataLength(frame_ptr),
   getFrameHistoryLength(frame_ptr)
   );
  */  
  
  total_frame_size = getFrameLength(frame_ptr);
  
  return total_frame_size;
  
}

int splitFakeFrame(
                   PHDWORD* frame_ptr,
                   int total_frame_size,
                   int length_of_buffer,
                   PHDWORD* start_of_buffer[]
                   )
{
  int i,j;
  int number_of_buffers;
  PHDWORD* write_ptr;
  PHDWORD* write_ptr_buf;
  write_ptr = frame_ptr;
  
  number_of_buffers = total_frame_size / length_of_buffer + 1;
  
  for (i=0;i<number_of_buffers;i++){
    start_of_buffer[i] = (PHDWORD*) malloc(4*length_of_buffer);
    write_ptr_buf = start_of_buffer[i];
    
    for (j=0 ; j<length_of_buffer ; j++){
      *write_ptr_buf = *write_ptr; 
      write_ptr++;
      write_ptr_buf++;
    }
  }
  
  /*
  printf("\nFrame header from start_buffer[0]\n");
  dumpFrameHdr(start_of_buffer[0]);
  printf("\n");
  */
  
  return number_of_buffers;
  
}



