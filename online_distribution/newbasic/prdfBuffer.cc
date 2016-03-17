#include <prdfBuffer.h>
#include <EventTypes.h>
#include <A_Event.h>
#include <Cframe.h>
#include <stdio.h>

// the constructor first ----------------

prdfBuffer::prdfBuffer ()
{
  is_good =1;
}

prdfBuffer::prdfBuffer (PHDWORD *array , const int length )
{
  is_good =1;
  bptr =  (buffer_ptr) array;
  data_ptr = &(bptr->data[0]);
  max_length = length;
  current_index = 0;

  if (bptr->ID != BUFFERMARKER) 
  {
    //COUT << " will swap the buffer " << std::endl;
    unsigned int id = u4swap(bptr->ID);
    if (id != BUFFERMARKER) 
      {
	COUT << " wrong buffer" << std::endl;
	is_good =0;
	return;
      }
    if ( buffer_swap())
      {
	COUT << "problem in buffer swap" << std::endl;
	is_good = 0;
      }
  }
  buffer_size = array[0];
  //COUT << " done... " << std::endl;
 
}


prdfBuffer::~prdfBuffer()
{}

int prdfBuffer::buffer_swap()
{
 
  int i;
  unsigned int evtindex, frameindex;
  evtdata_ptr evtptr;
  PHDWORD *frameptr;

  //   swap the buffer header
  bptr->Length =  i4swap ( bptr->Length);
  bptr->ID =  i4swap ( bptr->ID);
  bptr->Bufseq =  i4swap ( bptr->Bufseq);
  bptr->Runnr =  i4swap ( bptr->Runnr);

  evtindex = 0;

  while (evtindex*4  < bptr->Length - BUFFERHEADERLENGTH)
    {
      //      COUT << "evt index " << evtindex << "  buffer length = " <<   bptr->Length << std::endl;

      // map event header on data
      evtptr = (  evtdata_ptr ) &bptr->data[evtindex];

      evtptr->evt_length = i4swap(evtptr->evt_length);
      evtptr->evt_type = i4swap(evtptr->evt_type);

      int written_length = (bptr->Length  + 8191)/8192;
      written_length *= 8192;

      if ( ( evtptr->evt_length + evtindex)*4 > written_length - BUFFERHEADERLENGTH  )
	{
	  std::cout << "Error: next event exceeds buffer length " << bptr->Length 
		    << " current index " << evtindex 
		    << "next event length " <<  evtptr->evt_length << std::endl;
 
	  // shit. that's wrong. terminate the buffer here.
	  bptr->data[evtindex] = 2;
	  bptr->data[evtindex] = 0;
	  return -1;
	}

      // see if we have got the end of buffer event
      if (evtptr->evt_length == 2 && evtptr->evt_type ==0) break;

      // swap the rest of the event header
      for (i=2; i<EVTHEADERLENGTH; i++) 
	bptr->data[evtindex+i] = i4swap(bptr->data[evtindex+i]);

      // mark first subevent
      frameindex = 0;

      int fl;
      
      while (frameindex < evtptr->evt_length - EVTHEADERLENGTH)
	{
	  // here we have a frame
	  frameptr =  &evtptr->data[frameindex];
	  i = frame_swap(frameptr, evtptr->evt_length - EVTHEADERLENGTH );   /* byte swap if wrong endianism */
	  if ( i )
	    {
	      evtptr->evt_type |= CORRUPTEVENTMASK;  
	      //  COUT << "about to call getFrameLength " << std::endl;
	      break;
	    }
	  fl =  getFrameLength(frameptr);
	  //  COUT << "framelength " << fl  << "  frameindex = " << frameindex << "  evt length = " <<  evtptr->evt_length << std::endl;

	  if (fl <=0 )
	    {
	      break;
	    }
	  frameindex += fl;
	}
      evtindex += evtptr->evt_length;
      if ( evtptr->evt_length ==0) 
	{
	  COUT << "0-length event found at " << evtindex << "  " << bptr->Length << std::endl; 
	  return -1;
	}
    }
  return 0;
}


int prdfBuffer::frame_swap(PHDWORD * fp, const int eventlength)
{
  int swapped_length = i4swap(*fp);
  if ( swapped_length > eventlength ||  swapped_length <0 )
    {
      //      std::cout << __FILE__ << "  " << __LINE__ << " corrupt frame length " << std::endl;
      return -1;
    }
  int i;
  for ( i = 0; i < swapped_length; i++)
    {
      fp[i] = i4swap(fp[i]);
    }

  return 0;
}


// ---------------------------------------------------------
Event * prdfBuffer::getEvent()

{

  if ( current_index < 0 ) return 0;
  if ( ! is_good ) return 0;
 
  Event *evt;
  evt =  0;

  // now is the new index pointing outside the allocated memory?
  if (current_index < 0 || current_index > max_length)
    {
      //COUT << "end of buffer r0 " << current_index << std::endl;
      current_index = -1;
      return evt;
    }


  int len = bptr->data[current_index];

  // maybe there is something wrong with it?
  if (len <= 0) 
    {
      //COUT << "end of buffer r1 " << std::endl;
      current_index = -1;
      return evt;
    }

  //  current_index += len;

  // are we pointing beyond the logical end of buffer?
  if (current_index >= (buffer_size/4) -6 )  //6 is 2 end-of-buffer + 4 buffer header
    {
      //COUT << "end of buffer r2 " << std::endl;
      current_index = -1;
      return evt;
    }

  // are we pointing to an end-of-buffer event? 
  if (bptr->data[current_index] == 2 && bptr->data[current_index+1] == 0)
    {
      //COUT << "end of buffer r3 " << std::endl;
      current_index = -1;
      return evt;
    }

  // none of the above, just return
  evt =  new A_Event( &bptr->data[current_index]);
  current_index += len;
  return evt;

}

// ---------------------------------------------------------
int * prdfBuffer::getEventData()

{

  if ( current_index < 0 ) return 0;
 

  int *evtData = 0;

  // now is the new index pointing outside the allocated memory?
  if (current_index < 0 || current_index > max_length)
    {
      //COUT << "end of buffer r1 " << current_index << std::endl;
      current_index = -1;
      return evtData;
    }


  int len = bptr->data[current_index];

  // maybe there is something wrong with it?
  if (len <= 0) 
    {
      current_index = -1;
      return evtData;
    }


  // are we pointing beyond the logical end of buffer?
  if (current_index >= (buffer_size/4) -6 )  //6 is 2 end-of-buffer + 4 buffer header
    {
      //COUT << "end of buffer r2 " << std::endl;
      current_index = -1;
      return evtData;
    }

  // are we pointing to an end-of-buffer event? 
  if (bptr->data[current_index] == 2 && bptr->data[current_index+1] == 0)
    {
      //COUT << "end of buffer r3 " << std::endl;
      current_index = -1;
      return evtData;
    }

  // none of the above, just return
  evtData = (int *) &bptr->data[current_index];
  current_index += len;
  return evtData;

}


