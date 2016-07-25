#include "oncsBuffer.h"
#include "oncsStructures.h"
#include "oncsSubConstants.h"

// the constructor first ----------------
oncsBuffer::oncsBuffer (PHDWORD *array , const PHDWORD length )
{
  bptr =  (buffer_ptr) array;
  data_ptr = &(bptr->data[0]);
  max_length = length;
  current_index = 0;

  if (bptr->ID != ONCSBUFFERID && bptr->ID != PRDFBUFFERID) // PRDFBUFFERID is for legacy data 
  {
    //    COUT << " will swap the buffer " << std::endl;
    unsigned int id = i4swap(bptr->ID);
    if (id != ONCSBUFFERID && id != PRDFBUFFERID ) 
      {
	COUT << " wrong buffer" << std::endl;
	return;
      }
    if ( buffer_swap()) 
      {
	COUT << "problem in buffer swap" << std::endl;
      }
  }
  buffer_size = array[0];
  
}

int oncsBuffer::buffer_swap()
{

  int i;
  unsigned int evtindex, sevtindex;
  oncsevtdata_ptr evtptr;
  subevtdata_ptr sevtptr;

  //   swap the buffer header
  bptr->Length =  i4swap ( bptr->Length);
  bptr->ID =  i4swap ( bptr->ID);
  bptr->Bufseq =  i4swap ( bptr->Bufseq);
  bptr->Runnr =  i4swap ( bptr->Runnr);

  evtindex = 0;

  while (evtindex < bptr->Length - BUFFERHEADERLENGTH)
    {
      //  COUT << "evt index " << evtindex << std::endl;

      // map event header on data
      evtptr = (  oncsevtdata_ptr ) &bptr->data[evtindex];

      evtptr->evt_length = i4swap(evtptr->evt_length);
      evtptr->evt_type = i4swap(evtptr->evt_type);

      // see if we have got the end of buffer event
      if (evtptr->evt_length == 2 && evtptr->evt_type ==0) break;

      // swap the rest of the event header
      for (i=2; i<EVTHEADERLENGTH; i++) 
	bptr->data[evtindex+i] = i4swap(bptr->data[evtindex+i]);

      // mark first subevent
      sevtindex = 0;

      while (sevtindex < evtptr->evt_length - EVTHEADERLENGTH)
	{
	  // map subevent structure on data
	  sevtptr = (subevtdata_ptr) &evtptr->data[sevtindex];

	  // swap buffer header
	  sevtptr->sub_length = i4swap(sevtptr->sub_length);
	  sevtptr->sub_id = i2swap(sevtptr->sub_id);
	  sevtptr->sub_type = i2swap(sevtptr->sub_type);
	  sevtptr->sub_decoding = i2swap(sevtptr->sub_decoding);
	  sevtptr->sub_padding = i2swap(sevtptr->sub_padding);
	  sevtptr->reserved[0] = i2swap(sevtptr->reserved[0]);
	  sevtptr->reserved[1] = i2swap(sevtptr->reserved[1]);
      
	  // now swap the data depending on the type
	  int *p = &sevtptr->data;

	  switch (sevtptr->sub_type)
	    {
	    case 1: break;

	    case 2: 
	      for (i=0; i<sevtptr->sub_length - SEVTHEADERLENGTH; i++)
		{
		  *p = i22swap(*p);
		  p++;
		}
	      break;

	    case 4:
	      for (i=0; i<sevtptr->sub_length - SEVTHEADERLENGTH; i++)
		{
		  *p = i4swap(*p);
		  p++;
		}
	      break;
      
	    default: 
	      COUT << "unknown data type " << sevtptr->sub_type << std::endl;
	      break;
	    }
	  sevtindex += sevtptr->sub_length;
	}
      evtindex += evtptr->evt_length;
    }
  return 0;
}

// ---------------------------------------------------------
int oncsBuffer::i4swap(const int in)
{
  union 
  {
    int i4;
    char c[4];
  } i,o;

  i.i4 = in;
  o.c[0] = i.c[3];
  o.c[1] = i.c[2];
  o.c[2] = i.c[1];
  o.c[3] = i.c[0];
  return o.i4;
}
// ---------------------------------------------------------
int oncsBuffer::i22swap(const int in)
{
  union 
  {
    int i4;
    char c[4];
  } i,o;

  i.i4 = in;
  o.c[0] = i.c[1];
  o.c[1] = i.c[0];
  o.c[2] = i.c[3];
  o.c[3] = i.c[2];
  return o.i4;
}

// ---------------------------------------------------------
short oncsBuffer::i2swap(const  short in)
{
  union 
  {
    short i2;
    char c[2];
  } i,o;

  i.i2 = in;
  o.c[0] = i.c[1];
  o.c[1] = i.c[0];

  return o.i2;
}


// ---------------------------------------------------------
Event * oncsBuffer::getEvent()
{
  if ( current_index < 0 ) return 0;

  Event *evt;
  evt =  new oncsEvent( &bptr->data[current_index]);

  
  int l =  evt->getEvtLength();
  if ( l<= 0) return 0;
  current_index += evt->getEvtLength();

  // now is the new index pointing outside the allocated memory?
  //  if (current_index < 0 || current_index > BUFFERSIZE)
  if (current_index < 0 || current_index >= buffer_size)
    {
      //COUT << "end of buffer r1" << current_index << std::endl;
      current_index = -1;
      return evt;
    }

  // are we pointing beyond the logical end of buffer?
  if (current_index > buffer_size/4 )
    {
      //COUT << "end of buffer r2" << std::endl;
      current_index = -1;
      return evt;
    }

  // are we pointing to an end-of-buffer event? 
  if (bptr->data[current_index] == 2 && bptr->data[current_index+1] == 0)
    {
      //COUT << "end of buffer r3" << std::endl;
      current_index = -1;
      return evt;
    }

  // none of the above, just return
  return evt;

}


