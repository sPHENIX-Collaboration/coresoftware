#include "oncsSub_id2evt.h"

oncsSub_id2evt::oncsSub_id2evt(subevtdata_ptr data)
  :oncsSubevent_w2 (data){}


int *oncsSub_id2evt::decode ( int *nwout)
{
  int *p,*k;
  int olength;
  int temp[MAX_OUTLENGTH];
  int dlength = ( getLength()-4)*2 - getPadding();
  int i;
  short *SubeventData = (short * ) &SubeventHdr->data;

  int status = decode_id2evt( temp, SubeventData, dlength
			  ,MAX_OUTLENGTH, &olength);
  if (status) return NULL;

  p = new int[olength];
  k = p;
  for (i =0; i<olength; i++) *k++ = temp[i];
  *nwout = olength;
  return p;
}

