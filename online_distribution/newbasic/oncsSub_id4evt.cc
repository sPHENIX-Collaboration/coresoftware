#include "oncsSub_id4evt.h"

oncsSub_id4evt::oncsSub_id4evt(subevtdata_ptr data)
  :oncsSubevent_w4 (data){}
  
int *oncsSub_id4evt::decode ( int *nwout)
{
  int *p,*k;
  int olength;
  int temp[MAX_OUTLENGTH];

  int dlength = ( getLength()-4) - getPadding();
  int i;
  int *SubeventData = &SubeventHdr->data;

  int status = decode_id4evt( temp, SubeventData, dlength
			  ,MAX_OUTLENGTH, &olength);

  if (status || olength<=0 ) return NULL;
 
  p = new int[olength];
  k = p;
  for (i =0; i<olength; i++) *k++ = temp[i];
  *nwout = olength;
  return p;
}






