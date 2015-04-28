#include "oncsSub_idsam.h"

oncsSub_idsam::oncsSub_idsam(subevtdata_ptr data)
  : oncsSubevent_w4 (data){}
  
int *oncsSub_idsam::decode ( int *nwout)
{
  int *p,*k;
  int olength;
  int temp[MAX_OUTLENGTH];
  int i;
  int dlength = ( getLength()-4) - getPadding();

  int *SubeventData = &SubeventHdr->data;

  int status = decode_idsam( temp, SubeventData, dlength
			  ,MAX_OUTLENGTH, &olength);

  if (status || olength<=0 ) return NULL;
 
  p = new int[olength];
  k = p;
  for (i =0; i<olength; i++) *k++ = temp[i];
  *nwout = olength;
  return p;
}









