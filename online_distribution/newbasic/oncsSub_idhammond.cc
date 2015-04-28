#include "oncsSub_idhammond.h"

oncsSub_idhammond::oncsSub_idhammond(subevtdata_ptr data)
  :oncsSubevent_w4 (data){}
  
int *oncsSub_idhammond::decode ( int *nwout)
{
  int *p,*k;
  int olength;
  int temp[MAX_OUTLENGTH];
  int dlength = ( getLength()-4) - getPadding();
  int i;
  int *SubeventData = &SubeventHdr->data;

  int status = decode_idhammond( temp, SubeventData, dlength
			  ,MAX_OUTLENGTH, &olength);
  if (status || olength<=0 ) return NULL;
 
  p = new int[olength];
  k = p;
  for (i =0; i<olength; i++) *k++ = temp[i];
  *nwout = olength;
  return p;
}

int   oncsSub_idhammond::iValue(const int ich, const int iy)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return 0;
    }

  // protect us from negative indexes
  if (ich < 0 || iy <0) return 0;

  // let's calculate the proper index
  int channel = ich * 8 + iy;

  // see if our array is long enough
  if (channel > data1_length) return 0;

  return decoded_data1[channel];
}

float   oncsSub_idhammond::rValue(const int ich, const int iy)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return 0;
    }

  // protect us from negative indexes
  if (ich < 0 || iy <0) return 0;

  // let's calculate the proper index
  int channel = ich * 8 + iy;

  // see if our array is long enough
  if (channel > data1_length) return 0;

  return float(decoded_data1[channel]);
}








