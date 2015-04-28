#include "oncsSub_idsis3300r.h"

oncsSub_idsis3300r::oncsSub_idsis3300r(subevtdata_ptr data)
  :oncsSub_idsis3300 (data)
{
  samples = 0;
  decoded = 0;

}
  
int *oncsSub_idsis3300r::decode ( int * p)
{
  if ( decoded) 
    {
      p = 0;
      return 0;
    }
  int c, i;

  decoded = 1;

  int *SubeventData = &SubeventHdr->data;

  samples = (*SubeventData) & 0xffff;
  wraparound = (*SubeventData) >> 16;

  for ( c = 0; c < 4; c++)
    {

      for ( i = 0; i< samples; i++)
	{
	  v[c*2].push_back ( (SubeventData[ c* samples + i+1] >> 16) & 0x3fff );
	  v[c*2+1].push_back ( SubeventData[ c* samples + i+1] & 0x3fff );
	}
    }

  p = 0;
  return 0;
}

int oncsSub_idsis3300r::iValue(const int ch ,const int s)
{
  if ( decoded == 0  ) decoded_data1 = decode(&data1_length);

  if ( ch < 0 || ch >7 ) return 0;
  if ( s < 0 || s >=samples ) return 0;

  return v[ch][s];

}

