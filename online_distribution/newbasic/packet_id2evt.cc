#include "packet_id2evt.h"

Packet_id2evt::Packet_id2evt(PACKET_ptr data)
  : Packet_w2 (data){}
  
int *Packet_id2evt::decode ( int *nwout)
{
  int *p,*k;
  int olength;
  int temp[MAX_OUTLENGTH];
  int i;
  int dlength = getDataLength();

  int status = decode_id2evt( temp
			      ,(short *)  findPacketDataStart(packet) 
			      ,dlength
			      ,MAX_OUTLENGTH, &olength);

  if (status || olength<=0 ) return NULL;
 
  p = new int[olength];
  k = p;
  for (i =0; i<olength; i++) *k++ = temp[i];
  *nwout = olength;
  return p;
}








