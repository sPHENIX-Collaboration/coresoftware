#include "oncsSub_idfnalmwpcv2.h"

// for ntohl and ntohs
#include <arpa/inet.h>
#include <cstring>


using namespace std;

oncsSub_idfnalmwpcv2::oncsSub_idfnalmwpcv2(subevtdata_ptr data)
  :oncsSub_idfnalmwpc (data)
{
}
  
oncsSub_idfnalmwpcv2::~oncsSub_idfnalmwpcv2() {}


int *oncsSub_idfnalmwpcv2::decode ( int *nwout)
{

  if ( _decoded) return 0;

  _decoded = 1;

  int i;

  int *d = (int *) &SubeventHdr->data;  // here begins the payload

  tdcmask = d[0];
  
  // count the number of "on" bits = number of TDCs
  n_tdcs=0;
  for ( i = 0; i < MAXNROFTDCS; i++)
    {
      if ( (1<< i) & tdcmask ) n_tdcs++;
    }

  // unsigned short *xx = (unsigned short *)&d[1];
  // for ( i = 0; i < 300; i++)
  //   {
  //     cout << setw(5) << i << hex << setw(8) << ntohs(*xx++) << dec << endl;
  //   }
  // cout << endl;


  int length = ntohl ( d[1] );
  spillinfo.spillwords = length; // in 16bit units

  //we now subtract 2 to remove the header length
  length -= 2;

  //  dlength is the remaining payload length without the 
  //  length word in 16bit units - we subtract 8 bytes for the mask and length
  int dlength = 2 * ( getLength()-getPadding() ) -4 ; 

  // length is the length of the entire spill stucture without the length 
  // in 16bit words; we do a consistency check
  if (length  > dlength ) return 0;
  
  short *s = ( short *) &d[2];
  int pos = 0;
  
  short v;

  spillinfo.spillcounter = ntohs ( s[pos++]);
  
  v =  ntohs ( s[pos++]);
  spillinfo.year  = ( v>> 8 ) & 0xff;
  spillinfo.month =  v & 0xff;

  v =  ntohs ( s[pos++]);
  spillinfo.day = ( v>> 8 ) & 0xff;
  spillinfo.hour = v & 0xff;

  v =  ntohs ( s[pos++]);
  spillinfo.minute = ( v>> 8 ) & 0xff;
  spillinfo.second = v & 0xff;

  // this is a bad kludge; the next two words should 
  // hold the spill trigger count but are 0.
  // so I jump ahead and getthem from the first TDC spill header. Oh well.
  int x = ntohs ( s[11]);  
  pos++;
  x <<=16;
  spillinfo.triggercount = x + ntohs ( s[12]);
  pos++;

  spillinfo.TDC_Status_Bits = ntohs ( s[pos++]);
  spillinfo.Spill_Link_Status_Bits = ntohs ( s[pos++]);

  //  cout << " triggers = " << spillinfo.triggercount << endl;

  TDCEvent *te;


  // we now push, up front, one event structure on the vector
  // for each trigger so we know that we have one for each event
  for ( i = 0; i< spillinfo.triggercount; i++)
    {
      te = new TDCEvent;
      memset ( te, 0, sizeof(TDCEvent) );
      TDCEventVector.push_back ( te);
    }

  
  for ( i = 0; i< n_tdcs; i++)
    {
      int *x = ( int *) &s[pos];  // we get a 32bit value
      tsh[i].spillwords = ntohl ( *x);
      pos+=2;  // pos counts in shorts

      tsh[i].TDC = ntohs ( s[pos++]) -1;

      x= ( int *) &s[pos];  // we get a 32bit value
      tsh[i].spilltriggercount = ntohl ( *x);
      pos +=2;

      tsh[i].spillstatus = ntohs ( s[pos++]);
    }

  int old_trigger_number = -1;


  while (pos < length)
    {
      short spillheaderlength = ntohs ( s[pos]);
      short tdcindex = ntohs ( s[pos+1]) - 1; // the controller counts like fortran
      if ( tdcindex >= MAXNROFTDCS ) // sanity check
	{
	  cout << __LINE__ << " " << __FILE__ 
	       << " wrong TDC index " << tdcindex << " pos= " << pos << endl;
	}

      // cout << __LINE__ << " " << __FILE__ 
      //  	   << " TDC Event header at pos " << pos << endl;

      // for ( i = 0 ; i < 8; i++)
      //  	{
      //  	  cout << "  " << i << "  " << ntohs ( s[pos+i] ) << endl;
      // 	}

      pos++;

      int *x= ( int *) &s[pos+2];  // we get a 32bit value
      int t = ntohl ( *x);
      pos +=2;
      //      cout << "   trigger number " << v << endl;

      if (  t <0 || t >= spillinfo.triggercount )
	{
	  cout << __LINE__ << " " << __FILE__ 
	       << " wrong Trigger number " << v << " pos= " << pos << endl;
	  return 0;
	}

 
      // we are only interested in this TDC if it has data (the header length is 8) 
      // (and passes the sanity check) 
      if ( spillheaderlength > 8 && tdcindex < MAXNROFTDCS )
	{
	  // we are lookin ahead two words to get the 
	  // trigger number. If this is different from the one
	  // we had before, we are done with the old event, and start 
	  // a new one. (intially we see -1, so it's always different).
	  // as a reminder - w.r.t. "pos":
	  // +0 TDC#
	  // +1 event status
	  // +2 "ordinal trigger numberMSB"
	  // +3 "ordinal trigger numberLSB"
	  // +4 trigger type
	  // +5 800MHz event time stamp (12 bits)
	  // +6 event time stamp (100Mhz) upper 16
	  // +7 event time stamp (100Mhz) lower 16
	  // +8 here start the hits

	  if ( t != old_trigger_number )  // ok, new event/trigger, so we get the next TDCEvent
	    {
	      old_trigger_number = t;
	      te = TDCEventVector[t];
	      te->trigger_nr = t;
	      te->evt_timestamp = ntohs (s[pos+5]); 
	      te->local_ts_upper = ntohs (s[pos+6]); 
	      te->local_ts_lower = ntohs (s[pos+7]); 

	      // quick sanity check: bits 3...11 of 
	      // the timstamp must be equal to the 
	      // 9 lsb's of te->local_ts_lower

	      // cout << " time stamps " << hex 
	      // 	   << ( ( te->evt_timestamp>>3) &0x1ff ) << "  " 
	      // 	   << ( te->local_ts_lower & 0x1ff) << dec << endl;


	      // now we construct the long 35 bit 800MHz time stamp 
	      // cout << "time stamps "  << hex 
	      // 	   << te->evt_timestamp << " " 
	      // 	   <<  te->local_ts_upper << " " 
	      // 	   <<  te->local_ts_lower << dec << endl;

	      te->absolute_time = te->local_ts_upper;
	      te->absolute_time <<= 16;
	      te->absolute_time += te->local_ts_lower;
	      te->absolute_time <<= 3;
	      te->absolute_time += ( te->evt_timestamp &0x3);
	      // and we clear the TDCData structures
	      memset ( te->td, 0, MAXNROFTDCS * sizeof( TDCData) );
	      
	      // we are done with filling in the new event "header" structure
	      // now we fill in the TCDData structure. some of this 
	      // is redundant to the header info, that's why we kept "pos" frozen
	    }
	      
	  te->td[tdcindex].words = spillheaderlength;
	  te->td[tdcindex].TDC = tdcindex;
	  te->td[tdcindex].EventStatus =ntohs (s[pos+1]); 
	  te->td[tdcindex].triggertype =ntohs (s[pos+4]); 

	  // this may look redundant with the code above, but is actually
	  // per TDC info - the abve was per event. Should match, and we fill it in
	  // so the user can perform a sanity check
	  te->td[tdcindex].evt_timestamp  = ntohs (s[pos+5]); 
	  te->td[tdcindex].local_ts_upper = ntohs (s[pos+6]); 
	  te->td[tdcindex].local_ts_lower = ntohs (s[pos+7]); 

	  // now we construct the long 35 bit 800MHz time stamp 
	  te->td[tdcindex].absolute_time = te->td[tdcindex].local_ts_upper;
	  te->td[tdcindex].absolute_time <<= 16;
	  te->td[tdcindex].absolute_time += te->td[tdcindex].local_ts_lower;
	  te->td[tdcindex].absolute_time <<= 3;
	  te->td[tdcindex].absolute_time += ( te->td[tdcindex].evt_timestamp &0x7);

	  // at long last, the business end of this TDC - the TDC hit 
	  int xpos = pos+8;  // that's where the hits start
	  for ( i = 0; i < te->td[tdcindex].words -9; i++)  // so many hits, header is 9 long
	    {
	      te->td[tdcindex].hits++;  // this is the same as the TDCHitlist.size()

	      TDC_hit *th = new TDC_hit;
	      te->td[tdcindex].TDCHitlist.push_back ( th );

	      short h = ntohs ( s[xpos+i]);  // this is now the TDC word
	      th->wire = ( h >> 10) & 0x3f;  // upper bits are "wire"
	      th->timestamp = h & 0x3ff;     // lower 10 bits are "timestamp"
	      
	    }
	  pos += spillheaderlength-1;  // -1 because we stepped "pos" one position in already
	}
      else
	{
	  //  short tdc = ntohs(s[pos]);
	  //  cout << " skipping TDC " << tdc << " words = " << l << " pos = " << pos << endl;
	  pos+= 8;
	}
      

    }

  return 0;
}


