#include "oncsSub_idfnalmwpc.h"
#include <cstring>

// for ntohl and ntohs
#include <arpa/inet.h>


using namespace std;

oncsSub_idfnalmwpc::oncsSub_idfnalmwpc(subevtdata_ptr data)
  :oncsSubevent_w4 (data)
{
  n_tdcs = 0;
  length= 0;
  _decoded = 0;
  
  // initialize all to 0
  memset (&spillinfo, 0 , sizeof(spillinfo) );
  memset (tsh, 0 , MAXNROFTDCS * sizeof(TDCspillheader) );
  

}
  

int *oncsSub_idfnalmwpc::decode ( int *nwout)
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

  spillinfo.triggercount = ntohs ( s[pos++]);
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
      tsh[i].spillwords = ntohs ( s[pos++]);
      tsh[i].TDC = ntohs ( s[pos++]) -1;
      tsh[i].spilltriggercount = ntohs ( s[pos++]);
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

      v = ntohs ( s[pos +2] ) -1;
      //      cout << "   trigger number " << v << endl;

      if (  v <0 || v >= spillinfo.triggercount )
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
	  // +2 "ordinal trigger number"
	  // +3 trigger type
	  // +4 800MHz event time stamp (12 bits)
	  // +5 event time stamp (100Mhz) upper 16
	  // +6 event time stamp (100Mhz) lower 16
	  // +7 here start the hits

	  if ( v != old_trigger_number )  // ok, new event/trigger, so we get the next TDCEvent
	    {
	      old_trigger_number = v;
	      te = TDCEventVector[v];
	      te->trigger_nr = v;
	      te->evt_timestamp = ntohs (s[pos+4]); 
	      te->local_ts_upper = ntohs (s[pos+5]); 
	      te->local_ts_lower = ntohs (s[pos+6]); 

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
	  te->td[tdcindex].triggertype =ntohs (s[pos+3]); 

	  // this may look redundant with the code above, but is actually
	  // per TDC info - the abve was per event. Should match, and we fill it in
	  // so the user can perform a sanity check
	  te->td[tdcindex].evt_timestamp  = ntohs (s[pos+4]); 
	  te->td[tdcindex].local_ts_upper = ntohs (s[pos+5]); 
	  te->td[tdcindex].local_ts_lower = ntohs (s[pos+6]); 

	  // now we construct the long 35 bit 800MHz time stamp 
	  te->td[tdcindex].absolute_time = te->td[tdcindex].local_ts_upper;
	  te->td[tdcindex].absolute_time <<= 16;
	  te->td[tdcindex].absolute_time += te->td[tdcindex].local_ts_lower;
	  te->td[tdcindex].absolute_time <<= 3;
	  te->td[tdcindex].absolute_time += ( te->td[tdcindex].evt_timestamp &0x7);

	  // at long last, the business end of this TDC - the TDC hit 
	  int xpos = pos+7;  // that's where the hits start
	  for ( i = 0; i < te->td[tdcindex].words -8; i++)  // so many hits, header is 8 long
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
	  pos+= 7;
	}
      

    }

  return 0;
}


int oncsSub_idfnalmwpc::iValue(const int ich,const char *what)
{

  decode (0);
  if ( strcmp(what,"NTDCS") == 0 )
  {
    return n_tdcs;
  }

  else if ( strcmp(what,"TRIGGERCOUNT") == 0 )
  {
    return spillinfo.triggercount;
  }

  else if ( strcmp(what,"SPILLCOUNTER") == 0 )
  {
    return spillinfo.spillcounter;
  }

  else if ( strcmp(what,"TDCSTATUSBITS") == 0 )
  {
    return   spillinfo.TDC_Status_Bits;
  }

  else if ( strcmp(what,"LINKSTATUSBITS") == 0 )
  {
    return spillinfo.Spill_Link_Status_Bits;
  }

  else if ( strcmp(what,"MONTH") == 0 )
  {
    return spillinfo.month;
  }

  else if ( strcmp(what,"DAY") == 0 )
  {
    return spillinfo.day;
  }

  else if ( strcmp(what,"YEAR") == 0 )
  {
    return spillinfo.year;
  }

  else if ( strcmp(what,"HOUR") == 0 )
  {
    return spillinfo.hour;
  }

  else if ( strcmp(what,"MINUTE") == 0 )
  {
    return spillinfo.minute;
  }

  else if ( strcmp(what,"SECOND") == 0 )
  {
    return spillinfo.second;
  }

  return 0;

}


int oncsSub_idfnalmwpc::iValue(const int trigger ,const int tdc, const int index)
{
  decode (0);
  if ( trigger < 0 || (unsigned int ) trigger >= TDCEventVector.size()) return 0;
  if ( tdc < 0 || tdc >= n_tdcs) return 0;
  if ( index < 0 || index >= TDCEventVector[trigger]->td[tdc].hits) return 0;

  return TDCEventVector[trigger]->td[tdc].TDCHitlist[index]->wire;
}

int oncsSub_idfnalmwpc::iValue(const int trigger ,const int tdc, const char *what)
{
  decode (0);
  if ( trigger < 0 || (unsigned int )trigger >= TDCEventVector.size()) return 0;
  if ( tdc < 0 || tdc >= n_tdcs) return 0;


  if ( strcmp(what,"HITS") == 0 )
  {
    if ( (unsigned int ) trigger >= TDCEventVector.size() )
      {
	cout << __LINE__ << " " << __FILE__ << " wrong trigger " <<  trigger << " max " << TDCEventVector.size() << endl;
	return 0;
      }

    return TDCEventVector[trigger]->td[tdc].hits;
  }
  else if ( strcmp(what,"TDC") == 0 )
  {
    return TDCEventVector[trigger]->td[tdc].TDC;
  }

  return 0;
}

int oncsSub_idfnalmwpc::iValue(const int trigger ,const int tdc, const int index, const char *what)
{
  decode (0);
  if ( trigger < 0 || (unsigned int )trigger >= TDCEventVector.size()) return 0;
  if ( tdc < 0 || tdc >= n_tdcs) return 0;
  if ( index < 0 || index >= TDCEventVector[trigger]->td[tdc].hits) return 0;


  if ( strcmp(what,"WIRE") == 0 )
  {
    return iValue(trigger, tdc, index);
  }

  else if ( strcmp(what,"TIMESTAMP") == 0 )
  {
    return TDCEventVector[trigger]->td[tdc].TDCHitlist[index]->timestamp;
  }

  
  return 0;
}



void  oncsSub_idfnalmwpc::dump ( OSTREAM& os )
{

  decode (0);
  identify(os);

  int i, j, k, l;


  os << " Date         " << iValue(0, "MONTH") 
     << "/"              << iValue(0,"DAY") 
     << "/20"            << iValue(0,"YEAR")  
     << " "              << iValue(0,"HOUR")
     << ":"              << iValue(0,"MINUTE") 
     << ":"              << iValue(0,"SECOND") << endl;
  
  os << " Spill Counter " << iValue(0,"SPILLCOUNTER") << endl;
  
  int tc = iValue(0,"TRIGGERCOUNT");

  os << " Trigger Count " << tc  << endl;
  os << " TDC Status Bits " <<  hex << iValue(0,"TDCSTATUSBITS") << dec << endl;
  os << " Link Status bits " << hex << iValue(0,"LINKSTATUSBITS") << dec << endl;

  for ( k = 0; k < tc; k++)
    {

      os << " new event " 
	 << "      trg#     tdc     hit    wire  timestamp" << endl;

      for ( l = 0; l < iValue(0,"NTDCS"); l++)
	{
	  int hits =iValue ( k, l, "HITS"); 
	  if ( hits)
	    {
	      int tdc = iValue(k,l,"TDC");

	      for ( j = 0; j < hits; j++)
		{
		  os << "              " << setw(6) << k 
		     << setw(8) <<  tdc  << setw(8) << j 
		     <<  setw (8)  << iValue(k,l,j)
		     << setw(8)  << iValue(k,l,j,"TIMESTAMP")  << endl;
		}
	    }
	}
    }
  os << endl;
  
  
  
}

oncsSub_idfnalmwpc::~oncsSub_idfnalmwpc()
{
  
  int i, j, k, l;
  TDCEvent *te;
  TDC_hit *th;

  for ( k = 0; k< TDCEventVector.size() ; k++)
    {
      te = TDCEventVector[k];

      for ( l = 0; l < n_tdcs; l++)
	{
	  for  ( i = 0; i < te->td[l].TDCHitlist.size() ; i++)
	    {
	      delete te->td[l].TDCHitlist[i];
	    }
	}
      delete te;
    }
}

