#include "oncsSub_idfnalmwpcv2.h"

// for ntohl and ntohs
#include <arpa/inet.h>
#include <cstring>


using namespace std;

oncsSub_idfnalmwpcv2::oncsSub_idfnalmwpcv2(subevtdata_ptr data)
  :oncsSubevent_w4 (data)
{
  n_tdcs = 0;
  length= 0;
  _decoded = 0;

  // initialize all to 0
  memset (&spillinfo, 0 , sizeof(spillinfo) );
  memset (tsh, 0 , MAXNROFTDCS * sizeof(TDCspillheader) );

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


  int length = ntohl ( d[1] );
  spillinfo.spillwords = length; // in 16bit units
  //  cout << __FILE__ << " " << __LINE__ << "Spillwords " << length << endl;

  //we now subtract 2 to remove the header length
  length -= 2;

  //  dlength is the remaining payload length without the 
  //  length word in 16bit units - we subtract 8 bytes for the mask and length
  int dlength = 2 * ( getLength()-getPadding() ) -8 ; 

  // length is the length of the entire spill stucture without the length 
  // in 16bit words; we do a consistency check
  if (length  > dlength ) return 0;
  
  unsigned short *s = ( unsigned short *) &d[1];

  //pre-swap the payload
  for ( i = 0; i < length; i++)
    {
      *s = ntohs(*s);
      s++;
    }

  // s = ( unsigned short *) &d[1];
  // for ( i = 100; i < 200; i++)
  //   {
  //     cout << setw(4) << i << "  " << hex << s[i] << dec << endl;
  //   }
  

  
  //and reset the "s" pointer
  // we make it so that "pos" point to the position in the payload
  // so we can think in terms of the counts in teh docs
  s = ( unsigned short *) &d[1];
  int pos = 2;  // we have digested 0 and 1
  
  short v;

  spillinfo.spillcounter = s[pos++];
  //cout << __FILE__ << " " << __LINE__ << "spillinfo.spillcounter " << spillinfo.spillcounter << endl;
  
  v =  s[pos++];
  spillinfo.year  = ( v>> 8 ) & 0xff;
  spillinfo.month =  v & 0xff;

  v =  s[pos++];
  spillinfo.day = ( v>> 8 ) & 0xff;
  spillinfo.hour = v & 0xff;

  v =  s[pos++];
  spillinfo.minute = ( v>> 8 ) & 0xff;
  spillinfo.second = v & 0xff;

  // Spill trigger count  32bits
  int x = s[pos++];    //MSB
  x <<=16;
  spillinfo.triggercount = x + s[pos++];  //LSB
  
  //  cout << __FILE__ << " " << __LINE__ << " spillinfo.triggercount " << spillinfo.triggercount << endl;

  spillinfo.TDC_Status_Bits = s[pos++];
  spillinfo.Spill_Link_Status_Bits = s[pos++];

  //cout << __FILE__ << " " << __LINE__ << " spillinfo.TDC_Status_Bits " << spillinfo.TDC_Status_Bits << endl;
  //cout << __FILE__ << " " << __LINE__ << " spillinfo.Spill_Link_Status_Bits " << spillinfo.Spill_Link_Status_Bits << endl;


  TDCEvent *te;


  // we now push, up front, one event structure on the vector
  // for each trigger so we know that we have one for each event
  for ( i = 0; i< spillinfo.triggercount; i++)
    {
      te = new TDCEvent;
      memset ( te, 0, sizeof(TDCEvent) );
      TDCEventVector.push_back ( te);
    }

  // fill the WireChamberTdcSpillHeader (per TDC)  
  for ( i = 0; i< n_tdcs; i++)
    {

      //      cout << __FILE__ << " " <<  __LINE__ << " pos, value " << setw(4) << pos << " " << hex << s[pos] << dec << endl;
      tsh[i].spillwords =  s[pos++] << 16;  // we get a 32bit value
      tsh[i].spillwords +=  s[pos++] ;  
      //cout << __FILE__ << " " <<  __LINE__ << " i, spillwords, " << i << " " << tsh[i].spillwords << endl;

      tsh[i].TDC = s[pos++] -1;
      
      //cout << __FILE__ << " " <<  __LINE__ << " i, TDC " << i << " " << tsh[i].TDC << endl;

      tsh[i].spilltriggercount = pos++ << 16;
      tsh[i].spilltriggercount += pos++;
      //cout << __FILE__ << " " <<  __LINE__ << " spilltriggercount " <<  tsh[i].spilltriggercount << endl;

      tsh[i].spillstatus = s[pos++];
    }

  int old_trigger_number = -1;


  // here we keep the pos at the same value for a while
  while (pos < length)
    {
      //      cout << __FILE__ << " " <<  __LINE__ << " Event hdr: pos, value, length " << setw(4) << pos << " " << hex << s[pos] << dec << "  " << length << endl;
      int EventWordCount = s[pos];
      //      cout << __FILE__ << " " <<  __LINE__ << " EventWordCount " << EventWordCount << endl;

      short tdcindex =  s[pos+1] - 1; // the controller counts like fortran
      //cout << __FILE__ << " " <<  __LINE__ << " tdcindex " << tdcindex << endl;
      if ( tdcindex >= MAXNROFTDCS ) // sanity check
	{
	  cout << __LINE__ << " " << __FILE__ 
	       << " wrong TDC index " << tdcindex << " pos= " << pos-1 << endl;
	}

      int EventStatus = s[pos+2];

      //cout << __FILE__ << " " <<  __LINE__ << " pos, value " << setw(4) << pos << " " << hex << s[pos] << dec << endl;
      int SpillTriggerCount = s[pos+3]<<16;  // we get a 32bit value
      SpillTriggerCount += s[pos+4] -1;

      if (  SpillTriggerCount <0 || SpillTriggerCount >= spillinfo.triggercount )
	{
	  cout << __LINE__ << " " << __FILE__ 
	       << " wrong Trigger number " << SpillTriggerCount << " pos= " << pos << endl;
	  return 0;
	}

      //cout << __FILE__ << " " <<  __LINE__ << " SpillTriggerCount " << SpillTriggerCount << endl;
      
      if ( SpillTriggerCount != old_trigger_number )  // ok, new event/trigger, so we get the next TDCEvent
	{
	  old_trigger_number = SpillTriggerCount;
	  te = TDCEventVector[SpillTriggerCount];
	  te->trigger_nr = SpillTriggerCount;
	  te->evt_timestamp = s[pos+5]; 
	  te->local_ts_upper = s[pos+6]; 
	  te->local_ts_lower = s[pos+7]; 

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
	      
      te->td[tdcindex].words = EventWordCount;
      te->td[tdcindex].TDC = tdcindex;
      te->td[tdcindex].EventStatus =s[pos+2]; 
      te->td[tdcindex].triggertype =s[pos+5]; 

      // this may look redundant with the code above, but is actually
      // per TDC info - the abve was per event. Should match, and we fill it in
      // so the user can perform a sanity check
      te->td[tdcindex].evt_timestamp  = s[pos+6]; 
      te->td[tdcindex].local_ts_upper = s[pos+7]; 
      te->td[tdcindex].local_ts_lower = s[pos+8]; 

      // now we construct the long 35 bit 800MHz time stamp 
      te->td[tdcindex].absolute_time = te->td[tdcindex].local_ts_upper;
      te->td[tdcindex].absolute_time <<= 16;
      te->td[tdcindex].absolute_time += te->td[tdcindex].local_ts_lower;
      te->td[tdcindex].absolute_time <<= 3;
      te->td[tdcindex].absolute_time += ( te->td[tdcindex].evt_timestamp &0x7);

      // at long last, the business end of this TDC - the TDC hit 
      int xpos = pos+9;  // that's where the hits start
      for ( i = 0; i < te->td[tdcindex].words -9; i++)  // so many hits, header is 9 long
	{
	  te->td[tdcindex].hits++;  // this is the same as the TDCHitlist.size()
	  
	  TDC_hit *th = new TDC_hit;
	  te->td[tdcindex].TDCHitlist.push_back ( th );

	  short h = s[xpos+i];  // this is now the TDC word
	  th->wire = ( h >> 10) & 0x3f;  // upper bits are "wire"
	  th->timestamp = h & 0x3ff;     // lower 10 bits are "timestamp"
	      
	}
      pos += EventWordCount;  // -1 because we stepped "pos" one position in already
    }

  

  return 0;
}

int oncsSub_idfnalmwpcv2::iValue(const int ich,const char *what)
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


int oncsSub_idfnalmwpcv2::iValue(const int trigger ,const int tdc, const int index)
{
  decode (0);
  if ( trigger < 0 || (unsigned int ) trigger >= TDCEventVector.size()) return 0;
  if ( tdc < 0 || tdc >= n_tdcs) return 0;
  if ( index < 0 || index >= TDCEventVector[trigger]->td[tdc].hits) return 0;

  return TDCEventVector[trigger]->td[tdc].TDCHitlist[index]->wire;
}

int oncsSub_idfnalmwpcv2::iValue(const int trigger ,const int tdc, const char *what)
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

int oncsSub_idfnalmwpcv2::iValue(const int trigger ,const int tdc, const int index, const char *what)
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

void  oncsSub_idfnalmwpcv2::dump ( OSTREAM& os )
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

