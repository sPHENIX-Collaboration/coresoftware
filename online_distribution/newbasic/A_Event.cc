#include "A_Event.h"

#include "EventTypes.h"

#include "packet_collection.h"

#include "packet.h"
#include "Cframe.h"
#include "framePackets.h"
#include "dataBlock.h"

#include "frameRoutines.h"
#include "frameHdr.h"


// the constructor first ----------------
A_Event::A_Event (PHDWORD *data)
{ 
  // set the framlist pointer to 0 just in case.
  framelist = 0;


  hasMap = 0;
  // we always make a pointer-based event.
  is_data_type = 0;
  errorcode = 0;

  // cast the pointer to the EventData pointer
  EventData = (evtdata_ptr) data;
  if ( getEvtType() & CORRUPTEVENTMASK)
    {
      errorcode=1;
    }
  else
    {
      errorcode = updateFramelist();
    }
}

A_Event::A_Event (int *data)
{ 
  // set the framlist pointer to 0 just in case.
  framelist = 0;

  hasMap = 0;
  // we always make a pointer-based event.
  errorcode = 0;
  is_data_type = 0;
  
  // cast the pointer to the EventData pointer
  EventData = (evtdata_ptr) data;
  if ( getEvtType() & CORRUPTEVENTMASK)
    {
      errorcode=1;
    }
  else
    {
      errorcode = updateFramelist();
    }

}

A_Event::~A_Event()
{
  delete [] framelist;
  if (is_data_type) delete [] (PHDWORD *) EventData;
  pmap.clear();

}


// the info-type calls
unsigned int 
A_Event::getEvtLength()
{
  return EventData->evt_length;
}

int 
A_Event::getEvtType()
{
  return EventData->evt_type;
}

int 
A_Event::getEvtSequence()
{
  return EventData->evt_sequence;
}

int 
A_Event::getRunNumber()
{
  return EventData->run_number;
}

#ifdef WIN32
const ULONGLONG  ticOffset = 35067168000000000UL;
const ULONGLONG  ticFactor = 10000000;
#else
const  unsigned long long  ticOffset = 35067168000000000ULL;
const   unsigned long long ticFactor = 10000000;
#endif


time_t A_Event::getTime() const 
{
  if ( EventData->time == -1)
     {
       return (time_t) EventData->date;
     }
  else
    {

      unsigned long long  x,y,z;
      x = (unsigned int)  EventData->time;
      y = (unsigned int)  EventData->date;
      z = y | ( x<<32 );
      time_t t =  (z - ticOffset) / ticFactor;
      return t;
    }
}


void A_Event::listFrame ( const int id, OSTREAM &os) const
{
  PHDWORD *fp;
  PHDWORD *pp;
  int i = 0;
  int j = 0;

  if (id == 0) 
    {
      
      while ( ( fp = framelist[i++]) )
	{
	  dumpFrame(fp, os);
	  //os << "Frame " << ++j << "align length: "<< getFrameAlignLength(fp) << std::endl;
	}
    }
  else // find the frame for a given packet
    {
      
      while ( (fp = framelist[i++]) )
	{
	  j++;
	  if ( ( pp = findFramePacketId (fp, id) ) !=  ptrFailure) 
	    {
	      dumpFrame(fp, os);
	      //os << "Frame " << j << "align length: "<< getFrameAlignLength(fp) << std::endl;
	      return;
	    }
	}
    }
}

void
A_Event::listHistory(const int id, OSTREAM &os) const
{
  os << "History Block: " << std::endl;

  int i = 0;
  if (id == 0) 
    {      
      while ( PHDWORD* fp = framelist[i++] ) 
	{
	  PHDWORD* h_ptr = findFrameHistoryStart(fp);
	  UINT len = getFrameHistoryLength(fp);
	  dumpBlock(h_ptr, len, os);
	}
    }
  else // find the frame for a given packet
    {      
      int j = 0;
      while ( PHDWORD* fp = framelist[i++] )
	{
	  j++;
	  if ( findFramePacketId(fp, id) != ptrFailure ) 
	    {
	      PHDWORD* h_ptr = findFrameHistoryStart(fp);
	      UINT len = getFrameHistoryLength(fp);
	      dumpBlock(h_ptr, len, os);
	      return;
	    }
	}
    }
  return;
}

void
A_Event::listError( const int id, OSTREAM& os ) const
{
  os << "Error Block: " << std::endl;
  int i = 0;
  if (id == 0) 
    {      
      while ( PHDWORD* fp = framelist[i++] ) 
	{
	  PHDWORD* ptr = findFrameErrorStart(fp);
	  UINT len = getFrameErrorLength(fp);
	  dumpBlock(ptr, len, os);
	  dumpErrorBlock(fp,os);
	}
    }
  else // find the frame for a given packet
    {      
      int j = 0;
      while ( PHDWORD* fp = framelist[i++] )
	{
	  j++;
	  if ( findFramePacketId(fp, id) != ptrFailure ) 
	    {
	      PHDWORD* ptr = findFrameErrorStart(fp);
	      UINT len = getFrameErrorLength(fp);
	      dumpBlock(ptr, len, os);
	      return;
	    }
	}
    }
}

void A_Event::dumpFrame(PHDWORD *fp, OSTREAM &os)
{
  // DLW: for SEQ number and code, there is no corresponding routine in Cframe.h, 
  // so I directly unpack the field.  This will break if the frame format ever changes...
  //
  os << "Frame length:       " << std::dec << getFrameLength(fp) << std::endl;
  os << "Frame mark:         " << std::hex << getFrameMark(fp) << std::dec << std::endl;
  os << "Frame Hdr version:  " << getFrameHdrVersion(fp) << std::endl;
  os << "Frame Hdr length:   " << getFrameHdrLength(fp) << std::endl ;
  os << "Frame Status:       " << getFrameStatus(fp) << std::endl;
  os << "Frame Seq Number:   " << (((*(fp+3))&0xff000000)>>24) << std::endl;
  os << "Frame Seq Code:     " << (((*(fp+3))&0x00ff0000)>>24) << std::endl;
  os << "Frame Source id:    " << getFrameSourceId(fp) << std::endl;
  os << "Frame data type:    " << getFrameDataType(fp) << std::endl;
  os << "Frame type:         " << getFrameType(fp) << std::endl;
  os << "Frame Error Length: " << getFrameErrorLength(fp) << std::endl;
  os << "Frame Hist Length:  " << getFrameHistoryLength(fp) << std::endl;
  os << "Frame align length: " << getFrameAlignLength(fp) << std::endl;
  os << "Frame padding:      " << getFramePadding(fp) << std::endl;
  unsigned int i = 0;
  PHDWORD *p = findFrameAlignBlock(fp);
  for (i = 0; i<  getFrameAlignLength(fp); i++ )
    {
      os << "  - Alignment word " << SETW(2) << i << ":  0x" ; 

      os.fill('0');
      os << SETW(8) << std::hex << *p++ << std::dec << std::endl;  
      os.fill (' ');
    }
  os << std::endl;
}

void
A_Event::dumpErrorBlock(PHDWORD *fp, OSTREAM &os)
{
  PHDWORD* ptr = findFrameErrorStart(fp);
  UINT len = getFrameErrorLength(fp);
  UINT nerr = calcNumErrorsV1(len);
  if ( nerr == 0 ) return;

  errorEntryV1* p = reinterpret_cast<errorEntryV1*>(ptr);
  for (UINT i=0; i<nerr; ++i)
    {
      errorEntryV1& e = *p;
      os << "ErrorEntry " << i << ": ";
      os << "severity: " << (int) e.severity << " "
	 << "deviceType: " << (int)e.deviceType << " "
	 << "deviceId: " << std::dec << e.deviceId << " "
	 << "errorCode: " << e.errorCode << " "
	 << "detectCode: " << e.detectCode << " "
	 << "addData: (" << std::hex << e.addData[0] << "," << std::hex << e.addData[1] << ")"
	 << std::dec
	 << std::endl;
      p++;
    }
  os << std::endl;
}

void
A_Event::dumpBlock(PHDWORD* p, UINT len, OSTREAM& os, const int how)
{
  if ( len == 0 ) { 
    os << "   (empty)\n" << std::endl;
    return; 
  }

  unsigned int j;
  switch (how)
    {
    case (EVT_HEXADECIMAL):
      j = 0;
      while (1)
	{
	  os << SETW(5) << j << " |  ";
	  for ( UINT l=0; l<4; l++ )
	    {
	      if ( j >= len ) break;
	      os << std::hex << SETW(8) << p[j++] << " " ;
	    }
	  if ( j >= len ) break;
	  os << std::dec << std::endl;
	}
      break;
		
    case (EVT_DECIMAL):
      j = 0;
      while (1)
	{
	  os << std::dec << SETW(5) << j << " |  ";
			 
	  for ( UINT l=0; l<6; l++ )
	    {
	      os << SETW(10) << p[j++] << " ";
	      if ( j >= len ) break;
	    }
	  if ( j >= len ) break;
	  os << std::endl;
	}
      break;
	 
    default: 
      break;
    }
  os << std::endl << std::endl;

}

unsigned int A_Event::getFrameEntry(const char *what, const int id, const int index) const
{
  
  PHDWORD *fp;
  PHDWORD *pp;
  int i = 0;
  int j = 0;

  if (id == 0) 
    {
      
      while ( ( fp = framelist[i++]) )
	{
	  return getFrameValue( what,fp,index);
	}
    }
  else // find the frame for a given packet
    {
      
      while ( (fp = framelist[i++]) )
	{
	  j++;
	  if ( ( pp = findFramePacketId (fp, id) ) !=  ptrFailure) 
	    {
	      return getFrameValue( what,fp,index);
	    }
	}
    }
  return 0;
}




unsigned int A_Event::getFrameValue(const char *what, PHDWORD *fp, const unsigned int index) const
{

  if ( strcmp(what,"FRAMELENGTH") == 0)     return getFrameLength(fp);
  else if  ( strcmp(what,"FRAMEMARK") ==0)  return  getFrameMark(fp);
  else if  ( strcmp(what,"FRAMEHDRVERSION") == 0)  return getFrameHdrVersion(fp);
  else if  ( strcmp(what,"FRAMEHDRLENGTH") == 0)  return  getFrameHdrLength(fp);
  else if  ( strcmp(what,"FRAMESTATUS") == 0) return  getFrameStatus(fp);
  else if  ( strcmp(what,"FRAMESOURCEID") ==0 )return  getFrameSourceId(fp);
  else if  ( strcmp(what,"FRAMEDATATYPE") == 0) return getFrameDataType(fp);
  else if  ( strcmp(what,"FRAMETYPE") == 0) return  getFrameType(fp);
  else if  ( strcmp(what,"FRAMEALIGNLENGTH") == 0) return getFrameAlignLength(fp);
  else if  ( strcmp(what,"FRAMEALIGNMENTWORD") == 0)
    {
      PHDWORD *p = findFrameAlignBlock(fp);
      if ( index >= getFrameAlignLength(fp) ) return 0;
      return p[index];
    }  

  return 0;
}


// getSubevent (int)
Packet* 
A_Event::getPacket (const int id)
{
  return getPacket(id,0);
}

#if !defined(SunOS) && !defined(OSF1)

int A_Event::createMap()
{
  int i = 0;
  PHDWORD *fp;
  PHDWORD *pp;

  if ( ! framelist ) 
    {
      errorcode = -3;
    }

  if (errorcode)  return 0;

  unsigned int pos_in_event;

  if (!hasMap)
    {
      
      while ( (fp = framelist[i++]) )
	{
	  pp = findFramePacketIndex (fp, 0);
	  

	  while ( pp  !=  ptrFailure) 
	    {

	      pos_in_event = ( int )(  pp - ( PHDWORD *) EventData );
	      
	      //	      std::cout << "pos in event = " << pos_in_event <<
	      //	"  packet length = " << *pp << " ev. length= " <<   getEvtLength() << std::endl;

	      if (pp && *pp > getEvtLength() - pos_in_event ) 
		{
		  std::cout << "Found wrong packet length " << *pp
			    << std::hex << "(0x" << *pp << ")" << std::dec 
			    << " packet Id: " <<  getPacketId(pp)  
			    << " Event: " << getEvtSequence() 
			    << " EvtLength: " << getEvtLength()
			    << " PosInEvent: " << pos_in_event
			    << std::endl;
		  errorcode =-2;
		  break;
		}
	      if ( pp != 0 && *pp == 0) 
		{
		  std::cout << "found 0-length packet" << std::endl;
		  errorcode =-1;
		  break;
		}


	      pmap[getPacketId(pp)] = pp;
	      // std::cout << "Packet id " << getPacketId(pp) << std::endl;

	      pp =  findNextFramePacket(fp, pp);

	    }
	}
      hasMap = 1;
    }
  return 0;
}


Packet* 
A_Event::getPacket (const int id, const int hitFormat)
{
  
  PHDWORD *pp;


  if (!hasMap) createMap();
  if ( errorcode) return 0;


  pp = pmap[id];
  if (!pp) return 0;

  return makePacket(pp,hitFormat);
}

#else

// the STL-Free solaris version

Packet* 
A_Event::getPacket (const int id, const int hitFormat)
{
  
  int i = 0;
  PHDWORD *fp;
  PHDWORD *pp;
  UINT ids = id;
  int wanted_hitformat;

  while ( fp = framelist[i++] )
    {
      if ( ( pp = findFramePacketId (fp, ids) ) !=  ptrFailure) 
	{
	  return makePacket(pp,hitFormat);
	}
    }
  return 0;
}
#endif

int 
A_Event::getPacketList( Packet* sl[], const int ne)
{
  int i = 0;
  PHDWORD *fp;
  PHDWORD *pp;

  if (!hasMap) createMap();
  if ( errorcode) return 0;

  int entries = 0;
  
  while ( (fp = framelist[i++]) )
    {

      pp = findFramePacketIndex (fp, 0);

      while ( pp  !=  ptrFailure) 
	{
	  if (getPacketStructure(pp) == Unstructured)
	    {
		  sl[entries++] = makePacket(pp,0);
		  //  sl[entries-1]->identify();
		  
	    }
	  
	  if (entries >= ne) return ne;
	  if ( (pp =  findNextFramePacket(fp, pp)) == ptrFailure)
            {
              break;
            }
	  if (*pp > getEvtLength()) 
	    {
	      std::cout << "Found wrong packet length " << *pp << std::endl;
	      break;
	    }
	  if ( pp != 0 && *pp == 0) 
	    {
	      std::cout << "found 0-length packet" << std::endl;
	      //  PHDWORD *x = pp - 10;
	      // std::cout << "--------------------------" << std::endl;
	      // for (i=0; i< 20; i++)
	      //	{
	      //	  std::cout << i << "  " << x << "  " << std::hex <<*x++ << std::dec << std::endl;
	      //	}
	      // std::cout << "--------------------------" << std::endl;
 

	      break;
	    }
	}
    }
  return entries;
}




Packet *A_Event::makePacket(PHDWORD *pp, const int hitFormat)
{

  int wanted_hitformat;

  if (getPacketStructure(pp) != Unstructured) return 0;
	  

  if (hitFormat)  wanted_hitformat = hitFormat;
  else wanted_hitformat  = getUnstructPacketHitFormat(pp);

  switch (wanted_hitformat)
    {

      // pbsc "32 channel format"

    case 50400:
    case IDHBD_FPGA:
    case IDHBD_FPGA0SUP:
    case IDFOCAL_FPGATEST:
      return new Packet_hbd_fpga(pp);
      break;

    case IDHBD_FPGASHORT:
    case IDHBD_FPGASHORT0SUP:
      return new Packet_hbd_fpgashort(pp);
      break;

    case IDCDEVPOLARIMETER:
      return new Packet_cdevpolarimeter(pp);
      break;

    case IDCDEVPOLARIMETERTARGET:
      return new Packet_cdevpoltarget(pp);
      break;

    case IDCDEVIR:
      return new Packet_cdevir(pp);
      break;

    case IDCDEVWCMHISTORY:
      return new Packet_cdevwcm(pp);
      break;

    case IDCDEVBPM:
      return new Packet_cdevbpm(pp);
      break;

    case IDCDEVDVM:
      return new Packet_cdevdvm(pp);
      break;

    case IDCDEVRING:
      return new Packet_cdevring(pp);
      break;

    case IDCDEVRINGPOL:
      return new Packet_cdevring(pp);
      break;

    case IDCDEVRINGFILL:
      return new Packet_cdevring(pp);
      break;

    case IDCDEVRINGNOPOL:
      return new Packet_cdevringnopol(pp);
      break;

    case IDCDEVBUCKETS:
      return new Packet_cdevbuckets(pp);
      break;

      //mlp 10/27/03 added this - the SIS is a straight array of numbers.
    case IDCDEVSIS:
      return new Packet_id4evt(pp);
      break;

    case IDCDEVMADCH:
      return new Packet_cdevmadch(pp);
      break;

    case ID4SCALER:
      return new Packet_id4scaler(pp);
      break;

    case IDGL1P:
      return new Packet_gl1p(pp);
      break;

    case IDGL1_EVCLOCK:
      return new Packet_gl1_evclocks(pp);
      break;

    case IDGL1PSUM:
    case IDGL1PSUMOBS:
      return new Packet_gl1psum(pp);
      break;

    case ID4EVT:
      return new Packet_id4evt(pp);
      break;

    case ID2EVT:
      return new Packet_id2evt(pp);
      break;
      
    case IDCSTR:
      return new Packet_idcstr(pp);
      break;
      
    case IDSTARSCALER:
      return new Packet_starscaler(pp);
      break;

    case IDCDEVDESCR:
      return new Packet_idcdevdescr(pp);
      break;


    default:
      switch (getUnstructPacketWordSize (pp))
	{
	case 1:
	  return new Packet_w1(pp);
	  break;
	case 2:
	  return new Packet_w2(pp);
	  break;
	case 4:
	  return new Packet_w4(pp);
	  break;
	default:
	  return new Packet_w4(pp); 
	}
    }
  
  return 0;
}



void  
A_Event::identify (OSTREAM &os) const
{ 
  os << " -- Event "   << EventData->evt_sequence;

  os << " Run: "  << SETW(5)  << EventData->run_number;

  os << " length: "     << SETW(5)  <<EventData->evt_length;

  os << " frames: " << SETW(3) << NumberFrames;

  os << " type: "      << SETW(2)  << EventData->evt_type ;

  os << " (" << get_evt_mnemonic(EventData->evt_type );


  if ( ( EventData->evt_type & CORRUPTEVENTMASK ) )
  {
    os << " *Corrupt* "; 
  }
  os << ") ";

  time_t x = getTime();
  os << x;
  
  os << std::endl;

};



int  
A_Event::existPacket (const int id)
{
#if !defined(SunOS) && !defined(OSF1)

  PHDWORD *pp;

  if (!hasMap) 
    {
      createMap();
    }

  pp = pmap[id];

  if (!pp) 
    {
      return 0;
    }
  return 1;
#else
  return 0;
#endif
}



// the Copy routine
int  
A_Event::Copy (int * array, const unsigned int length, int *nw, const char *what)
{
  if (length< getEvtLength() )
    {
      *nw = 0;
      return -1;
    }
  char *to = (char *) array;
  char *from;
  unsigned int l;
  if ( strcmp (what, "DATA") ==0  )
    {
      from= (char *)  &EventData->data[0];
      l = getEvtLength() - EVTHEADERLENGTH;
    }
  else
    {
      from= (char *) EventData;
      l = getEvtLength();
    }
  //  for (i=0; i<l ; i++) *to++ = *from++;
  //
  *nw = l;
  memcpy (to, from, l*4);
  return 0;
}

int A_Event::convert()
{
  if (is_data_type) return -1;

  PHDWORD *tmp;

  tmp = new PHDWORD[getEvtLength()];
  PHDWORD *from= (PHDWORD *)  EventData;
  PHDWORD *to=tmp;
  for (unsigned int k=0; k< getEvtLength(); k++) 
    { 
      *to++ = *from++;
    }

  delete [] framelist;
  EventData = (evtdata_ptr) tmp;
  updateFramelist();

  is_data_type = 1;
  pmap.clear();
  hasMap = 0;
  return 0;

}

int A_Event::is_pointer_type() const
{
  if (is_data_type) return 0;
  else return 1;
}


int A_Event::updateFramelist()
{
  // go through the data and see how may, 
  // if any, frames we have got.
  int max_index = EventData->evt_length - EVTHEADERLENGTH;
  int index = 0;
  int number_of_frames = 0;
  int flength;
  int cont;
  int i;
  // count the number of frames, and allocate a properly-sized vector
  // with pointers to it. 
  cont = 1;
  while (index < max_index && cont)
    {
      PHDWORD *f = &EventData->data[index];
      // so here we believe that we point to a frame start
      if (  validFrameHdr(f) )
	{
	  number_of_frames++;
	  flength = getFrameLength (f);
	  index += flength;
	}
      else 
	{
	  COUT << "invalid frame header, frame nr " << number_of_frames << " index = " << index <<  std::endl; //mlpd
	  COUT << " -- Event*"   << EventData->evt_sequence;

	  COUT << " Run: "  << SETW(5)  << EventData->run_number;

	  COUT << " length: "     << SETW(5)  <<EventData->evt_length;
	  
	  COUT << " type: "      << SETW(2)  << EventData->evt_type 
	       << " (" << get_evt_mnemonic(EventData->evt_type) << ") ";
	  COUT << std::endl;

	  for (i=0; i< 8; i++) 
	    COUT << i << "  " << std::hex << f[i] << std::dec << std::endl;

	  return -1;

	  // COUT << "I will continue anyway" << std::endl;
	  // cont = 0;
	}
    }

  // ok, so many frames. get the vector.
  framelist = new PHDWORD *[number_of_frames +1];
  
  // NumberFrames is a class data member.
  NumberFrames =  number_of_frames;

  // now we go through once more and remember the pointers to where the 
  // frames start
  index = 0;
  int ifn = 0;
  while (index < max_index && cont)
    {
      PHDWORD *f = &EventData->data[index];
      if (  validFrameHdr(f) )
	{
	  framelist[ifn++] = f;
	  flength = getFrameLength (f);
	  index += flength;
	}
      else 
	{
	  COUT << "invalid frame header, frame nr " << ifn << " index = " << index <<  std::endl; //mlpd
	  for (i=0; i< 8; i++) 
	    COUT << i << "  " << std::hex << f[i] << std::dec << std::endl;


	  cont = 0;
	}
    }

  //we terminate the list of frames with a 0.
  framelist[ifn] = 0;
  return 0;
}

int A_Event::getErrorCode()
{
#if !defined(SunOS) && !defined(OSF1)
  createMap();
#endif
  return errorcode;
}
