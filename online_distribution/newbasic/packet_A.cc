#include "packet_A.h"
#include "packetHeaders.h"
#include  <string.h>
#include "buffer.h"

Packet_A::Packet_A()
{

  data1_length = 0;
  data2_length = 0;
  data3_length = 0;
  data4_length = 0;
  data5_length = 0;
  data6_length = 0;
  data7_length = 0;
  data8_length = 0;
  data9_length = 0;
  data10_length = 0;

  decoded_data1 = NULL;
  decoded_data2 = NULL;
  decoded_data3 = NULL;
  decoded_data4 = NULL;
  decoded_data5 = NULL;
  decoded_data6 = NULL;
  decoded_data7 = NULL;
  decoded_data8 = NULL;
  decoded_data9 = NULL;
  decoded_data10 = NULL;

  is_data_type = 0;  //assume "p" type first.
}

Packet_A::Packet_A(PACKET_ptr packet_ptr)
{
  
  packet = packet_ptr;
  //  packetHdr = ( PACKETHDR_ptr ) packet;

  data1_length = 0;
  data2_length = 0;
  data3_length = 0;
  data4_length = 0;
  data5_length = 0;
  data6_length = 0;
  data7_length = 0;
  data8_length = 0;
  data9_length = 0;
  data10_length = 0;

  decoded_data1 = NULL;
  decoded_data2 = NULL;
  decoded_data3 = NULL;
  decoded_data4 = NULL;
  decoded_data5 = NULL;
  decoded_data6 = NULL;
  decoded_data7 = NULL;
  decoded_data8 = NULL;
  decoded_data9 = NULL;
  decoded_data10 = NULL;


  is_data_type = 0;  //assume "p" type first.
}

// ----------------------------------------------

Packet_A::~Packet_A()
{

  if (decoded_data1 != NULL) 
    {
      delete [] decoded_data1;
      decoded_data1 = 0;
    }
  if (decoded_data2 != NULL) 
    {
      delete [] decoded_data2;
      decoded_data2 = 0;
    }
  if (decoded_data3 != NULL) 
    {
      delete [] decoded_data3;
      decoded_data3 = 0;
    }
  if (decoded_data4 != NULL) 
    {
      delete [] decoded_data4;
      decoded_data4 = 0;
    }
  if (decoded_data5 != NULL) 
    {
      delete [] decoded_data5;
      decoded_data5 = 0;
    }
  if (decoded_data6 != NULL) 
    {
      delete [] decoded_data6;
      decoded_data6 = 0;
    }
  if (decoded_data7 != NULL) 
    {
      delete [] decoded_data7;
      decoded_data7 = 0;
    }
  if (decoded_data8 != NULL) 
    {
      delete [] decoded_data8;
      decoded_data8 = 0;
    }
  if (decoded_data9 != NULL) 
    {
      delete [] decoded_data9;
      decoded_data9 = 0;
    }
  if (decoded_data10 != NULL) 
    {
      delete [] decoded_data10;
      decoded_data10 = 0;
    }

  if (is_data_type) delete [] packet;
}


int Packet_A::is_pointer_type() const
{
  if (is_data_type) return 0;
  else return 1;
}

// ----------------------------------------------

int Packet_A::convert()
{
  if (is_data_type) return -1;

  PHDWORD *tmp;

  tmp = new PHDWORD[getLength()];
  PHDWORD *from=  packet;
  PHDWORD *to=tmp;
  for (int k=0; k< getLength(); k++) 
    { 
      *to++ = *from++;
    }
  packet = tmp;
  is_data_type = 1;
  return 0;

}

// ------------------------------------------------------
int Packet_A::getLength() const
{
  return getPacketLength (packet) ;
}

// ------------------------------------------------------ 
int Packet_A::getErrorLength() const { return
getPacketErrorLength(packet); } 

//------------------------------------------------------ 
int Packet_A::getDebugLength() const { return
getPacketDebugLength(packet); } 

//------------------------------------------------------ 
int Packet_A::getIdentifier() const { return getPacketId (packet); } 

//------------------------------------------------------ 
int Packet_A::getPadding() const { 
return 
getPacketPadding(packet); 
}

// ------------------------------------------------------
int Packet_A::getStructure() const
{
  return getPacketStructure (packet);
}

// ------------------------------------------------------
int Packet_A::getHitFormat() const
{
  return getUnstructPacketHitFormat (packet);
}

// ------------------------------------------------------
//int Packet_A::getWordSize() const
//{
//  return ???  // ** action
//
//}

// ------------------------------------------------------
int Packet_A::getDataLength() const
{
  return getPacketDataLength (packet);
}

// ------------------------------------------------------

void  Packet_A::fullIdentify( OSTREAM& out ) const
{
  out << std::dec;
  out << "*** Packet id         " <<  getPacketId (packet) << std::endl;
  out << "    Hdr version       " <<  getPacketHdrVersion (packet) << std::endl;
  out << "    Hdr length        " <<  getPacketHdrLength (packet) << std::endl;
  out << "    Packet Length     " <<  getPacketLength (packet) << std::endl;
  out << "    Packet Status     " <<  getPacketStatus (packet) << std::endl;
  out << "    Debug Length      " <<  getPacketDebugLength (packet) << std::endl;
  out << "    Error Length      " <<  getPacketErrorLength (packet) << std::endl;
  out << "    Structure         " <<  getStructure() << " (";

  switch (getStructure()) 
    {
      
    case Unstructured:
      out << "Unformatted)  "; 
      break;
      
    case Packets:
      out << "Unformatted)  "; 
      out << "Packets    )  "; 
      break;
      
    case hitArray:
      out << "HitArray   )  "; 
      break;
      
    case hitList:
      out << "HitList    )  "; 
      break;
    }
  out << std::endl;
  
  out << "    Descriptor Words: " <<  getPacketDataDescrLength (packet) << std::endl;
  out << "    Endianism:        " <<  getPacketEndianism (packet);
  
  switch (getPacketEndianism (packet) )
    {
      
    case 1:
      out << " (Little Endian) " << std::endl;
      break;
      
    case 2:
      out << " (Big Endian) " << std::endl;
      break;
      
    default:
      out << " (Unknown) " << std::endl;
      break;
    }

  out << "    Padding:          " <<  getPacketPadding (packet) << std::endl;


}

void  Packet_A::identify( OSTREAM& out ) const
{
  // out << "identify packet type A" << std::endl;
  out << std::dec
      << "Packet " << SETW(6) <<  getIdentifier() 
      << " " << SETW(5) << getLength() 
      << " " << SETW(2) << getStructure() << " (";

  switch (getStructure()) 
    {

    case Unstructured:
      out << "Unformatted)  "; 
      break;
      
    case Packets:
      out << "Unformatted)  "; 
      out << "Packets    )  "; 
      break;
      
    case hitArray:
      out << "HitArray   )  "; 
      break;

    case hitList:
      out << "HitList    )  "; 
      break;
    }

  out << SETW(4) << getHitFormat() 
      << " (" << get_mnemonic( getStructure(), getHitFormat()) << ")";

  if (getErrorLength() ) out << " ** Error Block";

  out << std::endl;
}


// ------------------------------------------------------

int Packet_A::iValue(const int ich)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return 0;
    }

  // see if our array is long enough
  if (ich < 0 || ich >= data1_length) return 0;

  return decoded_data1[ich];
}

// ------------------------------------------------------

int   Packet_A::iValue(const int ich, const char *what)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return 0;
    }

  // see if our array is long enough
  if (ich >= data1_length) return 0;

  return decoded_data1[ich];
}

// ------------------------------------------------------

int   Packet_A::iValue(const int ich, const int iy)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return 0;
    }

  // see if our array is long enough
  if (ich > data1_length) return 0;

  return decoded_data1[ich];
}

// ------------------------------------------------------

float Packet_A::rValue(const int ich)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return 0;
    }

  // see if our array is long enough
  if (ich > data1_length) return 0;

  return float(decoded_data1[ich]);
}

// ------------------------------------------------------

float Packet_A::rValue(const int ich, const char *what)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return 0;
    }

  // see if our array is long enough
  if (ich > data1_length) return 0;

  return float(decoded_data1[ich]);
}

// ------------------------------------------------------

float Packet_A::rValue(const int ich, const int iy)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return 0;
    }

  // see if our array is long enough
  if (ich > data1_length) return 0;

  return float(decoded_data1[ich]);
}

// ----------------------------------------------

int   
Packet_A::getArraylength(const char *what)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return -1;
    }

  return data1_length;
}

// ----------------------------------------------

int
Packet_A::fillIntArray (int iarr[],
			const int nlen, int *nwout,
			const char *what)
{
  return standardIntArray(iarr, nlen, nwout, what);

}



// ----------------------------------------------

int
Packet_A::standardIntArray (int iarr[],
			const int nlen, int *nwout,
			const char *what)
{

  int *from;
  int howmuch;

  *nwout = 0;

  if (strcmp(what,"") == 0)
    {
      // now let's derefence the proxy array. If we didn't decode
      // the data until now, we do it now
      if (decoded_data1 == NULL )
	{
	  if ( (decoded_data1 = decode(&data1_length))==NULL)
	    {
	      *nwout=0;
	      return -1;
	    }
	}
      howmuch = data1_length;
      from = decoded_data1;
    }

  else if (strcmp(what,"RAW") == 0)
    {
      howmuch = getLength();
      from = (int *) packet;
    }

  else if (strcmp(what,"DATA") == 0)
    {
      howmuch = getDataLength();
      from = (int *) findPacketDataStart(packet);
    }

  else 
    {
      *nwout = 0;
      return 0;
    }

  // see if by any chance we got a negative length (happens)
  if (howmuch < 0)
    {
      *nwout = 0;
      return -3;
    }

  // see if our array is long enough
  if (nlen < howmuch) 
    {
      *nwout = 0;
      return -2;
    }

  if ( from == 0) 
    { 
      *nwout = 0;
      return -1;
    }
  
  // and copy the data to the output array
  memcpy(iarr, from, 4*howmuch);

  // tell how much we copied
  *nwout = howmuch;
  return 0;
}

// ------------------------------------------------------


int Packet_A::copyMe(int dest[],  const int maxlength) const
{
  
  if ( getLength() > maxlength )
    {
      return 0;
    }

  memcpy((void *) dest, (void *) packet, 4*getLength() );
  return getLength();
}



// ------------------------------------------------------



int
Packet_A::fillFloatArray (float rarr[],
			  const int nlen, int *nwout,
			  const char *what)
{

  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	{
	  *nwout=0;
	  return -1;
	}
    }

  // see if our array is long enough
  if (nlen < data1_length) return -2;

  // and copy the data to the output array
  int *from = decoded_data1;
  for (int i=0; i<data1_length; i++) *rarr++ = float (*from++);

  // tell how much we copied
  *nwout = data1_length;
  return 0;
}

// ------------------------------------------------------
int*   
Packet_A::getIntArray (int *nwout, const char *what)
{

  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	{
	  *nwout=0;
	  return NULL;
	}
    }

  int *temp=new int[data1_length];
  int is = fillIntArray(temp, data1_length, nwout);
  if (is)
    {
      *nwout=0;
      return NULL;
    }
  return temp;
}

// ------------------------------------------------------
float*   Packet_A::getFloatArray (int *nwout, const char *what)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	{
	  *nwout=0;
	  return NULL;
	}
    }

  float *temp=new float[data1_length];
  int is = fillFloatArray(temp, data1_length, nwout);
  if (is)
    {
      *nwout=0;
      return NULL;
    }
  return temp;
}

// ------------------------------------------------------
void  Packet_A::dumpErrorBlock ( OSTREAM& out ) 
{
  int el = getErrorLength();
  if (!el ) 
    {
      //     out << "--- No error block" << std::endl;
      return;
    }

  int j,l;

  int * data = (int *) findPacketErrorStart (packet);

  // here we make sure that the packet isn't obviously bogus,
  // although it's not a comprehensive check. 
  if ( el > getLength()  || data < (int *) packet ) return;


  j = 0;
  COUT << "--- Error block:" << std::endl;
  while (1)
    {
      out << SETW(5) << j << " |  ";
      for (l=0;l<4;l++)
	{
	  out << std::hex << SETW(8) << data[j++] << " " ;
	  if (j >= el ) break;
	}
      if (j >= el ) break;
      out << std::dec<< std::endl;
    }
  out << std::dec << std::endl;
}



// ------------------------------------------------------
void   Packet_A::dumpDebugBlock ( OSTREAM& out ) 
{
  int el = getDebugLength();
  if (!el ) 
    {
      //      out << "--- No Debug block" << std::endl;
      return;
    }

  int j,l;

  int * data = (int *) findPacketDebugStart (packet);

  if ( el > getLength()  || data < (int *) packet ) return;

  j = 0;
  COUT << "--- Debug block:" << std::endl;
  while (1)
    {
      out << SETW(5) << j << " |  ";
      for (l=0;l<4;l++)
	{
	  out << std::hex << SETW(8) << data[j++] << " " ;
	  if (j >= el ) break;
	}
      if (j >= el ) break;
      out << std::dec<< std::endl;
    }
  out <<  std::dec<<  std::endl;

  out << "--- DCM Checksum: " << getCheckSumStatus() << std::endl;

}

#ifdef LVL2_WINNT
  void fix_endianess ( LONGLONG *x)
#else
  void fix_endianess ( long long *x)
#endif
{
    unsigned int *i = (unsigned int *) x;
    unsigned int r = i[0];
    i[0] = i[1];
    i[1] = r;
}


void Packet_A::fix_endianess ( double *x)
{
    unsigned int *i = (unsigned int *) x;
    unsigned int r = i[0];
    i[0] = i[1];
    i[1] = r;
}

void Packet_A::fix_endianess ( char *str, const int length)
{
  int j;
  unsigned int *i = (unsigned int *) str;
  for (j = 0; j< length/4; j++) i[j] = buffer::i4swap(i[j]);

}



int Packet_A::getCheckSumStatus() const
{

  // if we don't have a debug block, wearedone (but don't
  // know what the answer is...

  if ( getDebugLength() == 0 ) return 0;
  int i;

  int *d = (int *) packet;

  int checksum_recalc = 0x0;
  for (i = 0;i < getLength() - 2;i++)
    {
      checksum_recalc ^= *d++;
    }
  d++;

  //  COUT << "Checksum : " << std::hex << checksum_recalc << "   " << *d << std::dec << ENDL;

  if (checksum_recalc != *d) return -1;
  return 1;
}


