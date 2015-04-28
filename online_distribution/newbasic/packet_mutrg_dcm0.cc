#include "packet_mutrg_dcm0.h"

Packet_mutrg_dcm0::Packet_mutrg_dcm0()
{
}

Packet_mutrg_dcm0::Packet_mutrg_dcm0(PACKET_ptr data)
  : Packet_w4 (data)
{

}

Packet_mutrg_dcm0::~Packet_mutrg_dcm0()
{
}

// ------------------------------------------------------


int  Packet_mutrg_dcm0::iValue(const int ich, const char *what)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now


  if ( strcmp(what,"EVTNR")==0)
    {
      if (decoded_data2 == NULL )
	{
	  if ( (decoded_data2 = decode_misc(&data2_length))==NULL)
	    return 0;
	}
      return decoded_data2[0];
    }

  else if ( strcmp(what,"FLAG")==0)
    {
      if (decoded_data2 == NULL )
	{
	  if ( (decoded_data2 = decode_misc(&data2_length))==NULL)
	    return 0;
	}
      
      return decoded_data2[1];
    }

  else if ( strcmp(what,"DETID")==0)
    {
      if (decoded_data2 == NULL )
	{
	  if ( (decoded_data2 = decode_misc(&data2_length))==NULL)
	    return 0;
	}
      
      return decoded_data2[2];
    }


  else if ( strcmp(what,"MODADDR")==0)
    {
      if (decoded_data2 == NULL )
	{
	  if ( (decoded_data2 = decode_misc(&data2_length))==NULL)
	    return 0;
	}
      
      return decoded_data2[3];
    }


  else if ( strcmp(what,"BCLCK")==0)
    {
      if (decoded_data2 == NULL )
	{
	  if ( (decoded_data2 = decode_misc(&data2_length))==NULL)
	    return 0;
	}
      
      return decoded_data2[4];
    }


  else if ( strcmp(what,"USERWORD")==0)
    {
      if (decoded_data2 == NULL )
	{
	  if ( (decoded_data2 = decode_misc(&data2_length))==NULL)
	    return 0;
	}
      
      return decoded_data2[5];
    }

  else if ( strcmp(what,"PARITY")==0)
    {
      if (decoded_data2 == NULL )
	{
	  if ( (decoded_data2 = decode_misc(&data2_length))==NULL)
	    return 0;
	}
      
      return decoded_data2[6];
    }



  return 0;
}



// ------------------------------------------------------

int Packet_mutrg_dcm0::iValue (const int strip)
{
    if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return 0;
    }


  // see if our array is long enough
  if (strip >= data1_length || strip < 0) return 0;

  return decoded_data1[strip];
}



int *Packet_mutrg_dcm0::decode (int *nwout)
{
  int *p,*k,*m;

    
  int dlength = getDataLength();

  m = (int *) findPacketDataStart(packet);
  if (m == 0) 
    {
      *nwout = 0;
      return 0;
    }
  k = &m[5]; // skip the header info

  if ( dlength < 235 ) // that's the length... 
    {
      *nwout = 0;   
      return 0;
    }

  
  p = new int [ 228];  

  int i;
  for ( i=0; i<228; i++)
    {
      p[i] = m[i] & 0xffff;

    }
  *nwout = 228;
  return p;
}

// ------------------------------------------------------

int *Packet_mutrg_dcm0::decode_misc (int *nwout)
{
  int *p,*k;


  k = (int *) findPacketDataStart(packet);
  if (k == 0) 
    {
      *nwout = 0;
      return 0;
    }
  
  p = new int[7];

  p[0] = k[0] & 0xffff;  // event number
  p[1] = k[1] & 0xffff;  // flag
  p[2] = k[2] & 0xffff;  // Detector ID
  p[3] = k[3] & 0xffff;  // module address
  p[4] = k[4] & 0xffff;  // beam clock counter
  p[5] = k[233] & 0xffff;  // userword
  p[6] = k[234] & 0xffff;  // pariyy

  *nwout = 5;
  return p;

}

// ------------------------------------------------------

void Packet_mutrg_dcm0::dump ( OSTREAM &os) 
{

  
  this->identify(os);
  
  os << " Event number:                    " << iValue(0,"EVTNR") << std::endl; 
  os << " Flag word                        " << std::hex << iValue(0,"FLAG") <<std::dec << std::endl;
  os << " Det id                           " << iValue(0,"DETID") << std::endl;

  os << " module address                   " << iValue(0,"MODADDR") << std::endl;
  os << " beam clock counter               " << iValue(0,"BCLCK") << std::endl;


  os << " Userword                         " << std::hex << iValue(0,"USERWORD") << std::dec <<std::endl;
  os << " Parity                           " << std::hex << iValue(0,"PARITY") << std::dec <<std::endl;

  int l ;
  int strip;
  int station;
  int octant;
  int bclk;
  int pos = 0;

  os << "   strip                                   value   pos"<< std::endl;

  while (pos<228)
    {
      for (bclk = 0; bclk < 3; bclk++)
	{
      
	  station = 1;
	  for ( strip=0; strip < 6; strip++) 
	    {
	      os << "BCLK " << bclk << " Octant A station " << station  
		 << " strips " << std::setw(3) << strip *16 +1 << " - " << std::setw(3) << strip*16 +16 << ": "   
		 << std::hex << std::setw(4) << iValue ( pos) << std::dec << "   (" << pos++ << ")" <<std::endl;
	    }
	  os << std::endl;


	  station = 2;
	  for ( strip=0; strip < 12; strip++) 
	    {
	      os << "BCLK " << bclk << " Octant A station " << station  
		 << " strips " << std::setw(3) << strip *16 +1 << " - " << std::setw(3) << strip*16 +16 << ": "   
		 << std::hex << std::setw(4) << iValue ( pos) << std::dec << "   (" << pos++ << ")" <<std::endl;
	    }
	  os << std::endl;

	  station = 3;
	  for ( strip=0; strip < 20; strip++) 
	    {
	      os << "BCLK " << bclk << " Octant A station " << station  
		 << " strips " << std::setw(3) << strip *16 +1 << " - " << std::setw(3) << strip*16 +16 << ": "   
		 << std::hex << std::setw(4) << iValue ( pos) << std::dec << "   (" << pos++ << ")" <<std::endl;
	    }
	  os << std::endl;

	  station = 1;
	  for ( strip=0; strip < 6; strip++) 
	    {
	      os << "BCLK " << bclk << " Octant B station " << station  
		 << " strips " << std::setw(3) << strip *16 +1 << " - " << std::setw(3) << strip*16 +16 << ": "   
		 << std::hex << std::setw(4) << iValue ( pos) << std::dec << "   (" << pos++ << ")" <<std::endl;
	    }
	  os << std::endl;

	  station = 2;
	  for ( strip=0; strip < 12; strip++) 
	    {
	      os << "BCLK " << bclk << " Octant B station " << station  
		 << " strips " << std::setw(3) << strip *16 +1 << " - " << std::setw(3) << strip*16 +16 << ": "   
		 << std::hex << std::setw(4) << iValue ( pos) << std::dec << "   (" << pos++ << ")" <<std::endl;
	    }
	  os << std::endl;

	  station = 3;
	  for ( strip=0; strip < 20; strip++) 
	    {
	      os << "BCLK " << bclk << " Octant B station " << station  
		 << " strips " << std::setw(3) << strip *16 +1 << " - " << std::setw(3) << strip*16 +16 << ": "   
		 << std::hex << std::setw(4) << iValue ( pos) << std::dec << "   (" << pos++ << ")" <<std::endl;
	    }
	  os << std::endl;

	}
      
    }
  

}








