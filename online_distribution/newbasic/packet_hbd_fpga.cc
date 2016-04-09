#include "packet_hbd_fpga.h"


#define HBD_MAX_MODULES 4

Packet_hbd_fpga::Packet_hbd_fpga(PACKET_ptr data)
  : Packet_w4 (data)
{
  nr_modules = 0;
  if ( getHitFormat() == IDHBD_FPGA3SAMPLES)
    {
      HBD_NSAMPLES = 3;
    }
  else
    {
      HBD_NSAMPLES = 24;
    }
  std::cout << "Samples = " << HBD_NSAMPLES << std::endl;
}
  
int *Packet_hbd_fpga::decode ( int *nwout)
{

  int *k;

  int dlength = getDataLength();


  k = (int *) findPacketDataStart(packet);
  if (k == 0) 
    {
      *nwout = 0;
      return 0;
    }


  int *iarr = new int[ 48 * HBD_NSAMPLES * HBD_MAX_MODULES  ];


  decoded_data2 = new int[ HBD_MAX_MODULES ]; // trig nr
  decoded_data3 = new int[ HBD_MAX_MODULES ]; // beam clock
  decoded_data4 = new int[ HBD_MAX_MODULES ]; // module number in crate list

  memset( iarr, 0, 48*HBD_NSAMPLES*4*HBD_MAX_MODULES);
  memset( decoded_data2, 0, 4*HBD_MAX_MODULES);
  memset( decoded_data3, 0, 4*HBD_MAX_MODULES);
  memset( decoded_data4, 0, 4*HBD_MAX_MODULES);

  int pos = 0;

  while ( pos < dlength - 10)
    {
      int i;
      for ( i = 0; i< 10; i++)
	{
	  if   ( (k[pos] & 0xF0000000) !=  0x80000000   )
	    {
	      pos++;
	    }
	  else 
	    {
	      break;
	    }
	}
	  
      if   ( (k[pos] & 0xF00FF000) ==  0x800FF000   )  // begin of header is now pos
	{
	  nr_modules++;

	  //	  std::cout << pos << "  " << std::hex << k[pos] << std::dec <<   std::endl;

	  int log_mod_nr  = k[pos] & 0xf;
	  int l1_trig_nr  = k[pos+1] & 0xfff;
	  int beam_clock  = k[pos+2] & 0xfff;
	  int phys_mod_nr = k[pos+3] & 0x1f;

	  //std::cout << "module nr: " <<  log_mod_nr 
	  //	    << "  trig nr: "  <<  l1_trig_nr 
	  //	    << "  beam clock " <<  beam_clock 
	  //	    << "  phys. mod nr " <<  phys_mod_nr << std::endl;


	  if (  log_mod_nr >= HBD_MAX_MODULES)
	    {
	      std::cout << __FILE__ << " " <<  __LINE__ 
		    << " wrong logical module number " <<   log_mod_nr << std::endl;
	    }
	  else
	    {
	      decoded_data2[log_mod_nr] =  l1_trig_nr;
	      decoded_data3[log_mod_nr] =  beam_clock;
	      decoded_data4[log_mod_nr] =  phys_mod_nr;
	    }
	  pos += 4;
	}
      else
	{
	  delete [] iarr;
	  delete [] decoded_data2;
	  delete [] decoded_data3;
	  delete [] decoded_data4;
	  decoded_data2 = 0;
	  decoded_data3 = 0;
	  decoded_data4 = 0;
	  *nwout = 0;
	  return 0;
	}
	  

      while ( (k[pos] & 0xF0002000) ==  0x40002000 )  // we have a first adc
	{
	  int adc          =  ( k[pos] & 0xfff);
	  int clockphase   = (( k[pos] >>12 ) & 0x1);
	  int f_channr     = (( k[pos] >>16 ) & 0x3f);
	  int f_modnr      = (( k[pos] >>22 ) & 0x3);
	  int slot = 48*HBD_NSAMPLES*f_modnr + HBD_NSAMPLES*f_channr ; 
	  //std::cout <<  pos << "  " << std::hex << k[pos] << std::dec 
	  //	    << "   " << f_modnr << "  " << f_channr  << "  "
	  //	    << firstadc  << "  "<< clockphase  << "  "<< adc <<  "  " << slot << std::endl;
	  iarr[slot++] = adc;
	  pos++;
	  while ( (k[pos] & 0xF0002000) ==  0x40000000 )  // the remaining ones
	    {
	      int channr     = (( k[pos] >>16 ) & 0x3f);
	      int modnr      = (( k[pos] >>22 ) & 0x3);
	      if ( channr != f_channr ||  f_modnr != modnr ) 
		{
		  //	  std::cout <<  pos << " chan/mod nr changed unexpectedly " 
		  //	    << std::hex << k[pos] << std::dec << std::endl; 
		  break;
		}
	      adc        =  ( k[pos] & 0xfff);  
	      clockphase = (( k[pos] >>12 ) & 0x1);
	      // std::cout <<  pos << "  " << std::hex << k[pos] << std::dec 
	      //		<< "   " << f_modnr << "  " << f_channr  << "  "
	      //		<< firstadc  << "  "<< clockphase  << "  "<< adc <<  "  " << slot << std::endl;
	      iarr[slot++] = adc;
	      pos++;
	    }

	}
      if ( (k[pos] & 0xF0000000) ==  0x20000000 )
	{
	  // the parity data
	  pos++;
	}
    }
  *nwout = 48 * HBD_NSAMPLES * HBD_MAX_MODULES;
  return iarr;

}


int  Packet_hbd_fpga::iValue(const int ich, const int is)
{
  if (ich < 0 || ich >= nr_modules *HBD_NSAMPLES*48) return 0;
  if (is < 0 || is >= HBD_NSAMPLES) return 0;

  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return 0;
    }

  return decoded_data1[ich*HBD_NSAMPLES + is]; 
}

// ------------------------------------------------------


int  Packet_hbd_fpga::iValue(const int ich, const char *what)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now


  if (strcmp(what,"TRIGGER") == 0)  // user requested TRIGGER
    {			
      if (decoded_data1 == NULL ) // no mistake, we decode this to get that info
	{
	  if ( (decoded_data1 = decode(&data1_length))==NULL)
	    return 0;
	}
      if (ich < 0 || ich >= nr_modules) return 0;

      return decoded_data2[ich];  // ich really refers to sample index
    }

  else if (strcmp(what,"BCLK") == 0)  // user requested beam clock
    {			

      if (decoded_data1 == NULL ) // no mistake, we decode this to get that info
	{
	  if ( (decoded_data1 = decode(&data1_length))==NULL)
	    return 0;
	}
      if (ich < 0 || ich >= HBD_MAX_MODULES) return 0;

      return decoded_data3[ich];  // ich really refers to sample index
    }


  else if (strcmp(what,"MODULEID") == 0)  // user requested the module number
    {			

      if (decoded_data1 == NULL ) // no mistake, we decode this to get that info
	{
	  if ( (decoded_data1 = decode(&data1_length))==NULL)
	    return 0;
	}
      if (ich < 0 || ich >= HBD_MAX_MODULES) return 0;

      return decoded_data4[ich];  // ich really refers to sample index
    }

  else if (strcmp(what,"NRMODULES") == 0)  // user requested the module number
    {			

      if (decoded_data1 == NULL ) // no mistake, we decode this to get that info
	{
	  if ( (decoded_data1 = decode(&data1_length))==NULL)
	    return 0;
	}

      return nr_modules;
    }

  return 0;
}




// ------------------------------------------------------

void Packet_hbd_fpga::dump ( OSTREAM &os)
{
  int i,j;

  this->identify(os);

  os << "  Number of Modules:  "  << SETW(8) << iValue(0,"NRMODULES") <<  std::endl; 

  for (i = 0; i < iValue(0,"NRMODULES"); i++)
    {
      os << "  Module # " << std::setw(2) << i 
	 <<  " Trigger: "  << SETW(8) << iValue(i,"TRIGGER") 
	 <<  " Beam Clock: "  << SETW(8) << iValue(i,"BCLK") 
	 <<  " Module Id: "  << SETW(8) << iValue(i,"MODULEID") 
	 <<  std::endl; 
    }

  for (i = 0; i < iValue(0,"NRMODULES") * 48 ; i++)
    {
      os << std::setw(5) << i << " | " ;
      for ( j = 0; j < HBD_NSAMPLES; j++)
	{
	  os <<  std::setw(5) << iValue(i,j);
	} 

      os << std::endl;
    }
   dumpErrorBlock(os);
   dumpDebugBlock(os);

}
