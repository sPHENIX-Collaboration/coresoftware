#include "packet_hbd_fpgashort.h"


#define HBD_MAX_MODULES 4
//#define HBD_NSAMPLES 20

Packet_hbd_fpgashort::Packet_hbd_fpgashort(PACKET_ptr data)
  : Packet_w4 (data)
{
  nr_modules = 0;
  HBD_NSAMPLES = 24;
  hbd_parity=0;
}
  
int *Packet_hbd_fpgashort::decode ( int *nwout)
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

  memset( iarr, 0, 48*HBD_NSAMPLES*sizeof(int)*HBD_MAX_MODULES);
  memset( decoded_data2, 0, 4*HBD_MAX_MODULES);
  memset( decoded_data3, 0, 4*HBD_MAX_MODULES);
  memset( decoded_data4, 0, 4*HBD_MAX_MODULES);

  int pos = 0;
  int log_mod_nr;
  int l1_trig_nr;
  int beam_clock;
  int phys_mod_nr;

  while ( pos < dlength - 10)
    {
      //int i;
	  
      // the 0xF indicates the start of a new module

      if   ( (k[pos] & 0xF0000000) ==  0xF0000000   )  // and this is now pos
	{
	  nr_modules++;

	  //	  std::cout << pos << "  " << std::hex << k[pos] << std::dec <<   std::endl;

	  log_mod_nr  = k[pos] & 0xf;
	  l1_trig_nr  = k[pos+1] & 0xfff;
	  beam_clock  = k[pos+2] & 0xfff;
	  phys_mod_nr = k[pos+3] & 0x1f;

	  //  std::cout << "module nr: " <<  log_mod_nr 
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
	  
      // now we are supposed to get an 0xE which holds the 
      // channel address 
      while ( (k[pos] & 0xF0000000) ==  0xE0000000 )  // we have a first adc
	{

	  int f_channr     = (( k[pos] >>16 ) & 0x3f);
	  pos++;
	  int slot = 48*HBD_NSAMPLES*log_mod_nr + HBD_NSAMPLES*f_channr ; 

	  // here we mask out the two highest bits which need to be 0x4 
	  while  ( (k[pos] & 0xC0000000) ==  0x40000000 )  // the remaining ones
	    {
	      int adc0        =  ( k[pos] & 0xfff);  
	      iarr[slot++] = adc0;
	      int adc1        =  ( (k[pos]>>16) & 0xfff);  
	      iarr[slot++] = adc1;
	      pos++;

/*
	      	      std::cout <<  pos << "  " << std::hex << k[pos] << std::dec 
	      		<< "   " << log_mod_nr << "  " << f_channr  << "  "
	      		<< adc0  <<  "  " << adc1 << "  " << slot << std::endl;
*/


	    }


	}

    }
  if ( (k[pos] & 0xE0000000) ==  0xA0000000 )
    {
      // the parity data
      hbd_parity = (k[pos] & 0xfff);
      pos++;
    }

  
  *nwout = 48 * HBD_NSAMPLES * HBD_MAX_MODULES;
  return iarr;
}



int  Packet_hbd_fpgashort::iValue(const int ich, const int is)
{
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return 0;
    }

  if (ich < 0 || ich >= nr_modules *48) return 0;

  if ( is == 100) return ( iValue(ich,0) + iValue(ich,1));
  if ( is == 101) return ( iValue(ich,4));
  if ( is == 102) return ( iValue(ich,10) + iValue(ich,11));

  if (is < 0 || is >= HBD_NSAMPLES) return 0;

  return decoded_data1[ich*HBD_NSAMPLES + is]; 
}

// ------------------------------------------------------


int  Packet_hbd_fpgashort::iValue(const int ich, const char *what)
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

  else if (strcmp(what,"NRSAMPLES") == 0)  // how many samples does this format have?
    {			
      return HBD_NSAMPLES;
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

  else if (strcmp(what,"PARITY") == 0)  // user requested the module number
    {			

      if (decoded_data1 == NULL ) // no mistake, we decode this to get that info
	{
	  if ( (decoded_data1 = decode(&data1_length))==NULL)
	    return 0;
	}
      return hbd_parity;
	
    }

  return 0;
}




// ------------------------------------------------------

void Packet_hbd_fpgashort::dump ( OSTREAM &os)
{
  int i,j;

  this->identify(os);

  os << "  Number of Modules:  "  << SETW(8) << iValue(0,"NRMODULES") <<  "  Parity: "  << SETW(8) << std::hex << iValue(0,"PARITY") << std::dec  <<  std::endl; 
  os << "  Number of Samples:  "  << SETW(8) << iValue(0,"NRSAMPLES") << std::endl;


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
      os << "   ";
/*
      for ( j = 100; j < 103; j++)
	{
	  os <<  std::setw(5) << iValue(i,j);
	} 
*/
      
      os << std::endl;
    }
   dumpErrorBlock(os);
   dumpDebugBlock(os);

}
