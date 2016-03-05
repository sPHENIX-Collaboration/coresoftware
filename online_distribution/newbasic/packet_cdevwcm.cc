#include <packet_cdevwcm.h>


Packet_cdevwcm::Packet_cdevwcm(PACKET_ptr data)
  : Packet_w4 (data)
{
  ps = 0;
  numberofreadings= 0;
}
  
int *Packet_cdevwcm::decode ( int *nwout)
{

  int i;
  if (ps != 0) return 0; 

  int *k = (int *) findPacketDataStart(packet);
  if (k == 0) 
    {
      ps = 0;
      *nwout = 0;
      return 0;
    }
  numberofreadings = k[0];

  ps = ( struct cdevWCMHistory *) k;

  if (numberofreadings < 0 || numberofreadings >100) return 0;

  for ( i = 0; i< numberofreadings; i++) 
    {
       fix_endianess (&(ps->reading[i].beamcurrent));

    }

  *nwout = 0;

  return 0;
}


// ------------------------------------------------------


void Packet_cdevwcm::dump ( OSTREAM &os) 
{

  int i,j;

  decode (&i);

  this->identify(os);

  os << "Number of samples  " << numberofreadings  << std::endl;
  for ( i = 0; i < numberofreadings; i++)
    {
      os << "Beam current for Sample "  << std::setw(3) << i << " " 
	 << std::setw(12) << ps->reading[i].beamcurrent << 
	 " Time: " << std::setw(8)  << ps->reading[i].cdevCaptureTimeStamp << std::endl;
    }


  for ( i = 0; i < numberofreadings; i++)
    {
      os << " --- Bunch current values for Sample " << std::setw(3) << i << std::endl;

      for (j=0; j< 360; j++)
	{
	  if ( ps->reading[i].bunchcurrent[j])
	    {
	      os << std::setw(3) << i << "  " << std::setw(4) << j << "  " << 
		ps->reading[i].bunchcurrent[j] << std::endl;
	    }
	}
    }
  
  
  dumpErrorBlock(os);
  dumpDebugBlock(os);
}

int   Packet_cdevwcm::iValue(const int ich, const char *what)
{


  int i;
  decode (&i);

  if ( strcmp(what, "SAMPLES") == 0 ) return numberofreadings ;
  if ( strcmp(what, "TIMESTAMP") == 0 ) 
    {

      if (ich < 0 || ich >= numberofreadings) return 0;
      return ps->reading[ich].cdevCaptureTimeStamp;
  
    }

  std::cout << "packet_cdevwcm::iValue error unknown datum: " << what << std::endl;
  return 0;

}


float   Packet_cdevwcm::rValue(const int ich, const char *what)
{


  int i;
  decode (&i);

  if ( strcmp(what, "SAMPLES") == 0 ) return numberofreadings ;
  
  if ( strcmp(what, "BEAMCURRENT") == 0 ) 
    {

      if (ich < 0 || ich >= numberofreadings) return 0;
      return ps->reading[ich].beamcurrent;
  
    }

  std::cout << "packet_cdevwcm::rValue error unknown datum: " << what << std::endl;
  return 0;

}

float   Packet_cdevwcm::rValue(const int ich, const int y)
{


  int i;
  decode (&i);

  if (ich < 0 || ich >= numberofreadings) return 0;
  if (y < 0 || y >= 360) return 0;

  return ps->reading[ich].bunchcurrent[y];
  
}




