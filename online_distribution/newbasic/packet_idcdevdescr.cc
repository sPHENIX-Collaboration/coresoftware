#include "packet_idcdevdescr.h"
#include <stdio.h>
#include <vector>
#include <sstream>
#include <string.h>


Packet_idcdevdescr::Packet_idcdevdescr(PACKET_ptr data)
  : Packet_w1 (data)
{

  is_decoded = 0;
}
  
 

double Packet_idcdevdescr::dValue(const int ich, const char *what)
{

  if ( ! is_decoded )
    {

      decode();

    }
  
  static namedVector search;
  strcpy (  search.name, what );
  
  packetidcdevdescrSetiter it = nvset.find(&search);

  if ( it != nvset.end()) 
    {
      namedVector *x = *it;
      //      std::cout << what  << " x " << x->values.size() << std::endl;
      if ( ich < 0 || ich >= x->values.size()) return 0;
      return  (x->values)[ich];

    }
  
  return 0;


}



Packet_idcdevdescr::~Packet_idcdevdescr( )
{
  
   packetidcdevdescrSetiter it ;
   for (it =  nvset.begin(); it !=  nvset.end(); it++)
    {
      delete (*it);
      // nvset.erase(it); // commented out by K.Aoki
    }
   nvset.clear(); // added by K.Aoki

}

int *Packet_idcdevdescr::decode ( )
{
  is_decoded = 1;

  PHDWORD  *x = findPacketDataStart(packet);
 
  char* from = ( char  *) &x[1];

  int dlength = 4 * getDataLength() - 4;

  char *t = new char[dlength+1];
  strncpy ( t, from, dlength);
  
  char *token = strtok(t, ";");
  while (token)
    {
      struct namedVector *xx = new struct namedVector;
      double value;
      char * start = strstr( token, "=");
      if ( start)
	{
	  *start = '\0';
	  //	  std::cout << "---" << token << std::endl;
	  strcpy ( xx->name, token);

	  start++;
	  char *copy = NULL;
	  copy = new char[strlen(start)+1];
	  strcpy ( copy, start);
	
	  {
	    char *c;
	  
	    for ( c= copy; *c;)
	      {
		if ( *c == ',') *c = ' ';
		c++;
	      }
	  }

	  {
            std::stringstream iss(copy);
            while ( iss >> value )
	      {
		//		  std::cout  << " " << value << "  ";
		xx->values.push_back(value);
	      }
	    //	  os << copy << std::endl;
	    if ( nvset.find(xx) != nvset.end() ) // added by K.Aoki
	      {
		std::cout << "multiple entries\n"
			  << "   " << xx->name <<  std::endl;
		delete xx;
	      }
	    else
	      {
		nvset.insert ( xx);
	      }
	    //std::cout << "Token " << token << " vector length = " << xx->values.size() << std::endl;
	  }
	  delete[] copy;
	}else{
	  delete xx; /// added by K.Aoki
	}

      token = strtok(0, ";");
    }
  delete []t;  // added by K.Aoki

  return 0;
}

void Packet_idcdevdescr::dump (  OSTREAM& os)
{
 
  PHDWORD  *x = findPacketDataStart(packet);
 
  char* from = ( char  *) &x[1];

  int dlength = 4 * getDataLength() - 4;

  os << "Length = " << *x << std::endl;
 
  char *t = new char[dlength+1];
  strncpy ( t, from, dlength);
  
  char *token = strtok(t, ";");
  while (token)
    {
      os << token << std::endl;
      token = strtok(0, ";");
    }

  /*  
  os << std::endl;
  os << dValue(0, "bi8-rot3-2.3") << std::endl;
  os << dValue(0, "rbpm.bi8-bh1avgOrbTimeStamp") << std::endl;
  os << dValue(0, "pol_countUpRightY") << std::endl;
  os << dValue(0, "pol_downCountsY") << std::endl;
  os << dValue(0, "pol_maxTimeB") << std::endl;
  */

}
