#include "TPCMap.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

TPCMap::~TPCMap()
{
}


TPCMap::TPCMap( const std::string r1, const std::string r2, const std::string r3 )
{
  
  _broken = 0;
  int status ; 
  status = digest_map(r1, 0);
  status = digest_map(r2, 1);
  status = digest_map(r3, 2);

}

int  TPCMap::digest_map(const std::string s, const unsigned int section_offset)
{
  ifstream infile( s, ios::in);

  if (! infile.is_open() )
    {
      cerr << "Could not open file " << endl;
      _broken = 1;
      return -1;
    }
  
  string line;
  getline ( infile, line); // throwaway - we skip the first line
  //  cout << __FILE__<< " " << __LINE__ << ": " << line << endl;

  int abs_pad;
  int Radius;
  // int Pad;
  // int  U;
  // int  G;
  // string  Pin;
  // int PinColID;
  // int PinRowID;
  // string PadName;
  int  FEE;
  //  int FEE_Connector;
  int FEE_Chan;
  // double phi;
  // double x;
  // double y;
  // double PadX;
  // double PadY;
  double PadR;
  double PadPhi;

  while (  getline ( infile, line) )
    {
      //      cout << line<< endl;
      stringstream ss(line);
      string next;

  //  0  26    26
  //  1   0    0
  //  2   0    0
  //  3   1    1
  //  4   1    1
  //  5   0    C5
  //  6   2    2
  //  7   5    5
  //  8   0    ZZ.00.000
  //  9   5    5.0
  // 10   0    J2
  // 11 147    147
  // 12   0    0.005570199740407434
  // 13  69    69.99891405342764
  // 14   0    0.38991196551332985
  // 15  77    77.86043476908294
  // 16 305    305.14820499531316
  // 17 314    314.9248709046211
  // 18   0    0.24982557215053805




      
      int index = 0;
      while ( ss.good())
        {
          getline(ss,next, ',');
	  if ( index == 0)
	    {
	      abs_pad = stoul(next);
	    }
	  else if ( index == 1)
	    {
	      Radius = stoul(next);
	    }
	  else if ( index == 9)
	    {
	      FEE = stoul(next);
	    }
	  else if  ( index == 11)
	    {
	      FEE_Chan = stoul(next);
	    }
	  else if  ( index == 17)
	    {
	      PadR = stod(next);
	    }
	  else if  ( index == 18)
	    {
	      PadPhi = stod(next);
	    }
	  index++;
	}

      if ( section_offset == 1) FEE+=6;
      if ( section_offset == 2) FEE+=14;
	

      
      struct tpc_map x;
      x.padnr = abs_pad;
      x.layer = Radius + section_offset*16 + 7;
      x.FEE = FEE;
      x.FEEChannel = FEE_Chan;
      x.PadR = PadR;
      x.PadPhi = PadPhi;

      unsigned int key = 256*( FEE)  + FEE_Chan;
      tmap[key] = x;
      //      cout << setw(5) << key << setw(5) << FEE << setw(5) << FEE_Chan << " " << PadR << "  " << PadPhi << endl;
     
    }
  return 0;
}

unsigned int TPCMap::getLayer(const unsigned int FEE, const unsigned int FEEChannel, const unsigned int packetid) const
{
  if ( FEE>=26 || FEEChannel >255) return 0.;
  unsigned int key = 256*FEE + FEEChannel;

  std::map <unsigned int, struct tpc_map>::const_iterator itr = tmap.find(key);
  if ( itr == tmap.end()) return 0;
  return  itr->second.layer;

}

double TPCMap::getR(const unsigned int FEE, const unsigned int FEEChannel, const unsigned int packetid) const
{
  if ( FEE>=26 || FEEChannel >255) return 0.;
  unsigned int key = 256*FEE + FEEChannel;

  std::map <unsigned int, struct tpc_map>::const_iterator itr = tmap.find(key);
  if ( itr == tmap.end()) return -100;
  return  itr->second.PadR;

}

double  TPCMap::getPhi(const unsigned int FEE, const unsigned int FEEChannel, const unsigned int packetid) const
{
  if ( FEE > 25 || FEEChannel > 255) return 0.;
  unsigned int key = 256*FEE + FEEChannel;

  std::map <unsigned int, struct tpc_map>::const_iterator itr = tmap.find(key);
  if ( itr == tmap.end()) return -100;
  return  itr->second.PadPhi;

}

