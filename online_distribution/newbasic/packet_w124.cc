#include "packet_w124.h"

Packet_w1::Packet_w1()
{}


Packet_w1::Packet_w1(PACKET_ptr packet_ptr)
  : Packet_A(packet_ptr)
{
}

Packet_w2::Packet_w2()
{}

Packet_w2::Packet_w2(PACKET_ptr packet_ptr)
  : Packet_A(packet_ptr)
{
}

Packet_w4::Packet_w4()
{
}

Packet_w4::Packet_w4(PACKET_ptr packet_ptr)
  : Packet_A(packet_ptr)
{
}

// ---- and the dump routines, which just call the 
// ---- generic dump routines, but may be overloaded
// ---- by the subclasses.

void Packet_w1::dump(OSTREAM& out)  
{ 
  gdump(2,out);
  dumpErrorBlock(out);
  dumpDebugBlock(out);
}

void Packet_w2::dump(OSTREAM& out)  
{ 
  gdump(2,out);
  dumpErrorBlock(out);
  dumpDebugBlock(out);
}

void Packet_w4::dump(OSTREAM& out)  
{ 
  gdump(2,out);
  dumpErrorBlock(out);
  dumpDebugBlock(out);
}


void Packet_w4::gdump(const int i, OSTREAM& out) const
{
  int j,l;

  int * packetData = (int *) findPacketDataStart (packet);
  identify(out);

  if (packetData == 0) 
    {
      return;
    }

  switch (i)
    {
    case (EVT_HEXADECIMAL):
      j = 0;
      while (1)
	{
	  out << SETW(5) << j << " |  ";
	  for (l=0;l<4;l++)
	    {
	      out << std::hex << SETW(8) << packetData[j++] << " " ;
	      if (j>=getDataLength()) break;
	    }
	  if (j>=getDataLength() ) break;
	  out << std::dec<< std::endl;
	} 
      out << std::dec<< std::endl;
      break;
     
    case (EVT_DECIMAL):
      j = 0;
      while (1)
	{
	  out << std::dec << SETW(5) << j << " |  ";
			 
	  for (l=0;l<6;l++)
	    {
	      out << SETW(10) << packetData[j++] << " ";
	      if (j>=getDataLength() ) break;
	    }
	  if (j>=getDataLength() ) break;
	  out << std::endl;
	}
      out << std::dec<< std::endl;
      break;
	 
    default: 
      break;
    }
  out << std::dec << std::endl << std::endl;
 
}

// ---------------------------------------------------------------------

void Packet_w2::gdump(const int i, OSTREAM& out) const
{
  short * packetData = (short *) findPacketDataStart (packet);

  identify(out);
  if (packetData == 0) 
    {
      return;
    }

  //      << SETW(6)  << packetHdr->sub_length
  //     << SETW(5)  << packetHdr->sub_id
  //     << SETW(4)  << packetHdr->sub_type
  //     << SETW(4)  << packetHdr->sub_std::decoding
  //     << "  (" << get_mnemonic(packetHdr->sub_std::decoding) << ")";

  int j,l;

  switch (i)
    {
    case (EVT_HEXADECIMAL):
      j = 0;
      while (1)
	{
	  out << SETW(5) << j << " |  ";
	  for (l=0;l<8;l++)
	    {
	      out << std::hex << SETW(4) << packetData[j++] << " " ;
	      if (j>=2*getDataLength() ) break;
	    }
	  if (j>=2*getDataLength() ) break;
	  out << std::dec <<  std::endl;
	}
      out << std::dec <<  std::endl;
      break;

    case (EVT_DECIMAL):
      j = 0;
      while (1)
	{
	  out << SETW(5) << j << " |  ";
	  for (l=0;l<8;l++)
	    {
	      out << SETW(6) << packetData[j++] << " ";
	      if (j>=2*getDataLength()) break;
	    }
	  if (j>=2*getDataLength() ) break;
	  out << std::endl;
	}
      out << std::dec <<  std::endl;
      break;

    default: 
      break;
    }
  out << std::dec << std::endl << std::endl;

}
// ---------------------------------------------------------------------

void Packet_w1::gdump(const int i, OSTREAM& out) const
{

  char * packetData = (char *) findPacketDataStart (packet);

  identify(out);
  if (packetData == 0) 
    {
      return;
    }

  int j,l;
  char cstring[20];
  char *c;

  j = 0;
  switch (i)
    {
    case (EVT_HEXADECIMAL):
      while (1)
	{
	  c = cstring;
	  out << std::dec << SETW(5) << j << " |  ";
	  for (l=0;l<16;l++)
	    {
	      if (j < 4*(getDataLength()) ) 
		{
		  *c++ = packetData[j];
		  out << std::hex << SETW(2) << packetData[j++] << " ";
		}
	      else
		{
		  *c++ = 0x20;
		  out << std::dec << "   ";
		}
	    }
	  *c = 0;
	  out << "  | " << cstring;
	  if (j >= 4*getDataLength() ) break;
	  out << std::endl;
	}
      break;

    case (EVT_DECIMAL):
      while (1)
	{
	  c = cstring;
	  out << std::dec << SETW(5) << j << " |  ";
	  for (l=0;l<12;l++)
	    {
	      if (j < 4*getDataLength() ) 
		{
		  *c++ = packetData[j];
		  out << std::hex << SETW(4) << packetData[j++] << " ";
		}
	      else
		{
		  *c++ = 0x20;
		  out << std::dec << "   ";
		}
	    }
	  *c = 0;
	  out << "  | " << cstring;
	  if (j >= 4*getDataLength() ) break;
	  out << std::endl;
	}
      break;

    default: break;
    }
  out << std::dec << std::endl << std::endl;
}
