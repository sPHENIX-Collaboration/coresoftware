#include "intt_pool.h"

#include <Event/packet.h>

#include <algorithm>  // for max
#include <cstring>
#include <iomanip>  // for operator<<, setw, setfill

using namespace std;

#define coutfl std::cout << __FILE__ << "  " << __LINE__ << " "
#define cerrfl std::cerr << __FILE__ << "  " << __LINE__ << " "

enum ITEM
{
  F_BCO = 1,
  F_FEE,
  F_CHANNEL_ID,
  F_CHIP_ID,
  F_ADC,
  F_FPHX_BCO,
  F_FULL_FPHX,
  F_FULL_ROC,
  F_AMPLITUDE,
  F_EVENT_COUNTER,
  F_DATAWORD
};

intt_pool::intt_pool(const unsigned int depth, const unsigned int low_mark)
  : _required_depth(depth)
  , _low_mark(low_mark)
{
  // last_index.fill(0);
  for (int fee = 0; fee < MAX_FEECOUNT; fee++)
  {
    last_index[fee] = 0;
  }
  packetData = new unsigned int[depth + 2*16384];
  _allocated_size=depth + 2*16384;

}

intt_pool::~intt_pool()
{
  if (_allocated_size) { delete [] packetData;
}
}

int intt_pool::addPacket(Packet *p)
{
//  int fee, i;

  if (_myPacketid == -1)
  {
    _myPacketid = p->getIdentifier();
  }
  else
  {
    if (_myPacketid != p->getIdentifier())
    {
      cerrfl << " received packet " << p->getIdentifier() << " for pool for id " << _myPacketid << std::endl;
      return -1;
    }
  }

  //coutfl << " adding packet ";
  //p->identify();


  //currentpos = packetData.size();
  // coutfl << "at position " << currentpos << " adding " << p->getDataLength() << " words" << endl;

  //  packetData.resize(packetData.size() + p->getDataLength(), 0);

  

  int nw;
//  int status  = p->fillIntArray( (int *) &packetData[writeindex], p->getDataLength(), &nw,  "DATA");
  p->fillIntArray( (int *) &packetData[writeindex], p->getDataLength(), &nw,  "DATA");

  writeindex += nw;
  //coutfl << "status = " << status << " nw: " << nw << " writeindex = " << writeindex << endl;


  return 0;
}

unsigned int intt_pool::rawValue(const int fee, const int index)
{
  if (fee < 0 || fee >= MAX_FEECOUNT)
  {
    return 0;
  }
  if (index < 0 || (unsigned int) index >= fee_data[fee].size())
  {
    return 0;
  }
  return fee_data[fee][index];
}


int intt_pool::iValue(const int i, const int j, const char * what)
{

  // so here we have the index i of the BCO, and the position j of that FEE
  if ( strcmp(what,"FEELIST") == 0)
    {

      unsigned long long BCO = lValue(i,"BCOLIST");
      if ( BCO == 0) 
	{
	  return -1;
	}
      unsigned int uj = j;
      if ( j < 0 || uj >=FEEs_by_BCO[BCO].size() )
	{
	  return -1;
	}
      auto it = FEEs_by_BCO[BCO].cbegin();
      for (unsigned int k = 0; k< uj; k++) 
	{
	  ++it;
	}
      return *it;
    }

  int fee = i;
  int index=j;

  if ( fee < 0 || fee >= MAX_FEECOUNT) 
    {
      return 0;
    }

  if ( index < 0 || (unsigned int) index >= fee_data[fee].size() ) 
    {
      return 0;
    }

  intt_decode();
  return fee_data[fee][index];
}




int intt_pool::iValue(const int fee, const char *what)
{

  intt_decode();

  if (strcmp(what, "FEE_LENGTH") == 0)
  {
    if (fee < 0 || fee >= MAX_FEECOUNT)
    {
      return 0;
    }
    return fee_data[fee].size();
  }

  unsigned int ibco = fee; // it's not a fee, it's just an index for this one
  if ( strcmp(what,"NR_FEES") == 0)
    {
      unsigned long long BCO = lValue(ibco, "BCOLIST");
      if ( BCO == 0) 
	{
	  return 0;
	}
      return FEEs_by_BCO[BCO].size();
    }
  
  if (strcmp(what, "UNIQUE_FEES") == 0)
    {
      return FEE_List.size();
    }
  
  if (strcmp(what, "FEE_ID") == 0)
    {
      unsigned int ufee = fee;
      if (ufee > FEE_List.size())
	{
	  return -1;
	}
      auto it = FEE_List.begin();
      for (unsigned int k = 0; k < ufee; k++)
	{
	  ++it;
	}
      return *it;
    }
  
  if (strcmp(what, "FEE_BCOS") == 0)
    {
      if (fee < 0 || fee >= MAX_FEECOUNT)
	{
	  return 0;
	}
      return BCOs_by_FEE[fee].size();
    }



  int hit = fee;

  if (strcmp(what, "NR_HITS") == 0)
  {
    return intt_hits.size();
  }

  if ( strcmp(what,"NR_BCOS") == 0)
    {
      return BCO_List.size();
    }

  if (strcmp(what, "ADC") == 0)
  {
    return iValue(hit, F_ADC);
  }

  else if (strcmp(what, "AMPLITUDE") == 0)
  {
    return iValue(hit, F_AMPLITUDE);
  }

  if (strcmp(what, "CHIP_ID") == 0)
  {
    return iValue(hit, F_CHIP_ID);
  }

  if (strcmp(what, "CHANNEL_ID") == 0)
  {
    return iValue(hit, F_CHANNEL_ID);
  }

  if (strcmp(what, "FULL_FPHX") == 0)
  {
    return iValue(hit, F_FULL_FPHX);
  }

  if (strcmp(what, "FEE") == 0)
  {
    return iValue(hit, F_FEE);
  }

  if (strcmp(what, "FPHX_BCO") == 0)
  {
    return iValue(hit, F_FPHX_BCO);
  }

  //--if (strcmp(what, "FULL_FPHX") == 0)
  //--{
  //--  return iValue(hit, F_FULL_FPHX);
  //--}

  if (strcmp(what, "FULL_ROC") == 0)
  {
    return iValue(hit, F_FULL_ROC);
  }

  if (strcmp(what, "EVENT_COUNTER") == 0)
  {
    return iValue(hit, F_EVENT_COUNTER);
  }

  if (strcmp(what, "DATAWORD") == 0)
  {
    return iValue(hit, F_DATAWORD);
  }

  return 0;
}

long long intt_pool::lValue(const int hit, const int field)
{
  intt_decode();
  if (hit < 0 || hit >= (int) intt_hits.size())
  {
    return 0;
  }

  // NOLINTNEXTLINE(hicpp-multiway-paths-covered)
  switch (field)
  {
  case F_BCO:
    return intt_hits[hit]->bco;
    break;

  default:
    coutfl << "Unknown field " << field << std::endl;
    break;
  }

  return 0;
}

long long intt_pool::lValue(const int hit, const char *what)
{
  intt_decode();

  if (strcmp(what, "BCO") == 0)
  {
    return lValue(hit, F_BCO);
  }

  unsigned int i= hit; //  size() is unsigned
  if ( strcmp(what,"BCOLIST") == 0)
    {
      if ( hit < 0 || i >= BCO_List.size()) { return 0;
}
      auto it = BCO_List.cbegin();
      for (unsigned int j = 0; j< i; j++) { ++it;
}
      return *it;
    }

  return 0;
}


long long intt_pool::lValue(const int fee, const int i, const char *what)
{
  unsigned int ui= i; //  size() is unsigned
  if ( strcmp(what,"BCOVAL") == 0)
    {
      if (BCOs_by_FEE[fee].size() == 0)
	{
	  return -1;
	}
      if (ui > BCOs_by_FEE[fee].size())
	{
	  return -1;
	}
      auto it = BCOs_by_FEE[fee].cbegin();
      for (unsigned int j = 0; j < ui; j++)
	{
	  ++it;
	}
      return *it;
    }
  return 0;
}



unsigned int intt_pool::min_depth() const
{
  return writeindex - currentpos;
}

int intt_pool::iValue(const int hit, const int field)
{
  intt_decode();
  if (hit < 0 || hit >= (int) intt_hits.size())
  {
    return 0;
  }

  switch (field)
  {
  case F_FEE:
    return intt_hits[hit]->fee;
    break;

  case F_CHANNEL_ID:
    return intt_hits[hit]->channel_id;
    break;

  case F_CHIP_ID:
    return intt_hits[hit]->chip_id;
    break;

  case F_ADC:
    return intt_hits[hit]->adc;
    break;

  case F_FPHX_BCO:
    return intt_hits[hit]->FPHX_BCO;
    break;

  case F_FULL_FPHX:
    return intt_hits[hit]->full_FPHX;
    break;

  case F_FULL_ROC:
    return intt_hits[hit]->full_ROC;
    break;

  case F_AMPLITUDE:
    return intt_hits[hit]->amplitude;
    break;

  case F_EVENT_COUNTER:
    return intt_hits[hit]->event_counter;
    break;

  case F_DATAWORD:
    return intt_hits[hit]->word;
    break;

  default:
    coutfl << "Unknown field " << field << std::endl;
    break;
  }

  return 0;
}

bool intt_pool::depth_ok() const
{
  if (verbosity > 5)
    {
      std::cout << "current Pool depth " << writeindex - currentpos
		<< " required depth: " << _required_depth
		<< std::endl;
    }
  return ( min_depth()   >= _required_depth);
}

int intt_pool::next()
{
  _is_decoded = 0;


  // for ( int fee = 0 ; fee < MAX_FEECOUNT ; fee++)
  //   {
  //     fee_data[fee].clear();
  //   }




  std::vector<intt_hit *>::const_iterator hit_itr;

  //  coutfl << "deleting " << intt_hits.size() << " hits"  << std::endl;

  for (hit_itr = intt_hits.begin(); hit_itr != intt_hits.end(); ++hit_itr)
  {
    //      coutfl << "deleting 0x" << std::hex << (*hit_itr)->bco << std::dec << std::endl;
    delete (*hit_itr);
  }

  intt_hits.clear();
  BCO_List.clear();

  for(auto& [x, feelist] : FEEs_by_BCO)
    {
      feelist.clear();
    }
  for(auto& [x,bcolist] : BCOs_by_FEE)
    {
      bcolist.clear();
    }

  FEEs_by_BCO.clear();
  BCOs_by_FEE.clear();

  // now move the remaining stuff to the start of the array
  //coutfl << " current_pos: " << currentpos << " writeindex " << writeindex << endl;

  unsigned int l = (writeindex - currentpos) * sizeof(unsigned int);
  memcpy(packetData, &packetData[currentpos], l);
  writeindex -= currentpos;
  currentpos = 0;

  //coutfl << " current_pos: " << currentpos << " writeindex " << writeindex << endl;

  return 0;
}


int intt_pool::intt_decode ()
{

  if ( _is_decoded) { return 0;
}
  _is_decoded = 1;




  unsigned int payload_length = writeindex - _low_mark;
//  coutfl << " payload_length:  " << payload_length << " writeindex " << writeindex << endl;

  unsigned int index = currentpos;
  

  unsigned int *buffer =   &packetData[currentpos];

  while ( index < payload_length -1)
    {
      // find the "ba" index
      while ( (buffer[index] & 0xff00ffff ) !=  0xf000caf0 )
	{
	  coutfl << "skipping  at " << index << " values " << hex << buffer[index] << dec << endl;
	  index++;
	  if (index >= payload_length)
	    {
	      coutfl << " reached end at  " << index << " values " << hex << buffer[index] << dec << endl;
	      //_broken = 1;
	      return -1;
	    }
	}

      
      unsigned short fee = ( buffer[index] >> 20U ) & 0xfU;
      unsigned short len = ( (buffer[index] >> 16U) & 0xfU) >>1U;
      // coutfl << "found start at index " << index << " values " << hex << buffer[index+1] << dec << " fee: " << fee << " len: " << len << " BCO: 0x" << hex << calcBCO(&buffer[index+1]) << dec <<endl;
      index++;



      for ( int i = 0; i < len ; i++)
	{
	  //coutfl << "adding to ladder " << fee << "  " << hex << buffer[index] << dec << endl;
	  fee_data[fee].push_back(buffer[index++]);
	}

      if ( payload_length - index < 100  && buffer[index-1] == 0xcafeff80 )
	{
	  //coutfl << " found end at index " << index-1 << " 0x" << hex << buffer[index-1] << dec << endl;
	  break;
	}

    }




  //  coutfl << "fee_data: ";


  // for ( int fee = 0 ; fee < MAX_FEECOUNT ; fee++)
  //   {
  //     cout << fee_data[fee].size() << " ";
  //   }
  // cout << endl;

  //coutfl << "payload_length; " << payload_length << " index: " << index << endl;

  currentpos = index;


  // for ( int j = 0;  j <  fee_data[0].size(); j++)
  //   {
  //     coutfl << "feedata[0]  0x" << hex << fee_data[0][j] << dec << endl;
  //   }


  // for ( int fee = 0 ; fee < MAX_FEECOUNT ; fee++)
  //   {
  //     for ( int j = fee_data[fee].size() -10 ;  j <  fee_data[fee].size(); j++) 
  // 	{
  // 	  coutfl << "feedata[" << fee << "][" << j << "]  0x" << hex << fee_data[fee][j] << dec << endl;
  // 	}
  //     cout << endl;
  //   }


   for ( int fee = 0 ; fee < MAX_FEECOUNT ; fee++)
    {

      int end_here = fee_data[fee].size() -1;

      for (; end_here > 0; end_here--) 
  	{

	  if ( fee_data[fee][end_here] == 0xcafeff80 )
	    {
	      //coutfl << " found feedata[" << fee << "] end at index " << end_here <<  "  0x" << hex << fee_data[fee][end_here] << dec << endl;
	      break;
	    }
	}

      // for ( j = 0;  j <  fee_data[fee].size(); j++)
      //  	{
      //  	  coutfl << "fee " << fee << "  " << j << " found code 0x" << hex << fee_data[fee][j] << dec << endl;
      // 	}


      //      int go_on = 0;
      int header_found = 0;
      
      std::vector<unsigned int> hitlist;
      int j = 0;
      

      while ( j < end_here )
	{
	  
	  //skip until we have found the first header
	  if (! header_found &&  (fee_data[fee][j] & 0xff00ffff )!= 0xad00cade )
	    {
	      j++;
	      continue;
	    }
	  header_found = 1;

	  // here j points to a "cade" word

	  // push back the cdae word, the BCO, and event counter
	  if ( end_here -j >=3 )
	    {
	      for ( int k = 0; k < 3; k++) { hitlist.push_back(fee_data[fee][j++]);
}
	    }
	  else
	    {
	      coutfl << " Warning - index is " << j << " and size is " << end_here << endl;
	      j+= end_here -j +1;
	    }

	  // ok, now let's go until we hit the end, or hit the next header, or a footer
	  while ( j <= end_here )
	    {
	      
	      // we break here if find the next header or a footer
	      if ( ( fee_data[fee][j] & 0xff00ffff ) == 0xad00cade )
		{
		  header_found  = 0;
		  j--;
		  // we have a full hitlist in the vector here
		  coutfl << "calling decode for FEE " << fee << " with size " << hitlist.size() << endl;
		  intt_decode_hitlist (hitlist, fee);
		  hitlist.clear();
		  break;
		}
	      
	      
	      if ( fee_data[fee][j] == 0xcafeff80 )
		{
		  // we have a full hitlist in the vector here
		  //		  coutfl << "calling decode for FEE " << fee << " with size " << hitlist.size() << endl;
		  intt_decode_hitlist (hitlist, fee);
		  hitlist.clear();
		  j++;
		  break;
		}
	      
	      hitlist.push_back(fee_data[fee][j]);

	      j++;
	    }
	}

      hitlist.clear();
      //coutfl << " end of fee_data for FEE " << fee << " size: " << fee_data[fee].size() << " position : " << j << endl;

      fee_data[fee].erase(fee_data[fee].begin(), fee_data[fee].begin() + j);

      //coutfl << " Fee_data size now: " << fee_data[fee].size() << endl;



      
    }
   //coutfl << " data size: " << writeindex - currentpos << " current pos: " << currentpos << endl;
  return 0;
}

unsigned long long intt_pool::calcBCO(unsigned int *hitlist) const
{
  unsigned long long BCO = 0;
  unsigned long long l = 0;
  
  l = hitlist[0];
  BCO |= (((l >> 16U) & 0xffU) << 32U);
  l = hitlist[1];
  BCO |= ((l & 0xffffU) << 16U);
  BCO |= ((l >> 16U) & 0xffffU);
  return BCO;
}


int intt_pool::intt_decode_hitlist(std::vector<unsigned int> &hitlist, const int fee)
{
   //  coutfl << " next hitlist, size " << hitlist.size() << " :" << std::endl;

   // for ( unsigned int i = 0; i < hitlist.size(); i++)
   //   {
   //     coutfl << i << " " << std::hex << hitlist[i] << std::dec << std::endl;
   //   }
   // std::cout << std::endl;

  if (hitlist.size() < 3)
  {
    coutfl << "hitlist too short " << std::endl;
    return 1;
  }

  unsigned long long BCO = 0;
  unsigned long long l = 0;

  l = hitlist[0];
  BCO |= (((l >> 16U) & 0xffU) << 32U);
  l = hitlist[1];
  BCO |= ((l & 0xffffU) << 16U);
  BCO |= ((l >> 16U) & 0xffffU);
  // unsigned int event_counter = hitlist[2];
  unsigned int event_counter = 0;
  l = hitlist[2];
  event_counter |= ((l & 0xffffU) << 16U);
  event_counter |= ((l >> 16U) & 0xffffU);

  FEE_List.insert(fee);
  BCO_List.insert(BCO);
  FEEs_by_BCO[BCO].insert(fee);
  BCOs_by_FEE[fee].insert(BCO);

//  int count = 0;
  for (unsigned int i = 3; i < hitlist.size(); i++)
  {
    unsigned int x = hitlist[i];
    intt_hit *hit = new intt_hit;
    hit->event_counter = event_counter;
    hit->fee = fee;
    hit->bco = BCO;
    hit->channel_id = (x >> 16U) & 0x7fU;  // 7bits
    hit->chip_id = (x >> 23U) & 0x3fU;     // 6
    hit->adc = (x >> 29U) & 0x7U;          // 3

    hit->FPHX_BCO = x & 0x7fU;
    hit->full_FPHX = (x >> 7U) & 0x1U;   // 1
    hit->full_ROC = (x >> 8U) & 0x1U;    // 1
    hit->amplitude = (x >> 9U) & 0x3fU;  // 1
    hit->word = x;

    // if (verbosity > 1)
    //   {
    // 	if (BCO > last_bco[fee])
    // 	  {
    // 	    coutfl << "fee " << fee << " old bco : 0x" << std::hex
    // 		   << last_bco[fee] << ", current: 0x" << BCO
    // 		   << std::dec << std::endl;
    // 	  }
    // 	// std::cout << Name() << " pushing back hit for FEE " << fee << " with BCO 0x" << std::hex << BCO << std::dec
    // 	//      << " chip " << hit->chip_id << " channel " << hit->channel_id << " hit length now " << intt_hits.size() << ", last bco: 0x" << std::hex << last_bco[fee] << std::dec << std::endl;
    // 	last_bco[fee] = BCO;
    //   }


    //    coutfl << "count " << count << "  " << hit->bco << endl;  
    intt_hits.push_back(hit);
//    count++;
  }
  // coutfl << "pushed back " << count  << " hits for FEE " << fee << " with BCO 0x" << std::hex << BCO << dec
  // 	 << " size of hitlist now " << intt_hits.size() << std::endl;



  auto it = BCO_List.end();
  it--;

  //  coutfl << " last BCO value: 0x" << hex << *(it) << dec << endl;




  return 0;
}


void  intt_pool::dump ( OSTREAM& os )
{

  //  os << "number_of_hits: " << iValue(0, "NR_HITS") << endl;
  intt_decode();
  // identify(os);


  os << " Number of unique BCOs: " << iValue(0, "NR_BCOS") << endl;
  for ( int b = 0; b < iValue(0, "NR_BCOS"); b++)
    {
      os << " BCO " << setw(3) << b << ":  0x" << hex << lValue(b, "BCOLIST") << dec << "    number of FEEs for this BCO " << setw(3) << iValue(b,"NR_FEES") <<  endl;
      os << "           Number of unique FEEs: ";

      for ( int i = 0; i < iValue(b, "NR_FEES"); i++)
      	{
      	  os << " " << setw(3) << iValue(b, i, "FEELIST");
      	}
      os << endl;
    }

  os << " Number of hits: " << iValue(0, "NR_HITS") << endl;


  os << "   #    FEE    BCO      chip_BCO  chip_id channel_id    ADC  full_phx full_ROC Ampl." << endl;

  for ( int i = 0; i < iValue(0, "NR_HITS"); i++)
    {
      os << setw(4) << i << " "
	 << setw(5) <<             iValue(i, F_FEE)     << " "
	 <<  hex <<  setw(11) <<   lValue(i, F_BCO)  << dec << "   " 
	 <<  hex <<  setw(2) << "0x" <<  iValue(i,F_FPHX_BCO)  << dec  << "   " 
	 << setw(5) <<             iValue(i,F_CHIP_ID)    << " " 
	 << setw(9) <<             iValue(i,F_CHANNEL_ID) << "     "
	 << setw(5) <<             iValue(i,F_ADC)        << " " 
	 << setw(5) <<             iValue(i,F_FULL_FPHX) << " "
	 << setw(9) <<             iValue(i,F_FULL_ROC)
	 << setw(8) <<             iValue(i,F_AMPLITUDE) 
	 << "     " 
	 << "0x" << setw(8) <<  hex << setfill('0') << iValue(i,F_DATAWORD)
	 <<  setfill(' ') << dec << endl;
      
    }
  os << endl;
  
}

