#include "InttPool.h"
#include <string.h>

using namespace std;

#define coutfl std::cout << __FILE__<< "  " << __LINE__ << " "
#define cerrfl std::cerr << __FILE__<< "  " << __LINE__ << " "

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
 F_DATAWORD
};




InttPool::InttPool( const unsigned int depth)
{
  _depth = depth;
  
}

int InttPool::addPacket( Packet *p)
{
  int fee, i;

  for ( fee = 0; fee < 16; fee++)
    {
      for ( i = 0; i < p->iValue(fee, "FEE_LENGTH"); i++)
	{
	  fee_data[fee].push_back(p->iValue(fee,i, "") );
	}
    }
  return 0;
}

unsigned int InttPool::rawValue(const int fee, const int index)
{
  if ( fee < 0 || fee >= MAX_FEECOUNT) return 0;
  if ( index < 0 || (unsigned int) index >= fee_data[fee].size() ) return 0;
  return fee_data[fee][index];
}

int InttPool::iValue(const int fee, const char *what)
{
  intt_decode();
  int hit = fee;
  if ( strcmp(what,"NR_HITS") == 0)
    {
      return intt_hits.size();
    }

  if ( strcmp(what,"FEE_LENGTH") == 0)
    {
      if ( fee < 0 || fee >= MAX_FEECOUNT) return 0;
      return fee_data[fee].size();
    }

  else if ( strcmp(what,"ADC") == 0)
    {
      return iValue(hit,F_ADC);
    }
    
  else if ( strcmp(what,"AMPLITUDE") == 0)
    {
      return iValue(hit,F_AMPLITUDE);
    }
      
  if ( strcmp(what,"CHIP_ID") == 0)
    {
      return iValue(hit,F_CHIP_ID);
    }

  if ( strcmp(what,"CHANNEL_ID") == 0)
        {
	  return iValue(hit,F_CHANNEL_ID);
	}

  if ( strcmp(what,"FULL_FPHX") == 0)
    {
      return iValue(hit,F_FULL_FPHX);
    }

  if ( strcmp(what,"FEE") == 0)
    {
      return iValue(hit,F_FEE);
    }

  if ( strcmp(what,"FPHX_BCO") == 0)
    {
      return iValue(hit,F_FPHX_BCO);
    }

  if ( strcmp(what,"FULL_FPHX") == 0)
    {
      return iValue(hit,F_FULL_FPHX);
    }

  if ( strcmp(what,"FULL_ROC") == 0)
    {
      return iValue(hit,F_FULL_ROC);
    }

  if ( strcmp(what,"DATAWORD") == 0)
    {
      return iValue(hit,F_DATAWORD);
    }

  return 0;
}


long long InttPool::lValue(const int hit, const int field)
{
  intt_decode();
  if ( hit < 0 || hit >= (int) intt_hits.size()) return 0;

  switch (field)
    {
    case F_BCO:
      return intt_hits[hit]->bco;
      break;
    }
  
  return 0;
}

long long  InttPool::lValue(const int hit, const char *what)
{
  intt_decode();

  if ( strcmp(what,"BCO") == 0)
    {
      return lValue(hit,F_BCO);
    }
  return 0;
}



unsigned int InttPool::min_depth() const
{
  unsigned int d = 0;

  for ( int fee = 0; fee < MAX_FEECOUNT ; fee++)
    {
      if ( fee_data[fee].size() > d) d = fee_data[fee].size();
    }

  return d;
}

int InttPool::iValue(const int hit, const int field)
{
  intt_decode();
  if ( hit < 0 || hit >= (int) intt_hits.size()) return 0;

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
      
    case F_DATAWORD:
      return intt_hits[hit]->word;
      break;

    }

  return 0;
}

bool InttPool::depth_ok() const
{
  return ( min_depth() >= _depth);
}

int InttPool::next() 
{
  _is_decoded = 0;
  std::vector<intt_hit*>::const_iterator hit_itr;

  for ( hit_itr = intt_hits.begin(); hit_itr != intt_hits.end(); ++hit_itr)
    {
      delete (*hit_itr);
    }
  intt_hits.clear();
  
  return 0;
}


int InttPool::intt_decode ()
{

    if (_is_decoded) return 0;
    _is_decoded = 1;
  
  for ( int i = 0 ; i < MAX_FEECOUNT ; i++)
    {

      // cout << endl;
      // for (  unsigned int j = 0;  j <  fee_data[i].size(); j++)
      // 	{
      // 	  coutfl << " fee " << i << "  word  " << j << "  "  << hex << fee_data[i][j] << dec << endl;
      // 	}

      if ( fee_data[i].size() < _depth) continue;
      
      int header_found = 0;
      unsigned long long old_BCO = 0;
      unsigned int last_index = 0;
      
      for (  unsigned int j = 0;  j <  fee_data[i].size(); j++)
	{
	  //skip until we have found the first header
	  if (! header_found && (fee_data[i][j] & 0xff00ffff )!= 0xad00cade )
	    {
	      continue;
	    }
	  header_found = 1;
	  last_index = j;  // remember that we found a header here
	  
	  coutfl << "fee " << i << " found code 0x" << hex << fee_data[i][j] << dec << " last_index " << last_index << endl;
	  
	  unsigned long long BCO = 0;
	  unsigned long long l = 0;
	  
	  // 1st word  --- cade add9 87ea 0fe3 cade add9
	  l = fee_data[i][j];
	  // coutfl << "fee " << i << " BCO MSB " << hex << l << dec << endl;
	  BCO |= ( ((l >> 16 ) & 0xff) << 32);
	  
	  j++;
	  if ( j >= fee_data[i].size() ) continue;
	  l = fee_data[i][j];
	  
	  // coutfl << "fee " << i << " BCO mid " << hex << l << dec << endl;
	  BCO |= ( (l & 0xffff) << 16);
	  BCO |= ( (l >> 16) & 0xffff);
	  
	  coutfl << "BCO for fee " << setw(3) << i << " : " << hex << BCO << dec << endl;

	  if ( !old_BCO ) old_BCO = BCO;

	  if ( old_BCO && BCO != old_BCO) // ok, we have reached a new BCO here
	    {
	      coutfl << "found a new BCO for fee " << setw(3) << i << " : " << hex << old_BCO << " - " << BCO << dec << endl;
	      old_BCO = BCO;
	      break;
	    }
	  
	  // ok, now let's go until we hit the end, or hit the next header
	  while ( ++j < fee_data[i].size() )
	    {
	      if ( ( fee_data[i][j] & 0xff00ffff ) == 0xad00cade )
		{
		  header_found = 0;
		  j--;
		  break;
		}
	      
	      uint32_t x = fee_data[i][j];
	      
	      // 0x 0301 0063     
	      intt_hit * hit= new intt_hit;
	      hit->fee  = i;
	      hit->bco        = BCO;
	      hit->channel_id = (x >> 16) & 0x7f;  // 7bits
	      hit->chip_id    = (x >> 23) & 0x3f;  // 6
	      hit->adc        = (x >> 29) & 0x7;   // 3

	      hit->FPHX_BCO   = x  & 0x7f;
	      hit->full_FPHX  = (x >> 7) & 0x1;   // 1
	      hit->full_ROC   = (x >> 8) & 0x1;   // 1
	      hit->amplitude  = (x >> 9) & 0x3f;   // 1
	      
	      hit->word      = x;
	      //	      coutfl <<  " pushing back fee " << i <<  "  " << hex << x << dec  << endl;
	      
	      intt_hits.push_back(hit);
	      //coutfl << "list size: " << intt_hits.size() << endl;
	      //	      ++fee_data_itr;
	    }
	  
	}
      coutfl << "done with BCO at index " << last_index  << " size is " << fee_data[i].size() << endl;
      for ( unsigned int j = 0; j < last_index; j++)
	{
	  fee_data[i].erase(fee_data[i].begin());
	}
      // coutfl << " size is now " << fee_data[i].size() << endl;
	      
      
    }
  return 0;
}
