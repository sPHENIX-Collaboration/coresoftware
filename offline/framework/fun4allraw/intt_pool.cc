#include "intt_pool.h"
#include <string.h>

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
  F_DATAWORD
};

intt_pool::intt_pool(const unsigned int depth, const unsigned int low_mark)
{
  _required_depth = depth;
  _low_mark = low_mark;
  // last_index.fill(0);
  for (int fee = 0; fee < MAX_FEECOUNT; fee++)
  {
    last_index[fee] = 0;
  }
}

int intt_pool::addPacket(Packet *p)
{
  int fee, i;

  if (_myPacketid == -1)
  {
    _myPacketid = p->getIdentifier();
  }
  else
  {
    if (_myPacketid != p->getIdentifier())
    {
      cerrfl << " received packet " << p->getIdentifier() << " for pool for id " << _myPacketid << endl;
      return -1;
    }
  }

  // coutfl << " adding packet ";
  // p->identify();

  for (fee = 0; fee < MAX_FEECOUNT; fee++)
  {
    for (i = 0; i < p->iValue(fee, "FEE_LENGTH"); i++)
    {
      // coutfl << " pushing back for FEE " << setw(2)<<i << "  " << hex << p->iValue(fee,i, "") << dec <<endl;
      fee_data[fee].push_back(p->iValue(fee, i, ""));
    }
    // coutfl << "fee " << fee << " size now " << fee_data[fee].size() << endl;
  }

  return 0;
}

unsigned int intt_pool::rawValue(const int fee, const int index)
{
  if (fee < 0 || fee >= MAX_FEECOUNT) return 0;
  if (index < 0 || (unsigned int) index >= fee_data[fee].size()) return 0;
  return fee_data[fee][index];
}

int intt_pool::iValue(const int fee, const char *what)
{
  if (strcmp(what, "FEE_LENGTH") == 0)
  {
    if (fee < 0 || fee >= MAX_FEECOUNT) return 0;
    return fee_data[fee].size();
  }

  intt_decode();
  int hit = fee;

  if (strcmp(what, "NR_HITS") == 0)
  {
    return intt_hits.size();
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

  if (strcmp(what, "FULL_FPHX") == 0)
  {
    return iValue(hit, F_FULL_FPHX);
  }

  if (strcmp(what, "FULL_ROC") == 0)
  {
    return iValue(hit, F_FULL_ROC);
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
  if (hit < 0 || hit >= (int) intt_hits.size()) return 0;

  switch (field)
  {
  case F_BCO:
    return intt_hits[hit]->bco;
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
  return 0;
}

unsigned int intt_pool::min_depth() const
{
  unsigned int d = 0;

  for (int fee = 0; fee < MAX_FEECOUNT; fee++)
  {
    if (fee_data[fee].size() > d) d = fee_data[fee].size();
  }

  return d;
}

int intt_pool::iValue(const int hit, const int field)
{
  intt_decode();
  if (hit < 0 || hit >= (int) intt_hits.size()) return 0;

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

bool intt_pool::depth_ok() const
{
  if (verbosity > 5)
  {
    std::cout << "current Pool depth " << min_depth()
              << " required depth: " << _required_depth
              << std::endl;
  }
  return (min_depth() >= _required_depth);
}

int intt_pool::next()
{
  _is_decoded = 0;
  std::vector<intt_hit *>::const_iterator hit_itr;

  //  coutfl << "deleting " << intt_hits.size() << " hits"  << endl;

  for (hit_itr = intt_hits.begin(); hit_itr != intt_hits.end(); ++hit_itr)
  {
    //      coutfl << "deleting 0x" << hex << (*hit_itr)->bco << dec << endl;
    delete (*hit_itr);
  }
  intt_hits.clear();

  return 0;
}

int intt_pool::intt_decode()
{
  //  coutfl << " pool depth too small still: " << min_depth() << " required " << _depth << endl;
  if (!depth_ok())
  {
    return 0;
  }

  if (_is_decoded) return 0;
  _is_decoded = 1;

  for (int fee = 0; fee < MAX_FEECOUNT; fee++)
  {
    unsigned int j = 0;

    // for ( j = 0;  j <  fee_data[fee].size(); j++)
    // 	{
    // 	  coutfl << "fee " << fee << "  " << j << " found code 0x" << hex << fee_data[fee][j] << dec << endl;
    // 	}

    //      int go_on = 0;
    int header_found = 0;

    std::vector<unsigned int> hitlist;
    j = 0;

    unsigned int remaining = fee_data[fee].size() - _low_mark;
    if (fee_data[fee].size() < _low_mark) remaining = 0;

    while (j < (remaining))
    {
      // skip until we have found the first header
      if (!header_found && (fee_data[fee][j] & 0xff00ffff) != 0xad00cade)
      {
        // coutfl<<"skip until header fee "<<fee<<" j="<<j<<" "<<hex<<fee_data[fee][j]<<dec<<endl;

        j++;
        last_index[fee] = j;
        if (j > fee_data[fee].size()) coutfl << "Warning " << j << " " << fee_data[fee].size() << endl;
        continue;
      }
      header_found = 1;

      // here j points to a "cade" word

      // push back the cdae word, the BCO, and event counter
      if (fee_data[fee].size() - j >= 3)
      {
        for (int k = 0; k < 3; k++) hitlist.push_back(fee_data[fee][j++]);
        last_index[fee] = j;
      }
      else
      {
        coutfl << " Warning - size is " << fee_data[fee].size() << " probably cut off" << endl;
        break;
      }
      last_index[fee] = j;
      if (j > fee_data[fee].size()) coutfl << "Warning " << j << " " << fee_data[fee].size() << endl;
      // ok, now let's go until we hit the end, or hit the next header, or a footer

      while (j < fee_data[fee].size())  // note we don't stop at the "leftover" amount here
      {
        // we break here if find the next header or a footer
        if ((fee_data[fee][j] & 0xff00ffff) == 0xad00cade)
        {
          header_found = 0;
          j--;
          last_index[fee] = j;
          if (j > fee_data[fee].size()) coutfl << "Warning " << j << " " << fee_data[fee].size() << endl;
          // we have a full hitlist in the vector here
          coutfl << "calling decode with size " << hitlist.size() << endl;
          intt_decode_hitlist(hitlist, fee);
          hitlist.clear();
          break;
        }

        if (fee_data[fee][j] == 0xcafeff80)
        {
          // we have a full hitlist in the vector here
          // x		  coutfl << "calling decode with size " << hitlist.size() << endl;
          // coutfl << "calling decode for FEE " << fee << " with size " << hitlist.size() << endl;
          intt_decode_hitlist(hitlist, fee);
          hitlist.clear();
          j++;
          last_index[fee] = j;
          if (j > fee_data[fee].size()) coutfl << "Warning " << j << " " << fee_data[fee].size() << endl;

          break;
        }

        hitlist.push_back(fee_data[fee][j]);

        j++;
        last_index[fee] = j;
      }
      last_index[fee] = j;

      remaining = fee_data[fee].size() - _low_mark;
      if (fee_data[fee].size() < _low_mark) remaining = 0;
    }

    //--coutfl<<"fee "<<fee
    //--      <<", remaining "<<remaining
    //--      <<", datasize "<<fee_data[fee].size()
    //--      <<", lmark "<< _low_mark
    //--      <<", last "<<last_index[fee]<< " "
    //--      <<( fee_data[fee].size()>0&& last_index[fee]>=fee_data[fee].size() ? "WARNING last index exceeds fee_data_size" : "")
    //--      <<endl;

    //--for(unsigned int ii=last_index[fee]; ii<fee_data[fee].size(); ii++){
    //--  coutfl<<"      data:"<<ii<<"  "<<hex<<fee_data[fee][ii]<<dec<<endl;
    //--}

    // all data is not in this pool. need to wait next pool data.
    // remaining data should be decoded in the next pool (after all data comes)
    // to do this, last_index goes back to the last header
    if (hitlist.size() > 0 && fee_data[fee].size() > 0 && last_index[fee] >= fee_data[fee].size())
    {
      last_index[fee] -= hitlist.size();
      hitlist.clear();
      // coutfl<<" last_index changed : fee "<<fee<<" "<<last_index[fee]
      //       <<" "<<hex<<fee_data[fee][last_index[fee]]<<" "<<fee_data[fee][last_index[fee]+1]<<dec<<endl;
    }

    if (hitlist.size())
    {
      // coutfl << "calling decode for FEE " << fee << " with size " << hitlist.size() << endl;
      intt_decode_hitlist(hitlist, fee);
      hitlist.clear();
    }
  }

  for (int fee = 0; fee < MAX_FEECOUNT; fee++)
  {
    // coutfl << "FEE " << fee << " erasing  " << last_index[fee] << " words, size is  " << fee_data[fee].size() << endl;
    for (unsigned int j = 0; j < last_index[fee]; j++)
    {
      if (fee_data[fee].size()) fee_data[fee].erase(fee_data[fee].begin());
    }
    // coutfl << "FEE " << fee << " size is now " << fee_data[fee].size() << endl;
  }

  //--for ( int fee = 0 ; fee < MAX_FEECOUNT ; fee++)
  //--  {
  //--    coutfl<< "decode end : "<<fee<<" size "<<fee_data[fee].size()<<endl;
  //--  }

  return 0;
}

// coutfl << "next fee: " << fee << " j = " << j << endl;

// while ( j < fee_data[fee].size() )
// 	{

// 	  //skip until we have found the first header
// 	  if (! header_found &&  (fee_data[fee][j] & 0xff00ffff )!= 0xad00cade )
// 	    {
// 	      j++;
// 	      continue;
// 	    }
// 	  header_found = 1;
// 	  last_index = j;  // remember that we found a header here

// 	  //	  coutfl << "fee " << fee << " found code 0x" << hex << fee_data[fee][j] << dec << " last_index " << last_index << endl;

// 	  unsigned long long l = 0;

// 	  // 1st word  --- cade add9 87ea 0fe3 cade add9
// 	  l = fee_data[fee][j];
// 	  // coutfl << "fee " << i << " BCO MSB " << hex << l << dec << endl;
// 	  BCO |= ( ((l >> 16 ) & 0xff) << 32);
// 	  l = fee_data[fee][j+1];
// 	  BCO |= ( (l & 0xffff) << 16);
// 	  BCO |= ( (l >> 16) & 0xffff);
// 	  //coutfl << "BCO for fee " << setw(3) << fee << " : " << hex << BCO << dec << endl;

// 	  if ( !old_BCO ) old_BCO = BCO;

// 	  if ( old_BCO && BCO != old_BCO) // ok, we have reached a new BCO here
// 	    {
// 	      // coutfl << "found a new BCO for fee " << setw(3) << fee << " : " << hex << old_BCO << " - " << BCO << dec << endl;
// 	      old_BCO = BCO;
// 	      break;
// 	    }

// 	  // here j points to a "cade" word

// 	  // push back the cdae word, the BCO, and event counter
// 	  for ( int k = 0; k < 3; k++) hitlist.push_back(fee_data[fee][j++]);

// 	  int go_on = 1;
// 	  // ok, now let's go until we hit the end, or hit the next header, or a footer
// 	  while ( j < fee_data[fee].size() && go_on)
// 	    {

// 	      // we break here if find the next header or a footer
// 	      if ( ( fee_data[fee][j] & 0xff00ffff ) == 0xad00cade )
// 		{
// 		  header_found  = 0;
// 		  j--;
// 		  // we have a full hitlist in the vector here
// 		  coutfl << "calling intt_decode_hitlist with size " << hitlist.size() << endl;
// 		  intt_decode_hitlist (hitlist, fee);
// 		  hitlist.clear();
// 		  go_on = 0;
// 		}

// 	      if ( fee_data[fee][j] == 0xcafeff80 )
// 		{
// 		  // we have a full hitlist in the vector here
// 		  //coutfl << "calling intt_decode_hitlist with size " << hitlist.size() << endl;
// 		  intt_decode_hitlist (hitlist, fee);
// 		  hitlist.clear();
// 		  j++;
// 		  go_on = 0;
// 		}

// 	      hitlist.push_back(fee_data[fee][j]);
// 	      // coutfl << "pos " << j << " fee length " << fee_data[fee].size()
// 	      // 	     << " hit length now " <<  hitlist.size() << " 0x" << hex  << fee_data[fee][j] << dec << endl;
// 	      j++;
// 	    }

// 	}

// //      coutfl << " erasing the first " << last_index << " entries from fee " << fee << endl;
// for ( int j = 0; j < last_index; j++)
// 	{
// 	  fee_data[fee].erase(fee_data[fee].begin());
// 	}

// coutfl << "done with BCO 0x" << hex << BCO << dec << " at index " << last_index  << " fee size is " << fee_data[fee].size() << endl;
// coutfl << "calling intt_decode_hitlist with size " << hitlist.size() << endl;
// intt_decode_hitlist (hitlist, fee);
// hitlist.clear();

// coutfl << " size is now " << fee_data[i].size() << endl;

//     }
//   return 0;
// }

int intt_pool::intt_decode_hitlist(std::vector<unsigned int> &hitlist, const int fee)
{
  //  coutfl << " next hitlist, size " << hitlist.size() << " :" << endl;

  // for ( unsigned int i = 0; i < hitlist.size(); i++)
  //   {
  //     coutfl << i << " " << hex << hitlist[i] << dec << endl;
  //   }
  // cout << endl;

  if (hitlist.size() < 3)
  {
    coutfl << "hitlist too short " << endl;
    return 1;
  }

  unsigned long long BCO = 0;
  unsigned long long l = 0;

  l = hitlist[0];
  BCO |= (((l >> 16) & 0xff) << 32);
  l = hitlist[1];
  BCO |= ((l & 0xffff) << 16);
  BCO |= ((l >> 16) & 0xffff);
  unsigned int event_counter = hitlist[2];

  int count = 0;
  for (unsigned int i = 3; i < hitlist.size(); i++)
  {
    unsigned int x = hitlist[i];
    intt_hit *hit = new intt_hit;
    hit->event_counter = event_counter;
    hit->fee = fee;
    hit->bco = BCO;
    hit->channel_id = (x >> 16) & 0x7f;  // 7bits
    hit->chip_id = (x >> 23) & 0x3f;     // 6
    hit->adc = (x >> 29) & 0x7;          // 3

    hit->FPHX_BCO = x & 0x7f;
    hit->full_FPHX = (x >> 7) & 0x1;   // 1
    hit->full_ROC = (x >> 8) & 0x1;    // 1
    hit->amplitude = (x >> 9) & 0x3f;  // 1
    hit->word = x;
    if (verbosity > 1)
    {
      if (last_bco[fee] > BCO)
      {
        coutfl << "fee " << fee << " old bco : 0x" << std::hex
               << last_bco[fee] << ", current: 0x" << BCO
               << std::dec << std::endl;
      }
      cout << Name() << " pushing back hit for FEE " << fee << " with BCO 0x" << hex << BCO << dec
           << " chip " << hit->chip_id << " channel " << hit->channel_id << " hit length now " << intt_hits.size() << ", last bco: 0x" << hex << last_bco[fee] << dec << endl;
      last_bco[fee] = BCO;
    }
    intt_hits.push_back(hit);
    count++;
  }
  // coutfl << "pushed back " << count  << " hits for FEE " << fee << " with BCO 0x" << hex << BCO << dec
  // 	 << " size of hitlist now " << intt_hits.size() << endl;

  return 0;
}

void intt_pool::dump(OSTREAM &os)
{
  //  os << "number_of_hits: " << iValue(0, "NR_HITS") << endl;
  intt_decode();
  //  identify(os);

  os << "  Number of hits: " << iValue(0, "NR_HITS") << endl;

  //  std::vector::<intt_hit*>::const_iterator hit_itr;

  os << "   #    FEE    BCO      chip_BCO  chip_id channel_id    ADC  full_phx full_ROC Ampl." << endl;

  for (int i = 0; i < iValue(0, "NR_HITS"); i++)
  {
    os << setw(4) << i << " "
       << setw(5) << iValue(i, F_FEE) << " "
       << hex << setw(11) << lValue(i, F_BCO) << dec << "   "
       << hex << setw(2) << "0x" << iValue(i, F_FPHX_BCO) << dec << "   "
       << setw(5) << iValue(i, F_CHIP_ID) << " "
       << setw(9) << iValue(i, F_CHANNEL_ID) << "     "
       << setw(5) << iValue(i, F_ADC) << " "
       << setw(5) << iValue(i, F_FULL_FPHX) << " "
       << setw(9) << iValue(i, F_FULL_ROC)
       << setw(8) << iValue(i, F_AMPLITUDE)
       << "     "
       << "0x" << setw(8) << hex << setfill('0') << iValue(i, F_DATAWORD)
       << setfill(' ') << dec << endl;
  }
}
