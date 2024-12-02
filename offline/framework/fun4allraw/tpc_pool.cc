#include "tpc_pool.h"

#include <stdint.h>
#include <string.h>

#include <Event/packet.h>

#define coutfl cout << __FILE__ << "  " << __LINE__ << " "
#define cerrfl cerr << __FILE__ << "  " << __LINE__ << " "

using namespace std;

tpc_pool::tpc_pool(const unsigned int depth, const unsigned int low_mark)
{
  if (depth < 4000000)
  {
    _required_depth = 4000000;
    cout << " Adjusted the depth to the minimum value " << endl;
  }
  else
  {
    _required_depth = depth;
  }

  if (low_mark > _required_depth)
  {
    _low_mark = _required_depth;
    cout << " Adjusted the low mark to not exceed minimum depth " << endl;
  }
  else
  {
    _low_mark = low_mark;
  }

  _is_decoded = 0;
  _last_requested_element = -1;  // impossible value to mark "unnused"
  _progress_index = 0;
}

tpc_pool::~tpc_pool()
{
  for (auto itr = waveforms.begin(); itr != waveforms.end(); ++itr)
  {
    delete (*itr);
  }
  waveforms.clear();

  for (auto itr = gtm_data.begin(); itr != gtm_data.end(); ++itr)
  {
    delete (*itr);
  }
  gtm_data.clear();
}

int tpc_pool::next()
{
  _is_decoded = 0;
  for (auto itr = waveforms.begin(); itr != waveforms.end(); ++itr)
  {
    delete (*itr);
  }
  waveforms.clear();

  for (auto itr = gtm_data.begin(); itr != gtm_data.end(); ++itr)
  {
    delete (*itr);
  }
  gtm_data.clear();

  for (int fee_id = 0; fee_id < MAX_FEECOUNT; fee_id++)
  {
    fee_data[fee_id].clear();
  }

  //  coutfl << "erasing " << _progress_index << " elements vector size is " << _thebuffer.size() << endl;
  _thebuffer.erase(_thebuffer.begin(), _thebuffer.begin() + _progress_index);
  // coutfl << " done, new size is " << _thebuffer.size() << endl;

  return 0;
}

bool tpc_pool::depth_ok() const
{
  return (min_depth() >= _required_depth);
}

int tpc_pool::addPacket(Packet *p)
{
  if (!p) return 0;

  p->identify();

  int l = p->getDataLength() + 16;
  int b[l];

  int l2;
  //  int l3 =  p->fillIntArray(b, l, &l2, "DATA");

  p->fillIntArray(b, l, &l2, "DATA");

  for (int i = 0; i < l2; i++)
  {
    unsigned short x = b[i] & 0xffff;
    _thebuffer.push_back(x);
    x = (b[i] >> 16) & 0xffff;
    _thebuffer.push_back(x);
  }

  return 0;
}

int tpc_pool::cacheIterator(const int n)
{
  if (n < 0) return 0;                         // not ok
  if (_last_requested_element == n) return 1;  // say "ok"

  unsigned int i = n;

  if (i >= waveforms.size())
  {
    _last_requested_element = -1;
    return 0;
  }

  auto it = std::next(waveforms.begin(), i);
  _cached_iter = it;
  _last_requested_element = i;
  return 1;
}

int tpc_pool::decode_gtm_data(unsigned short dat[16])
{
  unsigned char *gtm = (unsigned char *) dat;
  gtm_payload *payload = new gtm_payload;

  payload->pkt_type = gtm[0] | ((unsigned short) gtm[1] << 8);
  if (payload->pkt_type != GTM_LVL1_ACCEPT_MAGIC_KEY && payload->pkt_type != GTM_ENDAT_MAGIC_KEY)
  {
    delete payload;
    return -1;
  }

  payload->is_lvl1 = payload->pkt_type == GTM_LVL1_ACCEPT_MAGIC_KEY;
  payload->is_endat = payload->pkt_type == GTM_ENDAT_MAGIC_KEY;

  payload->bco = ((unsigned long long) gtm[2] << 0) | ((unsigned long long) gtm[3] << 8) | ((unsigned long long) gtm[4] << 16) | ((unsigned long long) gtm[5] << 24) | ((unsigned long long) gtm[6] << 32) | (((unsigned long long) gtm[7]) << 40);
  payload->lvl1_count = ((unsigned int) gtm[8] << 0) | ((unsigned int) gtm[9] << 8) | ((unsigned int) gtm[10] << 16) | ((unsigned int) gtm[11] << 24);
  payload->endat_count = ((unsigned int) gtm[12] << 0) | ((unsigned int) gtm[13] << 8) | ((unsigned int) gtm[14] << 16) | ((unsigned int) gtm[15] << 24);
  payload->last_bco = ((unsigned long long) gtm[16] << 0) | ((unsigned long long) gtm[17] << 8) | ((unsigned long long) gtm[18] << 16) | ((unsigned long long) gtm[19] << 24) | ((unsigned long long) gtm[20] << 32) | (((unsigned long long) gtm[21]) << 40);
  payload->modebits = gtm[22];

  this->gtm_data.push_back(payload);

  return 0;
}

int tpc_pool::tpc_decode()
{
  unsigned long long old_bco = 0;

  if (_is_decoded) return 0;
  _is_decoded = 1;

  unsigned int payload_length = _thebuffer.size();

  unsigned int index = 0;

  // demultiplexer
  while (index < payload_length)
  {
    // if ( index < 4) coutfl << index << " " << hex << _thebuffer[index] << dec << endl;

    // Length for the 256-bit wide Round Robin Multiplexer for the data stream
    const unsigned int datalength = 0xf;

    if ((_thebuffer[index] & 0xFF00) == FEE_MAGIC_KEY)
    {
      unsigned int fee_id = _thebuffer[index] & 0xff;
      // coutfl << " index = " << index << " fee_id = " << fee_id << " len = " << datalength << endl;
      ++index;
      if (fee_id < MAX_FEECOUNT)
      {
        for (unsigned int i = 0; i < datalength; i++)
        {
          // watch out for any data corruption
          if (index >= payload_length) break;

          fee_data[fee_id].push_back(_thebuffer[index++]);
        }
      }
    }

    else if ((_thebuffer[index] & 0xFF00) == GTM_MAGIC_KEY)
    {
      // if we hit the LVL1 BCO, we call it a day for this round
      if ((_thebuffer[index] & 0xFFFF) == GTM_LVL1_ACCEPT_MAGIC_KEY)
      {
        if (old_bco == 0)
        {
          old_bco = (_thebuffer[index + 3] & 0xff);
          old_bco <<= 8;
          old_bco |= (_thebuffer[index + 2] & 0xffff);
          old_bco <<= 16;
          old_bco |= (_thebuffer[index + 1] & 0xffff);
        }
        else
        {
          unsigned long long this_bco = 0;
          this_bco = (_thebuffer[index + 3] & 0xff);
          this_bco <<= 8;
          this_bco |= (_thebuffer[index + 2] & 0xffff);
          this_bco <<= 16;
          this_bco |= (_thebuffer[index + 1] & 0xffff);

          if (this_bco != old_bco)
          {
            // coutfl << "new GTM BCO, old " << hex << old_bco << " new: " << this_bco << dec << endl;
            //  coutfl << "erasing " << index << " elements vector size is " << _the_thebuffer.size() << endl;
            //  _thebuffer.erase(_thebuffer.begin(), _thebuffer.begin() + index );
            //  coutfl << " done, new size is " << _thebuffer.size() << endl;
            _progress_index = index;  // we are done until here
            break;
          }
        }
      }
      unsigned short buf[16];

      // memcpy?
      for (unsigned int i = 0; i < 16; i++)
      {
        buf[i] = _thebuffer[index++];
      }

      decode_gtm_data(buf);
    }
    else
    {
      // not FEE data, e.g. GTM data or other stream, to be decoded
      index += datalength + 1;
    }
  }

  _progress_index = index;  // we are done until here

  for (int ifee = 0; ifee < MAX_FEECOUNT; ifee++)
  {
    unsigned int pos;  // I'm too tired to play with iterators; pos is the position in a FEE vector

    // for ( fee_data_itr = fee_data[ifee].begin();  fee_data_itr != fee_data[ifee].end();  )
    for (pos = 0; pos < fee_data[ifee].size();)
    {
      // itr = fee_data_itr;
      int skip_amount = find_header(pos, fee_data[ifee]);
      if (skip_amount < 0) break;
      // for ( int i = 0; i < skip_amount; i++) ++fee_data_itr;  // skip that many words
      pos += skip_amount;

      // as we advance pos, let's remember where the start is
      unsigned int startpos = pos;

      // first the check if our vector cuts off before the fixed-length header, then we are already done
      if (startpos + HEADER_LENGTH >= fee_data[ifee].size() || startpos + fee_data[ifee][startpos] > fee_data[ifee].size())
      {
        pos = fee_data[ifee].size() + 1;  // make sure we really terminate the loop
      }
      else
      {
        // capture the header so we can easier get the bit shift stuff
        unsigned short header[HEADER_LENGTH];
        for (int i = 0; i < HEADER_LENGTH; i++) header[i] = (fee_data[ifee][pos++]);

        sampa_waveform *sw = new sampa_waveform;

        sw->fee = ifee;
        sw->pkt_length = header[0];
        sw->adc_length = header[5];
        sw->sampa_address = (header[1] >> 5) & 0xf;
        sw->sampa_channel = header[1] & 0x1f;
        sw->channel = header[1] & 0x1ff;
        sw->bx_timestamp = ((header[3] & 0x1ff) << 11) | ((header[2] & 0x3ff) << 1) | (header[1] >> 9);

        // now we add the actual waveform
        unsigned short data_size = header[5] - 1;

        //  coutfl << " Fee: " << ifee << " Sampa " << sw->sampa_address
        //  	 << " sampa channel: " << sw->sampa_channel
        //  	 << " channel: " << sw->channel
        //  	 << "  waveform length: " << data_size  << endl;

        for (int i = 0; i < data_size; i++)
        {
          sw->waveform.push_back(fee_data[ifee][pos++]);
        }

        // we calculate the checksum here because "pos" is at the right place
        unsigned short crc = crc16(ifee, startpos, header[0] - 1);
        // coutfl << "fee " << setw(3) << sw->fee
        // 	     << " sampla channel " << setw(3) <<  sw->channel
        // 	     << " crc and value " << hex << setw(5) << crc << " " << setw(5) << fee_data[ifee][pos] << dec;
        // if (  crc != fee_data[ifee][pos] ) cout << "  *********";
        // cout << endl;

        sw->checksum = crc;
        sw->valid = (crc == fee_data[ifee][pos]);

        waveforms.insert(sw);
      }

      // coutfl << "inserting at " << ifee*MAX_CHANNELS + sw->channel << " size is " << waveforms.size() << endl;
      // waveform_vector[ifee*MAX_CHANNELS + sw->channel].push_back(sw);
    }
  }

  return 0;
}

long long tpc_pool::lValue(const int n, const char *what)
{
  tpc_decode();

  const size_t i = n;

  if (strcmp(what, "N_TAGGER") == 0)  // the number of datasets
  {
    return gtm_data.size();
  }

  else if (strcmp(what, "TAGGER_TYPE") == 0)
  {
    if (i < gtm_data.size())
    {
      return gtm_data[i]->pkt_type;
    }
  }

  else if (strcmp(what, "IS_ENDAT") == 0)
  {
    if (i < gtm_data.size())
    {
      return gtm_data[i]->is_endat;
    }
  }

  else if (strcmp(what, "IS_LEVEL1_TRIGGER") == 0)
  {
    if (i < gtm_data.size())
    {
      return gtm_data[i]->is_lvl1;
    }
  }

  else if (strcmp(what, "BCO") == 0)
  {
    if (i < gtm_data.size())
    {
      return gtm_data[i]->bco;
    }
  }

  else if (strcmp(what, "LEVEL1_COUNT") == 0)
  {
    if (i < gtm_data.size())
    {
      return gtm_data[i]->lvl1_count;
    }
  }

  else if (strcmp(what, "ENDAT_COUNT") == 0)
  {
    if (i < gtm_data.size())
    {
      return gtm_data[i]->endat_count;
    }
  }

  else if (strcmp(what, "LAST_BCO") == 0)
  {
    if (i < gtm_data.size())
    {
      return gtm_data[i]->last_bco;
    }
  }

  else if (strcmp(what, "MODEBITS") == 0)
  {
    if (i < gtm_data.size())
    {
      return gtm_data[i]->modebits;
    }
  }

  return 0;
}

int tpc_pool::iValue(const int n, const int sample)
{
  if (n < 0) return 0;

  tpc_decode();

  if (cacheIterator(n))
  {
    unsigned int m = sample;
    if (m >= (*_cached_iter)->waveform.size()) return 0;
    return (*_cached_iter)->waveform[m];
  }
  return 0;
}

int tpc_pool::iValue(const int n, const char *what)
{
  tpc_decode();

  if (strcmp(what, "NR_WF") == 0)  // the number of datasets
  {
    return waveforms.size();
  }

  else if (strcmp(what, "MAX_FEECOUNT") == 0)
  {
    return MAX_FEECOUNT;
  }

  // see how many samples we have
  else if (strcmp(what, "SAMPLES") == 0)
  {
    if (cacheIterator(n))
    {
      return (int) (*_cached_iter)->waveform.size();
    }
    return 0;
  }

  else if (strcmp(what, "FEE") == 0)
  {
    if (cacheIterator(n))
    {
      return (int) (*_cached_iter)->fee;
    }
    return 0;
  }

  else if (strcmp(what, "SAMPAADDRESS") == 0)
  {
    if (cacheIterator(n))
    {
      return (int) (*_cached_iter)->sampa_address;
    }
    return 0;
  }

  else if (strcmp(what, "SAMPACHANNEL") == 0)
  {
    if (cacheIterator(n))
    {
      return (int) (*_cached_iter)->sampa_channel;
    }
    return 0;
  }

  else if (strcmp(what, "CHANNEL") == 0)
  {
    if (cacheIterator(n))
    {
      return (int) (*_cached_iter)->channel;
    }
    return 0;
  }

  else if (strcmp(what, "BCO") == 0)
  {
    if (cacheIterator(n))
    {
      return (int) (*_cached_iter)->bx_timestamp;
    }
    return 0;
  }

  else if (strcmp(what, "CHECKSUM") == 0)
  {
    if (cacheIterator(n))
    {
      return (int) (*_cached_iter)->checksum;
    }
    return 0;
  }

  else if (strcmp(what, "CHECKSUMERROR") == 0)
  {
    if (cacheIterator(n))
    {
      if ((*_cached_iter)->valid) return 0;
      return 1;
    }
    return 0;
  }

  return 0;
}

// int tpc_pool::find_header ( std::vector<unsigned short>::const_iterator &itr,  const std::vector<unsigned short> &orig)
int tpc_pool::find_header(const unsigned int yy, const std::vector<unsigned short> &orig)
{
  bool found = false;
  unsigned int pos = yy;
  std::vector<unsigned short> header_candidate;

  // we slide over the data and find the header, if any.
  // we calculate and return the amount of words we need to skip to find the vector.
  // if we find it right away, the amount returned is 0;
  // -1 for an error condition or end, like we hit the end without finding another header.
  // the complication here is that we look ahead a few words to recognize it.

  for (unsigned int i = 0; i < HEADER_LENGTH; i++)  // we start out with one "header candidate" = 7 entries
  {
    if (pos >= orig.size())  // if we reached the end, no more header here
    {
      return -1;
    }

    header_candidate.push_back(orig[pos]);
    pos++;
  }

  int skip_amount = 0;
  while (!found)
  {
    //      coutfl << " current pos is  " << pos  << "  vector length " << orig.size()  << endl;

    if (header_candidate[4] == MAGIC_KEY_0 && header_candidate[6] == MAGIC_KEY_1 && (header_candidate[0] - header_candidate[5] == HEADER_LENGTH))
    {
      // found it!
      found = true;
      //    coutfl << " found header skip value =  " << skip_amount << endl;
      break;
    }
    skip_amount++;
    if (pos >= orig.size())  // if we reached the end, no more header here
    {
      return -1;
    }

    //   coutfl << " next value " << pos << "  " << hex << orig[pos]  << dec << endl;
    header_candidate.erase(header_candidate.begin());  // delete the vector 1st element
    header_candidate.push_back(orig[pos]);             // push a new value, rinse and repeat
    pos++;
  }

  return skip_amount;
}

void tpc_pool::dump(std::ostream &os)
{
  tpc_decode();

  if (lValue(0, "N_TAGGER") == 0)
    os << "  No lvl1 and Endat taggers" << endl;
  else
  {
    os << "  TAGGER_TYPE    BCO          LEVEL1 CNT  ENDAT CNT     LAST_BCO     MODEBITS" << endl;

    for (int i = 0; i < lValue(0, "N_TAGGER"); ++i)  // go through the datasets
    {
      os << "  0x" << setw(4) << hex << lValue(i, "TAGGER_TYPE") << dec
         << " (" << (lValue(i, "IS_ENDAT") ? "ENDAT" : "") << (lValue(i, "IS_LEVEL1_TRIGGER") ? "LVL1 " : "")
         << ") "
         << setw(12) << lValue(i, "BCO") << " "
         << setw(10) << lValue(i, "LEVEL1_COUNT") << " "
         << setw(10) << lValue(i, "ENDAT_COUNT") << " "
         << setw(12) << lValue(i, "LAST_BCO")
         << "     0x" << std::setfill('0') << setw(2) << hex << lValue(i, "MODEBITS") << std::setfill(' ') << dec
         << endl;
    }

    os << endl;
  }
  waveform_set::const_iterator wf_iter;

  for (int i = 0; i < iValue(0, "NR_WF"); i++)  // go through the datasets
  {
    os << "  FEE   Channel   Sampachannel   Samples     BCO     CRC_ERR" << endl;

    os << setw(5) << iValue(i, "FEE") << " "
       << setw(9) << iValue(i, "CHANNEL") << " "
       << setw(9) << iValue(i, "SAMPACHANNEL") << " "
       << setw(12) << iValue(i, "SAMPLES") << " "
       << "     0x" << setw(5) << hex << iValue(i, "BCO") << dec
       //	 << " 0x" << setw(4) << hex << iValue(i, "CHECKSUM") << dec
       << setw(4) << iValue(i, "CHECKSUMERROR")
       << endl;

    for (int j = 0; j < iValue(i, "SAMPLES"); j += 10)
    {
      os << "                                                       ";
      for (int k = 0; k < 10; k++)
      {
        os << setw(4) << iValue(i, j + k) << " ";
      }
      os << endl;
    }
    os << endl;
  }
}

unsigned short tpc_pool::reverseBits(const unsigned short x) const
{
  unsigned short n = x;
  n = ((n >> 1) & 0x55555555) | ((n << 1) & 0xaaaaaaaa);
  n = ((n >> 2) & 0x33333333) | ((n << 2) & 0xcccccccc);
  n = ((n >> 4) & 0x0f0f0f0f) | ((n << 4) & 0xf0f0f0f0);
  n = ((n >> 8) & 0x00ff00ff) | ((n << 8) & 0xff00ff00);
  // n = (n >> 16) & 0x0000ffff | (n << 16) & 0xffff0000;
  return n;
}

unsigned short tpc_pool::crc16(const unsigned int fee, const unsigned int index, const int l) const
{
  unsigned short crc = 0xffff;

  for (int i = 0; i < l; i++)
  {
    unsigned short x = fee_data[fee][index + i];
    //      cout << "in crc " << hex << x << dec << endl;
    crc ^= reverseBits(x);
    for (unsigned short k = 0; k < 16; k++)
    {
      crc = crc & 1 ? (crc >> 1) ^ 0xa001 : crc >> 1;
    }
  }
  crc = reverseBits(crc);
  return crc;
}
