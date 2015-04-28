#include "oncsSub_idcaenv1742.h"
#include <cstring>

oncsSub_idcaenv1742::oncsSub_idcaenv1742(subevtdata_ptr data)
  :oncsSubevent_w4 (data)
{
  samples = 0;
  dlength = 0;

  evnr = 0;
  freq=3;
  group_mask=0;
  int i;
  for ( i=0; i< 4; i++)
    {
      index_cell[i] = 0;
      tr_present[i] = 0;
    }

}
  
int *oncsSub_idcaenv1742::decode ( int *nwout)
{
  int *p;


  int i,channel;
  int *SubeventData = &SubeventHdr->data;


  // the first word must have  0xa in the MSB, and
  // has the total payload size in bits 0-27  
  if ( (( *SubeventData >> 28)  & 0xf ) != 0xa )
    {
      std::cout << "error in data structure" << std::endl;
      return 0;
    }
  dlength  = *SubeventData  & 0x0fffffff; 
 
  // word 1 has the group enable pattern
  group_mask = SubeventData[1] & 0xf;
  
  // word 2 has the event counter in bits 0-21
  evnr = SubeventData[2] & 0x3fffff;
  //  std::cout << "Evnt nr: " << evnr << std::endl;

  // and we ignore word 3 for now

  // before we go through the groups indexed by group_index,
  // we take a peek into the first group to figure out how many samples we have.
  // the sample count is the same for all groups, although each group encodes this  
  // value again. We use this value here up front to allocate the right amount of memory
  // for the decoded waveforms. 

  int size = SubeventData[4] & 0xfff;
  samples = size / 3;

  p = new int [ samples * 8 * 4];
  memset(p, 0, samples * 8 * 4 * sizeof(int));

  // the trigger waveform, if any
  int *p_tr = new int [ samples  * 4];
  memset(p_tr, 0, samples * 4 * sizeof(int));
  decoded_data2 = p_tr;
  data2_length = samples *4;


  // we also extract the sampling frequency here (encoded in each group but 
  // the same for all)
  freq = (SubeventData[4] >> 16) & 3;

  
  //now we go through the groups
  int group_offset = 4;
  int group_nr;

  for ( group_nr=0; group_nr < 4; group_nr++)
    {
      int pos;

      if  ( (group_mask >> group_nr) &1) // we have that group present
	{
	  int *groupdata = &(SubeventData[group_offset]); // first group
	  tr_present[group_nr] = (groupdata[0] >> 12) & 1;
	  index_cell[group_nr] = (groupdata[0] >> 20) & 0x3ff;


	  // std::cout << "group " << group_nr << " size:  " << size << std::endl;
	  // std::cout << "contains_tr: " << tr_present[group_nr] << std::endl;
	  // std::cout << "frequency  : " << freq << std::endl;
	  // std::cout << "index cell : " << index_cell[group_nr] << std::endl;

	  int s, ch;
	  pos = 1;

	  for ( s = 0; s < samples; s++ )
	    {
	      ch = 0;
	      p[group_nr*samples*8 + samples * ch++ + s] =( groupdata[pos] & 0xfff);
	      p[group_nr*samples*8 + samples * ch++ + s] =(groupdata[pos] >> 12) & 0xfff;
	      p[group_nr*samples*8 + samples * ch++ + s] = ((groupdata[pos] >> 24) & 0xff) + ( (groupdata[pos+1] & 0xf)<<8);
	      
	      p[group_nr*samples*8 + samples * ch++ + s] =(groupdata[pos+1] >> 4) & 0xfff;
	      p[group_nr*samples*8 + samples * ch++ + s] =(groupdata[pos+1] >> 16) & 0xfff;
	      p[group_nr*samples*8 + samples * ch++ + s] =((groupdata[pos+1] >> 28) & 0xf) + ( (groupdata[pos+2] & 0xff)<<4);
	      
	      p[group_nr*samples*8 + samples * ch++ + s] =(groupdata[pos+2] >> 8) & 0xfff;
	      p[group_nr*samples*8 + samples * ch++ + s] =(groupdata[pos+2] >> 20) & 0xfff;
	      pos +=3;

	    }
	  if ( tr_present[group_nr])
	    {
	      s = 0;
	      while (s < samples)
		{
		  p_tr[group_nr*samples + s++] = ( groupdata[pos] & 0xfff);
		  p_tr[group_nr*samples + s++] =(groupdata[pos] >> 12) & 0xfff;
		  p_tr[group_nr*samples + s++] = ((groupdata[pos] >> 24) & 0xff) + ( (groupdata[pos+1] & 0xf)<<8);
	      
		  p_tr[group_nr*samples + s++] =(groupdata[pos+1] >> 4) & 0xfff;
		  p_tr[group_nr*samples + s++] =(groupdata[pos+1] >> 16) & 0xfff;
		  p_tr[group_nr*samples + s++] =((groupdata[pos+1] >> 28) & 0xf) + ( (groupdata[pos+2] & 0xff)<<4);
	      
		  p_tr[group_nr*samples + s++] =(groupdata[pos+2] >> 8) & 0xfff;
		  p_tr[group_nr*samples + s++] =(groupdata[pos+2] >> 20) & 0xfff;
		  pos +=3;
		}
	    }
	}
      group_offset += pos + 1;
    }
  *nwout = samples*8 *4;

  return p;

}

int oncsSub_idcaenv1742::iValue(const int ch)
{

  if ( decoded_data1 == 0 ) decoded_data1 = decode(&data1_length);

  if ( ch < 0 || ch >= data1_length ) return 0;

  return decoded_data1[ch];

}

int oncsSub_idcaenv1742::iValue(const int sample, const int ch)
{

  if ( decoded_data1 == 0 ) decoded_data1 = decode(&data1_length);

  if ( ch < 0 || ch >= 32 ) return 0;
  if ( sample < 0 || sample >= samples ) return 0;

  return decoded_data1[ch*samples + sample];

}

int oncsSub_idcaenv1742::iValue(const int n,const char *what)
{

  if ( decoded_data1 == 0 ) decoded_data1 = decode(&data1_length);

  if ( strcmp(what,"SAMPLES") == 0 )
  {
    return samples;
  }

  if ( strcmp(what,"EVNR") == 0 )
  {
    return evnr;
  }

  if ( strcmp(what,"TR0") == 0 )
  {
    if ( n <0 || n >=samples) return 0;

    return decoded_data2[n];
  }

  if ( strcmp(what,"TR1") == 0 )
  {
    if ( n <0 || n >=samples) return 0;

    return decoded_data2[2*samples + n];
  }

  if ( strcmp(what,"INDEXCELL") == 0 )
  {
    if ( n <0 || n >=4) return 0;

    return index_cell[n];
  }

  if ( strcmp(what,"TRPRESENT") == 0 )
  {
    if ( n <0 || n >=4) return 0;

    return tr_present[n];
  }

  if ( strcmp(what,"FREQUENCY") == 0 )
  {
    return freq;
  }

  return 0;

}


void  oncsSub_idcaenv1742::dump ( OSTREAM& os )
{
  int i,j;
  //  int *SubeventData = &SubeventHdr->data;
  
  os << "Samples: " << iValue(0,"SAMPLES") << std::endl;
  os << "Evt Nr:  " << iValue(0,"EVNR") << std::endl;
  int f = iValue(0,"FREQUENCY") ;
  os << "Sample Frequency ";
  switch (f)
    {
    case 0:
      os << " 5 GS/s ";
      break;

    case 1:
      os << " 2.5 GS/s ";
      break;

    case 2:
      os << " 1 GS/s ";
      break;

    default:
      os << " Unknown ";
      break;
    }
  os << "("<< f << ")" << std::endl;

  os << "contains trigger sample: ";
  for ( i = 0; i < 4; i++)
    {
      os<< iValue(i, "TRPRESENT") << "  ";
    }
  os << std::endl;

  os << "index cell:              ";
  for ( i = 0; i < 4; i++)
    {
      os<<  iValue(i, "INDEXCELL") << "  ";
    }
  os << std::endl;
  os << std::endl;

  
  for ( i = 0; i < samples ; i++)
    {
      os  << std::setw(4) << i << " |  ";

      for ( j = 0; j < 32 ; j++)
	{
	  
	  os << std::setw(4) << iValue(i,j) << " ";
	}
      os << " tr: " << iValue(i, "TR0") << std::endl;
    }
  
  os << std::endl;
}
