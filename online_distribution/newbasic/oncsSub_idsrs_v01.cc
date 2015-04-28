#include "oncsSub_idsrs_v01.h"
#include <cstring>

#include <arpa/inet.h>

using namespace std;

oncsSub_idsrs_v01::oncsSub_idsrs_v01(subevtdata_ptr data)
  :oncsSubevent_w4 (data)
{
  //  samples = 0;
  // nsamples = 0;
  //  hybrid = 0;
  nhybrids = 0;
  
}
  
oncsSub_idsrs_v01::~oncsSub_idsrs_v01()
{
  //  if ( samples )
  //  {
  //    delete [] samples;
  //    samples = 0;
  //    nsamples = 0;
  //  }

  std::vector<hybriddata*>::iterator it;
  std::vector<report *>::iterator reportit;


  for ( it = hybridlist.begin(); it != hybridlist.end(); ++it)
    {
      for ( reportit = (*it)->rowdata.begin(); reportit != (*it)->rowdata.end(); ++reportit)
	{
	  delete (*reportit);
	}
      delete (*it);
    }

}


int *oncsSub_idsrs_v01::decode ( int *nwout)
{
  int i;

  unsigned int *d = (unsigned int *) &SubeventHdr->data;  // here begins the payload

  int dlength = getLength()-4 - getPadding();

  for ( i = 0; i< dlength; i++)
    {
      hybriddata *hd = new hybriddata;
      if ( d[i] == 0xfafafafa ) break;  // end of payload
      unsigned int dh = ntohl(d[i]);
      int framecounter = (dh & 0xff);
      //cout << "framecounter: " << framecounter << endl;
      
      hd->framecounter = framecounter;
      
      dh = ntohl(d[i+1]);
      int header = dh;

      unsigned char c = (header >>24) & 0xff;
      hd->desc[0] = c;
      c = (header >>16) & 0xff;
      hd->desc[1] = c;
      c = (header >>8) & 0xff;
      hd->desc[2] = c;
      hd->desc[3] = 0;

      unsigned int channel = header & 0xff;
      //cout << " channel nr: " << channel << endl;
      hd->hdmi = channel;

      unsigned int w = ntohl(d[i+2]);
      int words = w & 0xffff;
      //cout << " words: " << words << endl;

      int j = 0;  // here begin the samples;
      int k = 0;  // here begin the samples;

      while ( d[j + i +3] != 0xf000f000 && j <  words )
	{
	  hd->adc.push_back(d[j+i+3] & 0x0fff);
	  hd->adc.push_back( (d[j+i+3] >> 16) & 0x0fff);
	  j++;
	}


      //cout << " EOP marker " << hex <<  d[j + i +3] << dec << endl;
      hd->words = hd->adc.size();
      //cout << " words: " << hd->words << "  " << words << endl;

      analyze ( hd);

      hybridlist.push_back(hd);
      nhybrids++;
      i+= 3+j;
    }

  *nwout = 0;
  return 0;
}


int oncsSub_idsrs_v01::analyze (hybriddata * hd)
{
  int i;

  int firstfound = 0;
  int onecount = 0;

  int current_report = 0;

  for ( i = 0; i < hd->words-( 128+8); i++)
    {

      if ( firstfound) 
	{
	  add_report ( hd, i, current_report++);
	  i += 11+128;
	}
      else
	{
	  if ( hd->adc[i] < 1220) onecount++;
	  if ( onecount == 3) 
	    {
	      // cout << " sequence start at " << i -2 << endl;
	      add_report ( hd, i-2, current_report++);
	      onecount = 0;

	      i += (  9 + 128);
	      firstfound = 1;
	    
	    }
	  if ( onecount && hd->adc[i] >= 1220) 
	    {
	      onecount = 0;
	    }
	}
    }
  return 0;
}

int oncsSub_idsrs_v01::add_report ( hybriddata * hd, const int start, const int nreport)
{
  int i;

  report *r = new report;

  r->error = 1;
  r->n = nreport;

  if ( hd->adc[start + 11] < 1220) 
    {
      r->error = 0;
    }
  else if ( hd->adc[start + 11] > 3000)
    {
      r->error = 1;
    }
  else 
    {
      r->error = -1;  // funny bit
    }

  int acode = 0;
  for ( i = 0; i < 8; i++)
    {
      if (hd->adc[start + 3 + i] < 1220)
	{
	  acode = ( acode << 1)+1;
	}
      else if (hd->adc[ start + 3 + i] > 2980)
	{
	  acode = ( acode << 1);
	}
      else 
	{
	  // cout << "bit expected but got " << hd->adc[start + 3 + i] << " as bit " << i << endl;
	  acode = ( acode << 1);
	}
    }
  r->address = acode;


   for (i = 0; i < 128; i++)
   {
     int channel = 32 * ( i % 4 ) + 8 * ( i >> 2) - 31 * ( i >> 4);
     r->adc[channel] = hd->adc[start + 12+i];
   }
  hd->rowdata.push_back(r);

  return 0;
}


  
int oncsSub_idsrs_v01::iValue(const int ich,const char *what)
{

  if ( nhybrids == 0 ) decoded_data1 = decode(&data1_length);

  if ( strcmp(what,"NHYBRIDS") == 0 )
  {
    return hybridlist.size();
  }

  else if ( strcmp(what,"NSAMPLES") == 0 )
  {
    if ( ich >= 0 && ich < hybridlist.size() )
      {
	return hybridlist[ich]->rowdata.size();
      }
  }

  return 0;

}


int oncsSub_idsrs_v01::iValue(const int ich,const int hybrid, const char *what)
{

  if ( nhybrids == 0 ) decoded_data1 = decode(&data1_length);

  std::vector<hybriddata*>::iterator it;
  std::vector<report *>::iterator reportit;

  if ( ich < 0) return 0;

  if ( strcmp(what,"RAWSAMPLES") == 0 )
  {
	    
    for ( it = hybridlist.begin(); it != hybridlist.end(); ++it)
      {
	if ( (*it)->hdmi == hybrid )
	  {
	    if ( ich >= (*it)->words) return 0;
	    return (*it)->adc[ich];
	  }
      }
  }

    return 0;
}

int oncsSub_idsrs_v01::iValue(const int ich, const int tsample, const int hybrid)
{

  if ( ich < 0 || ich > 127) return 0;
  if ( nhybrids == 0 ) decoded_data1 = decode(&data1_length);
 
  std::vector<hybriddata*>::iterator it;
  std::vector<report *>::iterator reportit;

  for ( it = hybridlist.begin(); it != hybridlist.end(); ++it)
    {

      if ( (*it)->hdmi == hybrid)
	{

	  //cout << __LINE__ << "  " << __FILE__ << " tsample " << tsample << "  " << (*it)->rowdata.size() << endl;

	  if ( (*it)->rowdata.size()  < tsample +1) return 0;

	  report *r  = (*it)->rowdata[tsample];
	  return r->adc[ich];

	}
    }
  return 0;
}



void  oncsSub_idsrs_v01::dump ( OSTREAM& os )
{

  identify(os);

  int i;
  int is =  iValue(0,"NHYBRIDS");


  os << "Number of Hybrids: " << is << std::endl;
 
  std::vector<hybriddata*>::iterator it;
  std::vector<report *>::iterator reportit;

  for ( it = hybridlist.begin(); it != hybridlist.end(); ++it)
    {
      os << std::endl;
      os << "Framecounter: " << (*it)->framecounter << std::endl;
      os << "HDMI Channel: " << (*it)->hdmi << std::endl;
      os << "Description:  " << (*it)->desc << std::endl;
      os << "Words:        " << (*it)->words << std::endl;
      os << "Time samples: " << (*it)->rowdata.size() << std::endl;

//       for ( i = 0; i < (*it)->words; i++) 
// 	{
// 	  os << setw(6) << i << setw(6) << (*it)->adc[i];
// 	  if ( (*it)->adc[i] < 1200 ) cout << "  D1";
// 	  else if ( (*it)->adc[i] > 3000 ) cout << "  D0";
// 	  os << endl;
// 	}

      if ( (*it)->rowdata.size() > 0) 
	{
	  os << "Sample ";
	  for ( reportit = (*it)->rowdata.begin(); reportit != (*it)->rowdata.end(); ++reportit)
	    {
	      os << setw(5) << (*reportit)->n;
	    }
	  os << endl;
	  os << "Addr   ";

	  for ( reportit = (*it)->rowdata.begin(); reportit != (*it)->rowdata.end(); ++reportit)
	    {
	      os << setw(5) << (*reportit)->address;
	    }
	  os << endl;
	  os << "Error  ";
	  
	  for ( reportit = (*it)->rowdata.begin(); reportit != (*it)->rowdata.end(); ++reportit)
	    {
	      os << setw(5) << (*reportit)->error;
	    }
	  os << endl;
	  os << "       ----------------------------------------------------------------------------------------" << endl;

	  
	  for ( i = 0; i < 128; i++)
	    { 
	      os << setw(4) << i << " | ";
	      
	      for ( reportit = (*it)->rowdata.begin(); reportit != (*it)->rowdata.end(); ++reportit)
		{
		  os << setw(5) << (*reportit)->adc[i];
		}
	      os << endl;
	    }
	}
    }
      
}

