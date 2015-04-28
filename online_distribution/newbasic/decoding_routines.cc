
// --------- all decoding routines -----------

// decode the ID4EVT data format


#include "decoding_routines.h"

int decode_id4evt( int iarr[]
		   ,int *SubeventData
		   ,int dlength
		   ,int nlen
		   ,int *olength)
{

  if ( SubeventData == 0) return -1;

  int i;

  // we clear the output vector if NLEN is not 0 
  if (nlen > 0) 
    {
      for (i=0; i<nlen;) iarr[i++]=0;
    }
 
  int *sptr = SubeventData;
 
  int nrl = 0;
  for (i=0; i < dlength ; i++)
    {
      /*
	test if i is greater than NLEN; we exceed the
	allowed space in ARR then
      */
      if (nlen > 0 &&  i >= nlen)
	{
	  *olength = i+1;
	  return -1;
	}
      iarr[i] = *sptr++;
      nrl = i;
    }
  *olength = nrl+1;
  return 0;
}



// decode the ID2EVT data format

int decode_id2evt( int iarr[]
		   ,short *SubeventData
		   ,int dlength
		   ,int nlen
		   ,int *olength)
{
  if ( SubeventData == 0) return -1;
  int i;
  // we clear the output vector if NLEN is not 0 
  if (nlen > 0) 
    {
      for (i=0; i<nlen;) iarr[i++]=0;
    }

  short *sptr = SubeventData;
 
  int nrl;
  for (i=0; i < dlength ; i++)
    {
      /*
	test if i is greater than NLEN; we exceed the
	allowed space in ARR then
      */
      if (nlen > 0 &&  i >= nlen)
	{
	  *olength = i+1;
	  return -1;
	}
      iarr[i] = *sptr++;
      nrl = i;
    }
  *olength = nrl+1;
  return 0;
}

// the hammond device format
int decode_idhammond( int iarr[]
		      ,int *SubeventData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  if ( SubeventData == 0) return -1;
  int i;

  //COUT << "in std::decode_idhammond" << std::endl;

  // we clear the output vector if NLEN is not 0 
  if (nlen > 0) 
    {
      for (i=0; i<nlen;) iarr[i++]=0;
    }
 
  int *sptr = SubeventData;
 
  int out_channel = 0;
  int is, k, h1, h2;

  for (i=0; i < dlength ; i+=2)
    {
      /*
	test if i is greater than NLEN; we exceed the
	allowed space in ARR then
      */


      h1 = *sptr++;
      h2 = *sptr++;

      for ( k = 0; k<30; k+=5)
	{
	  if (nlen > 0 &&  out_channel >= nlen)
	    {
	      *olength = out_channel+1;
	      return -1;
	    }
	  iarr[out_channel++] = ((h1>>k) & 0x1f) ;
	}
      
      is = ((h1>>30) & 3) | ( h2 << 2);

      for (k = 0; k<10; k+=5)
	{
	  if (nlen > 0 &&  out_channel >= nlen)
	    {
	      *olength = out_channel+1;
	      return -1;
	    }
	  iarr[out_channel++] = ((is>>k) & 0x1f) ;
	}

      //----------------------
    }
  *olength = out_channel;
  return 0;
}

// the hammondsetup device format
int decode_idhammondset( int iarr[]
			 ,int *SubeventData
			 ,int dlength
			 ,int nlen
			 ,int *olength)
{
  if ( SubeventData == 0) return -1;
  int i;

  // we clear the output vector if NLEN is not 0 
  if (nlen > 0) 
    {
      for (i=0; i<nlen;) iarr[i++]=0;
    }
 
  int *sptr = SubeventData;
 
  int nrl = 0;
  for (i=0; i < dlength ; i++)
    {
      /*
	test if i is greater than NLEN; we exceed the
	allowed space in ARR then
      */
      if (nlen > 0 &&  i >= nlen)
	{
	  *olength = i+1;
	  return -1;
	}
      iarr[i] = *sptr++;
      nrl = i;
    }
  *olength = nrl+1;
  return 0;
}

int decode_idtecfem( int iarr[]
		     ,int *SubeventData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  if ( SubeventData == 0) return -1;
  int i, istart;

  // we clear the output vector if NLEN is not 0 
  if (nlen > 0) 
    {
      for (i=0; i<nlen;) iarr[i++]=0;
    }
 
  int *sptr = SubeventData;
 
  i = 0;

  istart = 8;

  int j;
  int b20;
  int nrl = 0;
  for (i=istart; i < dlength ; i++)
    {
      /*
        test if i is greater than NLEN; we exceed the
        allowed space in ARR then
      */

      //      COUT << "in loop, i= " << i << "  " << nrl <<" " 
      //   << std::hex << sptr[i] << std::dec <<  std::endl;

      if (nlen > 0 &&  nrl >= nlen)
        {
          *olength = nrl;
          return -1;
        }

      b20 = sptr[i];  // these are the 20 bits in one 32 bit word.

      // then bits 0-4...
      j = 0;
      iarr[nrl++] = (b20 >>j ) & 0x1f;

      // and finally bits 5-9.
      j=5;      
      iarr[nrl++] = (b20 >>j ) & 0x1f;

      // get bits 10-14 first
      j =10;
      iarr[nrl++] = (b20 >>j ) & 0x1f;

      // then bits 15-19...
      j =15;
      iarr[nrl++] = (b20 >>j ) & 0x1f;

        
    }
  *olength = nrl;
  return 0;
}


#include <mizar.h>

int decode_idmiznhc( int parr[]
                     ,int *SubeventData
                     ,int dlength
                     ,int nlen
                     ,int *olength)
{
  if ( SubeventData == 0) return -1;
  int i;


  // we clear the output vector if NLEN is not 0 
  if (nlen > 0) 
    {
      for (i=0; i<nlen;) parr[i++]=0;
    }
 
  miz_irdg iarr = (  miz_irdg ) parr;

  int ipos,isdm,index,adr,len,len_adc;
  sdm_c_block sdm;
  
  ipos = 2;
  for(isdm=0;isdm<11;isdm++)
    {
      sdm = (sdm_c_block) &SubeventData[ipos];
      len = (sdm->sdmlength & 0x0000ffff);
      //      COUT << "sdm " << isdm << " len " << len << " ipos " << ipos << std::endl;
      
      if ((sdm->sdmlength & 0x20000000))   /* this indicates the byte error */
        {
          iarr->out[isdm].byte_err = 1;
        }
      else
        {
          iarr->out[isdm].write_cell = (sdm->conv1_info & 0xf0000000)>>28;
          iarr->out[isdm].c1_cell    = (sdm->conv1_info & 0x0f000000)>>24;
          iarr->out[isdm].c2_cell    = (sdm->conv2_info & 0x0f000000)>>24;
          iarr->out[isdm].board_adr  = (sdm->conv1_info & 0x0000ff00)>>8;
          iarr->out[isdm].ser_ret    = (sdm->conv1_info & 0x000000ff);
 
          /*   iarr->out[isdm].socket     = (sdm->dspmap & 0x000000ff);
               iarr->out[isdm].port       = (sdm->dspmap & 0x0000ff00)>>8;
               iarr->out[isdm].dsp        = (sdm->dspmap & 0x00ff0000)>>16;
          */
          iarr->out[isdm].words      = len;
 
          /*     len_adc = len -4;  */
          len_adc = len -6;
          index = 0;
          while (index<len_adc)
            {
              if (sdm->array[index] & 0x80000000)
                {

                  adr = ( (sdm->array[index]>>20) & 0x7ff ) ;
                  iarr->out[isdm].conv1_high[adr] = (sdm->array[index] & 0x3ff);
                  iarr->out[isdm].conv1_low[adr]  = ((sdm->array[index]>>10) & 0x3ff);
                  iarr->out[isdm].conv2_high[adr] = (sdm->array[index+1] & 0x3ff);
                  iarr->out[isdm].conv2_low[adr]  = ((sdm->array[index+1]>>10) & 0x3ff);
                  iarr->out[isdm].tac[adr] =        ((sdm->array[index+1]>>20) & 0x3ff);
                  index += 2;
                }
 
            }                   /* end while index < len_adc */
        }                   /* end if byte err is set */
      ipos += len;
    }                     /* end for isdm */

  *olength  =  sizeof(*iarr)/4;
  return 0;
}

// the SAM5305 device format
int decode_idsam( int iarr[]
		  ,int *SubeventData
		  ,int dlength
		  ,int nlen
		  ,int *olength)
{
  if ( SubeventData == 0) return -1;
  int i, istat;

  // the format is essentially the ID4EVT format, the only
  // difference is that we have to convert the number into
  // millivolts (divide by 20).
  // we call decode_id4evt to do the work for us
  istat = decode_id4evt(iarr, SubeventData, dlength, nlen, olength);

  // and then we scale by 1/20.
  for (i=0; i < *olength ; i++)
    {
      iarr[i] /= 20;
    }

  // we return the status of  decode_id4evt. 
  return istat;
}

int decode_iddcfem( int iarr[]
		    ,int *SubeventData
		    ,int dlength
		    ,int nlen
		    ,int *olength)
{
  int i, istart;

  // we clear the output vector if NLEN is not 0 
  if (nlen > 0) 
    {
      for (i=0; i<nlen;) iarr[i++]=0;
    }
 
  int *sptr = SubeventData;
 
  i = 0;
  while (sptr[i] != 0xdc111 &&  i++ < dlength ){} 
  //  COUT << "index: " << i << "value " << std::hex << sptr[i] << " " << sptr[i+1] 
  //   << std::dec << std::endl;
  //  COUT << "dlength = " << dlength << "nlen= " << nlen << std::endl;
  istart = i + 5;

  int j;
  int b20;
  int nrl = 0;
  for (i=istart; i < dlength ; i++)
    {
      /*
	test if i is greater than NLEN; we exceed the
	allowed space in ARR then
      */

      //      COUT << "in loop, i= " << i << "  " << nrl <<" " 
      //   << std::hex << sptr[i] << std::dec <<  std::endl;

      if (nlen > 0 &&  nrl >= nlen)
	{
	  *olength = nrl;
	  return -1;
	}

      if (sptr[i] == 0xff444)
	{
	  *olength = nrl;
	  return 0;
	}

      b20 = sptr[i];  // these are the 20 bits in one 32 bit word.

      // get bits 10-14 first
      j =10;
      iarr[nrl++] = (b20 >>j ) & 0x1f;
			
      // then bits 15-19...
      j =15;
      iarr[nrl++] = (b20 >>j ) & 0x1f;

      // then bits 0-4...
      j = 0;
      iarr[nrl++] = (b20 >>j ) & 0x1f;

      // and finally bits 5-9.
      j=5;	
      iarr[nrl++] = (b20 >>j ) & 0x1f;
	
    }

  *olength = nrl;
  return 0;
}


int decode_bbc_dcm0( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_mvd_dcm0( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_dch_dcm0( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_pc_dcm0( int iarr[]
		    ,int *packetData
		    ,int dlength
		    ,int nlen
		    ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_tec_dcm0( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{

  if ( packetData == 0) return -1;
  int i;

  // we clear the output vector if NLEN is not 0 
  if (nlen > 0) 
    {
      for (i=0; i<nlen;) iarr[i++]=0;
    }
 
  int *sptr = &packetData[5];
 

  //  COUT << "index: " << i << "value " << std::hex << sptr[i] << " " << sptr[i+1] 
  //   << std::dec << std::endl;
  //  COUT << "dlength = " << dlength << "nlen= " << nlen << std::endl;

  int j;
  int b20;
  int nrl = 0;

  i = 0;
  //COUT << "first word " << std::hex << sptr[i] << std::dec << std::endl;
  while ( (sptr[i] & 0x80000000)  == 0) 
    {
      //  COUT << i << "  " << std::hex << sptr[i] << std::dec << std::endl;
      if ( i >  dlength -5) break;
 	/*
	  test if i is greater than NLEN; we exceed the
	  allowed space in ARR then
	*/
	
	//      COUT << "in loop, i= " << i << "  " << nrl <<" " 
	//   << std::hex << sptr[i] << std::dec <<  std::endl;
	
	b20 = sptr[i] & 0xfffff;  // these are the 20 bits in one 32 bit word.
	
	int upperbits =  (sptr[i] >> 20 ) & 0xfff;
	int wordnr    = upperbits & 0x01f;
	int channel  = (upperbits >> 5 ) & 0x3f;
	
	//if (channel < 5 ) COUT << "found channel" << channel << " wordnr " << wordnr << std::endl;
	nrl = channel * 80 + wordnr *4;
	// get bits 0...4 first
	j = 0;
	iarr[nrl++] = (b20 >>j ) & 0x1f;
	//if (channel < 5 ) COUT << "-- " << channel << "  " << nrl-1<< "  " << iarr[nrl-1] << std::endl;
	
	// then bits 5,,,,9...
	j = 5;
	iarr[nrl++] = (b20 >>j ) & 0x1f;
	//if (channel < 5 ) COUT << "-- " << channel << "  " <<  nrl-1 << "  " << iarr[nrl-1] << std::endl;
	
	// then bits 10.. 14...
	j = 10;
	iarr[nrl++] = (b20 >>j ) & 0x1f;
	//if (channel < 5 ) COUT << "-- " << channel << "  " <<nrl-1 << "  " << iarr[nrl-1] << std::endl;
	
	// and finally bits 15-19.
	j = 15;	
	iarr[nrl++] = (b20 >>j ) & 0x1f;
	//if (channel < 5 ) COUT << "-- " << channel << "  " << nrl-1 << "  " << iarr[nrl-1] << std::endl;
	
	i++;
	
    
    }

  *olength = 64 * 80 ;
  return 0;

}


int decode_rich_dcm0( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_tof_dcm0( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

// ------------------------------------------------------------
int decode_pbsc_dcm0( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{

  if ( packetData == 0) return -1;
  typedef struct {
    int post;
    int pre;
    int timing;
  } *emchannel;

  int i;

  // we clear the output vector if NLEN is not 0 
  for (i=0; i<nlen;) iarr[i++]=0;
 
  emchannel emc = (emchannel) &packetData[8];
 
  int channel;

  for (i=0; i < 144 ; i++)
    {
      /*
	test if i is greater than NLEN; we exceed the
	allowed space in ARR then
      */
      
      channel= ( emc->timing >> 20) & 0xff;

      if (nlen > 0 &&  5*channel >= nlen)
	{
	  COUT << "too short" << ENDL;
	  *olength = (channel+1)*5;
	  return -1;
	}


      if ( (emc->timing & 0x90000  && emc->post &0xc0000 && emc->pre & 0xa0000 ) 
	   &&  channel >= 0 && channel < 144 ) 
	{
      
	  if ( emc->post & 0x8000) // high!
	    {
	      iarr[channel] = (emc->timing & 0xfff);
	      iarr[channel+1] = 0;
	      iarr[channel+2] = (emc->post & 0xfff);
	      iarr[channel+3] = 0;
	      iarr[channel+4] = (emc->pre & 0xfff);
	    }
	  
	  else
	    {
	      iarr[channel] = (emc->timing & 0xfff);
	      iarr[channel+1] = (emc->post & 0xfff);
	      iarr[channel+2] = 0;
	      iarr[channel+3] = (emc->pre & 0xfff);
	      iarr[channel+4] = 0;
	    }
	}

      emc++;
    }
  *olength = 144*5;
  return 0;
}
//------------------------------------------------------------
int decode_pbsc_dcm32( int iarr[]
		       ,int *packetData
		       ,int dlength
		       ,int nlen
		       ,int *olength)
{
  if ( packetData == 0) return -1;

  typedef struct {
    int timing;
    int higain_post;
    int logain_post;
    int higain_pre;
    int logain_pre;
  } *emchannel;

  int i;

  // we clear the output vector if NLEN is not 0 
  if (nlen > 0) 
    {
      for (i=0; i<nlen;) iarr[i++]=0;
    }
 
  emchannel emc = (emchannel) &packetData[9];
 
  int nrl = 0;
  int j = 0;

  for (i=0; i < 192 ; i++)
    {
      /*
	test if i is greater than NLEN; we exceed the
	allowed space in ARR then
      */
      if (nlen > 0 &&  j >= nlen)
	{
	  COUT << "too short" << ENDL;
	  *olength = j+1;
	  return -1;
	}
			
      iarr[j++] = (emc->timing & 0xfff);
      iarr[j++] = (emc->higain_post & 0xfff);
      iarr[j++] = (emc->logain_post & 0xfff);
      iarr[j++] = (emc->higain_pre & 0xfff);
      iarr[j++] = (emc->logain_pre & 0xfff);

      nrl = j;
      emc++;
    }
  *olength = nrl+1;
  return 0;
}

// ------------------------------------------------------------

int decode_pbgl_dcm0( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_muta_dcm0( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


// ------------------------------------------------------------
// ------------------------------------------------------------
// modefied for mutr 2023 words format
int decode_mutc_dcm0( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  if ( packetData == 0) return -1;

  int i;

  // we clear the output vector if NLEN is not 0 
  if (nlen > 0) 
    {
      for (i=0; i<nlen;) iarr[i++]=0;
    }
 
  //  int *ic =  &packetData[10];   
  int *ic =  &packetData[5];   // start from k[5], missing 4 cells and one header word   
 
  int j = 0;

  //  if (nlen > 0 &&  nlen < 4*128)
  if (nlen > 0 &&  nlen < 63*32)    // 63 samples , 32 channels
    {
      COUT << "too short" << ENDL;
      *olength = 0;
      return -1;
    }

  //  for (i=0; i < 4*128 ; i++)
  for (i=0; i < 63*32 ; i++)
    {
			
      iarr[j++] = ic[i];

    }
  //  *olength = 4*128;
    *olength = 63*32;
  return 0;

}



int decode_muid_dcm0( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_zdc_dcm0( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_rich_ll1( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_mvd_ll1( int iarr[]
		    ,int *packetData
		    ,int dlength
		    ,int nlen
		    ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_bbc_ll1( int iarr[]
		    ,int *packetData
		    ,int dlength
		    ,int nlen
		    ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_ntczdc_ll1( int iarr[]
		    ,int *packetData
		    ,int dlength
		    ,int nlen
		    ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_big_ll1( int iarr[]
		    ,int *packetData
		    ,int dlength
		    ,int nlen
		    ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_tof_ll1( int iarr[]
		    ,int *packetData
		    ,int dlength
		    ,int nlen
		    ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_muid_ll1( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_ert_ll1( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_pbgl_ll1( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_pbsc_ll1( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_gl1( int iarr[]
		,int *packetData
		,int dlength
		,int nlen
		,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_gl1p( int iarr[]
		,int *packetData
		,int dlength
		,int nlen
		,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_bbc_dcm1( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_mvd_dcm1( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_dch_dcm1( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_pc_dcm1( int iarr[]
		    ,int *packetData
		    ,int dlength
		    ,int nlen
		    ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_tec_dcm1( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_rich_dcm1( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_tof_dcm1( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_pbsc_dcm1( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_pbgl_dcm1( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_muta_dcm1( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}



// ----- modefied for mutr 532 words format, place holder, for MC
int decode_mutc_dcm1( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  if ( packetData == 0) return -1;
  int i;

  // we clear the output vector if NLEN is not 0 
  if (nlen > 0){
      for (i=0; i<nlen;) iarr[i++]=0;
    }
 
  //  int *ic =  &packetData[10];   
  int *ic =  &packetData[10];   // start from k[10]   
 
  int j = 0;

 
  //  if (nlen > 0 &&  nlen < 128)    // 1 sample , 128 channels
  if (nlen > 0 &&  nlen < 512){      // 4 samples, 128 channels  
  
    COUT << "too short" << ENDL;
    *olength = 0;
    return -1;
  }

  // for (i=0; i < 128 ; i++)
  for (i=0; i < 512 ; i++) {			
    iarr[j++] = ic[i];
  }
  *olength = 512;
  //  *olength = 128;
  return 0;
}

int decode_muid_dcm1( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_zdc_dcm1( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_bbc_dcm2( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_mvd_dcm2( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_dch_dcm2( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_pc_dcm2( int iarr[]
		    ,int *packetData
		    ,int dlength
		    ,int nlen
		    ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_pc_dcm3( int iarr[]
		    ,int *packetData
		    ,int dlength
		    ,int nlen
		    ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_tec_dcm2( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_rich_dcm2( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_tof_dcm2( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_pbsc_dcm2( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_pbgl_dcm2( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_muta_dcm2( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

// ----- modefied for mutr 532 words format, pass-through mode
int decode_mutc_dcm2( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  if ( packetData == 0) return -1;
  int i;

  // we clear the output vector if NLEN is not 0 
  if (nlen > 0){
    for (i=0; i<nlen;) iarr[i++]=0;
  }
 
  //  int *ic =  &packetData[10];   
  int *ic =  &packetData[10];   // start from k[10]   
 
  int j = 0;

  if (nlen > 0 &&  nlen < 512){    // 4 samples , 128 channels
    COUT << "too short" << ENDL;
    *olength = 0;
    return -1;
  }

  for (i=0; i < 512 ; i++) {	
    iarr[j++] = ic[i];
  }

  *olength = 512;
  return 0;

}

int decode_muid_dcm2( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}

int decode_zdc_dcm2( int iarr[]
		     ,int *packetData
		     ,int dlength
		     ,int nlen
		     ,int *olength)
{
  return decode_id4evt( iarr, packetData, dlength, nlen, olength);
}


int decode_pbsc_dcms( int iarr[]

		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{

  if ( packetData == 0) return -1;
  typedef struct 
  {
    int timing;
    int post;
    int pre;
  } *emchannel;
  
  int i;

  // we clear the output vector if NLEN is not 0 
  if (nlen > 0) 
    {
      for (i=0; i<nlen;) iarr[i++]=0;
    }
 
  emchannel emc = (emchannel) &packetData[9];
 
  int nrl = 0;
  int j = 0;

  for (i=0; i < 144 ; i++)
    {
      /*
	test if i is greater than NLEN; we exceed the
	allowed space in ARR then
      */
      if (nlen > 0 &&  j >= nlen)
	{
	  COUT << "too short" << ENDL;
	  *olength = j+1;
	  return -1;
	}
			
      if ( emc->post & 0x1000) // high!
	{
	  iarr[j++] = (emc->timing & 0xfff);
	  iarr[j++] = 0;
	  iarr[j++] = (emc->post & 0xfff);
	  iarr[j++] = 0;
	  iarr[j++] = (emc->pre & 0xfff);
	}

      else
	{
	  iarr[j++] = (emc->timing & 0xfff);
	  iarr[j++] = (emc->post & 0xfff);
	  iarr[j++] = 0;
	  iarr[j++] = (emc->pre & 0xfff);
	  iarr[j++] = 0;
	}

      nrl = j;
      emc++;
    }
  *olength = nrl+1;
  return 0;
}

// ----- modefied for mutr 532 words format, pass-through mode
int decode_mutc_dcm3( int iarr[]
		      ,int *packetData
		      ,int dlength
		      ,int nlen
		      ,int *olength)
{
  if ( packetData == 0) return -1;
  int i;

  // we clear the output vector if NLEN is not 0 
  if (nlen > 0){
    for (i=0; i<nlen;) iarr[i++]=0;
  }
 
  //  int *ic =  &packetData[10];   
  int *ic =  &packetData[10];   // start from k[10]   
 
  int j = 0;

  if (nlen > 0 &&  nlen < 512){    // 4 samples , 128 channels
    COUT << "too short" << ENDL;
    *olength = 0;
    return -1;
  }

  for (i=0; i < 512 ; i++) {	
    iarr[j++] = ic[i];
  }

  *olength = 512;
  return 0;

}




