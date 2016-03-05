#include <packet_gl1p.h>
#include <string.h>

Packet_gl1p::Packet_gl1p(PACKET_ptr data)
  : Packet_w4(data)
{
  sgl1p=0;
}

Packet_gl1p::~Packet_gl1p()
{
  if(sgl1p)
    {
      delete [] sgl1p;
    }
}

int *Packet_gl1p::decode ( int *nwout)
{
  int *p,*k;
  int i, olength;
  int temp[MAX_OUTLENGTH];
  
  int dlength = getDataLength();

  int status = decode_gl1p( temp
			      ,(int *)  findPacketDataStart(packet) 
			      ,dlength
			      ,MAX_OUTLENGTH, &olength);

  if (status || olength<=0 ) return NULL;
 
  p = new int[olength];
  k = p;
  for (i =0; i<olength; i++) *k++ = temp[i];
  *nwout = olength;
  return p;
}
int Packet_gl1p::iValue(const int ich, const int what)
{// int what is used hear instead of char* what to reduce data taking time. 
  if(!sgl1p) demangle();
  if(!sgl1p){
    COUT<<"Failed to fill data structure"<<std::endl;
    return 0;
  }
  switch (what)
    {
    case GL1P_HEADER:
      {
	if(ich<0 || ich>= nGL1Pboards) return 0;
	else return sgl1p[ich].header;
      }
    case GL1P_EVNUMBER:
      {
	if(ich<0 || ich>= nGL1Pboards) return 0;
	else return sgl1p[ich].ev_number;
      }
    case GL1P_MODEBIT:
      {
	if(ich<0 || ich>= nGL1Pboards) return 0;
	else return sgl1p[ich].modebit;
      }
    case GL1P_CLOCK:	
      {// it is supposed that the number of "clock counters" is 4*nGL1Pboards.
	//For example, to get Clock Counter for Scaler A of GL1P board 2 call iValue(5,GL1P_CLOCK)
	if(ich<0 || ich>= 4*nGL1Pboards) return -1;
	else return sgl1p[ich/4].clock[ich%4];
      }
    case GL1P_SCALER:	
      {//Scaler value. The same rule about ich as above.
	if(ich<0 || ich>= 4*nGL1Pboards) return 0;
	else return sgl1p[ich/4].scaler[ich%4];
      }
    default:
      return 0;
    }
}
int Packet_gl1p::iValue(const int ich, const char *what)
{

  
  if(!sgl1p) demangle();
  if(!sgl1p){
    COUT<<"Failed to fill data structure"<<std::endl;
    return 0;
  }

  if (!strcmp(what,"HEADER"))
    {
	if(ich<0 || ich>= nGL1Pboards) return 0;
	else return sgl1p[ich].header;
    }

  else if (!strcmp(what,"EVNUMBER"))
    {
	if(ich<0 || ich>= nGL1Pboards) return 0;
	else return sgl1p[ich].ev_number;
    }
  else if (!strcmp(what,"MODEBIT"))
    {
	if(ich<0 || ich>= nGL1Pboards) return 0;
	else return sgl1p[ich].modebit;
    }
  else if (!strcmp(what,"CLOCK"))
    {
	if(ich<0 || ich>= 4*nGL1Pboards) return -1;
	else return sgl1p[ich/4].clock[ich%4];
    }
  else if (!strcmp(what,"SCALER"))
    {
	if(ich<0 || ich>= 4*nGL1Pboards) return 0;
	else return sgl1p[ich/4].scaler[ich%4];
    }
  else return 0;

}
/** function fillIntArray fills destination[4*nGL1Pboards] if what = "CLOCK"
or "SCALER", and destination[nGL1Pboards] if what = "HEADER", "EVNUMBER" or
"MODEBIT" */
int Packet_gl1p::fillIntArray (int destination[],    // the data go here 
			       const int length,      // space we have in destination
			       int * nw,              // words actually used
			       const char * what) // type of data (see above)
{
  int i;
  if(!sgl1p) demangle();
  if(!sgl1p){
    COUT<<"Failed to fill data structure"<<std::endl;
    return -1;
  }
  if (!strcmp(what,"CLOCK"))
    {
      if(length>4*nGL1Pboards)*nw=4*nGL1Pboards;
      else *nw=length;
      for(i=0;i<*nw;i++)
	destination[i]=sgl1p[i/4].clock[i%4];
    }
  else if (!strcmp(what,"SCALER"))
    {
      if(length>4*nGL1Pboards)*nw=4*nGL1Pboards;
      else *nw=length;
      for(i=0;i<*nw;i++)
	destination[i]=sgl1p[i/4].scaler[i%4];
    }
  else if (!strcmp(what,"HEADER"))
    {
      if(length>nGL1Pboards)*nw=nGL1Pboards;
      else *nw=length;
      for(i=0;i<*nw;i++)
	destination[i]=sgl1p[i].header;
    }
  else if (!strcmp(what,"EVNUMBER"))
    {
      if(length>nGL1Pboards)*nw=nGL1Pboards;
      else *nw=length;
      for(i=0;i<*nw;i++)
	destination[i]=sgl1p[i].ev_number;
    }
  else if (!strcmp(what,"MODEBIT"))
    {
      if(length>nGL1Pboards)*nw=nGL1Pboards;
      else *nw=length;
      for(i=0;i<*nw;i++)
	destination[i]=sgl1p[i].modebit;
    }
  else
    {
      *nw=0;
    }
  return 0;
}

void Packet_gl1p::dump ( OSTREAM &os)
{
  int i,j,l,m,dlength;
  dlength = getDataLength();

	this->identify(os); 
	if(!sgl1p) demangle();
	if(!sgl1p){
	  os<<"Failed to fill sgl1p. Exit"<<std::endl;
	  return;
	}

	for(m=0;m<54;m++) os<<"=";
	os << std::endl;
	os << "GL1P data packet: ";
	os << "Detected " << std::dec << nGL1Pboards << " GL1-1P boards in the data packet." <<std::endl;
	for(i=0;i<nGL1Pboards;i++){
	  for(m=0;m<54;m++) os<<"-";
	  COUT<<std::endl;
	  os <<std::dec<<"GL1P["<<i<<"] header word    = 0x" << std::hex << SETW(4) << (unsigned int)sgl1p[i].header << std::endl;
	  os <<std::dec<<"Event = " <<SETW(4)<<(unsigned int)sgl1p[i].ev_number<<std::endl;
	  os <<"Modebit    = 0x" << std::hex << SETW(3) << sgl1p[i].modebit << std::endl;
	  os << "Beam Cross. Counters: ";
	  for(j=0;j<4;j++)
	    os <<std::dec<<SETW(12)<<(unsigned int)sgl1p[i].clock[j];
	  os<<std::endl;
	  os << "Scalers             : ";
	  for(j=0;j<4;j++)
	    os<<SETW(12)<<sgl1p[i].scaler[j];
	  os<<std::endl;
	}
	for(m=0;m<54;m++) os<<"-";
	// Raw data 
        j = 0;
	unsigned int* k=(unsigned int *) findPacketDataStart(packet);
	if (k == 0) 
	  {
	    os << std::endl;
	    return;
	  }

        while (1)
	{
	  os << std::endl << std::dec << SETW(5) << j << " |  ";
	  for (l=0;l<6;l++)
	  {
	      os << std::hex << SETW(8) << k[j++] << " ";
	      if (j>=dlength) break;
	  }
	  if (j>=dlength) break;
	}	
        os << std::endl;
	for(m=0;m<54;m++) os<<"=";
        os << std::dec << std::endl;

}

void Packet_gl1p::demangle()
{
  int sclk[]={1,0,3,2};
  /* 
     The job of this routine is to demangle the data into something
     that can be more properly addressed - the data is copied
     into a GL1 data structure that can be addressed properly.
  */
  int dlength = getDataLength();
  nGL1Pboards=dlength/6;
  if(sgl1p) delete sgl1p;
  if(nGL1Pboards<1) {
    sgl1p=0;
    return;
  }
  sgl1p=new GL1P_DATA[nGL1Pboards];
  if(!sgl1p){
    COUT<<"can't allocate memory for GL1P_DATA structure "<<std::endl;
    return;
  }
  unsigned int* buf = (unsigned int *) findPacketDataStart(packet);
  if (buf == 0) 
    {
      delete [] sgl1p;
      sgl1p = NULL;
      return;
    }
  // unpack the data into a GL1 structure:


  unsigned int  j, idx, idy;
  unsigned char bclock;
  int i;

  for(i=0;i<nGL1Pboards;i++){
    idx=i*6;
    sgl1p[i].header    = (unsigned char)((buf[idx]>>8)&0xFF);
    sgl1p[i].ev_number = (unsigned char)(buf[idx]&0xFF);
    sgl1p[i].modebit   = (unsigned short)((buf[idx]>>16)&0x1FF);
      /* Note the complicated way in which the 20-bit LUT inputs are packed */
    for(j=0;j<4;j++){
      idy=j*8;
      bclock = (unsigned char)((buf[idx+1]>>idy)&0x7F);
      sgl1p[i].clock[sclk[j]] =(bclock&0xF) + ((bclock>>4)&0x7)*15;
      sgl1p[i].scaler[j] = ((buf[idx+2+j]&0xFFFF)<<16) + ((buf[idx+2+j]>>16)&0xFFFF);
    }

  }

}
