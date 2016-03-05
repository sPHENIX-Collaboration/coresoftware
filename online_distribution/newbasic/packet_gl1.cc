#include <packet_gl1.h>
#include <string.h>

Packet_gl1::Packet_gl1(PACKET_ptr data)
  : Packet_w4 (data){sgl1=0;}

Packet_gl1::~Packet_gl1()
{
  if(sgl1)delete sgl1;
}

int *Packet_gl1::decode ( int *nwout)
{
  int *p,*k;
  int olength;
  int temp[MAX_OUTLENGTH];
  int i;
  int dlength = getDataLength();

  int status = decode_gl1( temp
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

int Packet_gl1::iValue(const int ich, const char *what)
{

  
  if(!sgl1) demangle();
  if(!sgl1){
    COUT<<"Failed to fill data structure"<<std::endl;
    return 0;
  }
  // GL1-3
  if (strcmp(what,"HEADER3")==0)
    {
      return sgl1->gl3_payload.gl3_header;
    }
  /* return time */
  else if ( strcmp(what,"YEAR")==0)
    {
      return sgl1->gl3_payload.timestamp.year;
    }
  else if ( strcmp(what,"MONTH")==0)
    {
      return sgl1->gl3_payload.timestamp.month;
    }
  else if ( strcmp(what,"DATE")==0)
    {
      return sgl1->gl3_payload.timestamp.date;
    }
  else if ( strcmp(what,"DAY")==0)
    {
      return sgl1->gl3_payload.timestamp.day;
    }
  else if ( strcmp(what,"HOUR")==0)
    {
      return sgl1->gl3_payload.timestamp.hour;
    }
  else if ( strcmp(what,"MIN")==0)
    {
      return sgl1->gl3_payload.timestamp.min;
    }
  else if ( strcmp(what,"SEC")==0)
    {
      return sgl1->gl3_payload.timestamp.sec;
    }

  else if (strcmp(what,"ALIGNMENT")==0)
    {
      return sgl1->gl3_payload.alignment;
    }
  else if ( strcmp(what,"CROSSCTR")==0)
    {
      return sgl1->gl3_payload.bunch_crossing_counter;
    }
  else if ( strcmp(what,"BEAMCTR0")==0)
    {
      return sgl1->gl3_payload.beam_crossing_counter[0];
    }
  else if ( strcmp(what,"BEAMCTR1")==0)
    {
      return sgl1->gl3_payload.beam_crossing_counter[1];
    }
  else if ( strcmp(what,"GACCEPT")==0)
    {
      return sgl1->gl3_payload.granule_accept_vector;
    }
  else if (strcmp(what,"ACPTORINP")==0)
    {
      return sgl1->gl3_payload.accept_or_input;
    }
  else if ( strcmp(what,"ACPTCTR")==0)
    {
      return sgl1->gl3_payload.gl1_accept_counter;
    }
  else if ( strcmp(what,"GRANCTR")==0)
    {
      if( (ich<0) || (ich>31) ) {
	return 0;
      }
      else{
	return sgl1->gl3_payload.granule_accept[ich];
      }
    }
  else if ( strcmp(what,"GDISABLE")==0)
    {
      return sgl1->gl3_payload.granule_disables;
    }
 else if ( strcmp(what,"FACCEPT")==0)
    {
      return sgl1->gl3_payload.forced_accepts;
    }
  // GL1-2 
  else if (strcmp(what,"HEADER2")==0)
    {
      return sgl1->gl2_payload.gl2_header;
    }
  else if ( strcmp(what,"PACCEPT")==0)
    {
      return sgl1->gl2_payload.partition_accept;
    }
  else if (strcmp(what,"MODEBITS")==0)
    {
      return sgl1->gl2_payload.mode_bits;
    }
  else if ( strcmp(what,"RBITS0")==0)
    {
      return sgl1->gl2_payload.reduced_bits[0];
    }
  else if ( strcmp(what,"RBITS1")==0)
    {
      return sgl1->gl2_payload.reduced_bits[1];
    }
  else if ( strcmp(what,"DCMFULL")==0)
    {
      return sgl1->gl2_payload.dcm_full_fem_busy;
    }
  else if ( strcmp(what,"FEMUNREL")==0)
    {
      return sgl1->gl2_payload.fem_unreliable;
    }
  else if ( strcmp(what,"GBUSY")==0)
    {
      return sgl1->gl2_payload.granule_busy;
    }
  else if ( strcmp(what,"PXBAR")==0)
    {
      return sgl1->gl2_payload.part_busy_xbar_out;
    }
  else if ( strcmp(what,"PBUSY")==0)
    {
      return sgl1->gl2_payload.part_busy_bus;
    }

  // GL1-1
  else if (strcmp(what,"HEADER1")==0)
    {
      if( (ich<0) || (ich>=sgl1->gl1_boards) ) {
	return 0;
      }
      else{
	return sgl1->gl1_payload[ich].gl1_header;
      }
    }
  else if ( strcmp(what,"LUTINPUT")==0)
    {// it is supposed that lut input is defined as: ich= i*8+lutinput,
      // were i is GL1-1 board number
      if( (ich<0) || (ich>=8*sgl1->gl1_boards) ) {
	return 0;
      }
      else{
        return sgl1->gl1_payload[ich/8].lut_input[ich%8];
      }
    }
  else if ( strcmp(what,"RAWTRIG")==0)
    {
      if( (ich<0) || (ich>=sgl1->gl1_boards) ) {
	return 0;
      }
      else{
        return sgl1->gl1_payload[ich].lut_output;
      }
    }
  else if ( strcmp(what,"TRIGBUSY")==0)
    {
      if( (ich<0) || (ich>=sgl1->gl1_boards) ) {
	return 0;
      }
      else{
        return sgl1->gl1_payload[ich].trigger_busy;
      }
    }
  else if ( strcmp(what,"LIVETRIG")==0)
    {
      if( (ich<0) || (ich>=sgl1->gl1_boards) ) {
	return 0;
      }
      else{
        return sgl1->gl1_payload[ich].live_trig_out;
      }
    }
  else if ( strcmp(what,"SCALEDTRIG")==0)
    {
      if( (ich<0) || (ich>=sgl1->gl1_boards) ) {
	return 0;
      }
      else{
        return sgl1->gl1_payload[ich].scaled_trig_out;
      }
    }
  else if ( strcmp(what,"TRIGPARXBAR")==0)
    {
      if( (ich<0) || (ich>=sgl1->gl1_boards) ) {
	return 0;
      }
      else{
        return sgl1->gl1_payload[ich].trig_part_xbar_out;
      }
    }


  else return 0;

}
int Packet_gl1::iValue(const int ich, const  int what)
{

  
  if(!sgl1) demangle();
  if(!sgl1){
    std::cout<<"Failed to fill data structure"<<std::endl;
    return 0;
  }
  switch (what) {
  // GL1-3
  case HEADER3:
    {
      return sgl1->gl3_payload.gl3_header;
    }
  /* return time */
  case YEAR:
    {
      return sgl1->gl3_payload.timestamp.year;
    }
  case MONTH:
    {
      return sgl1->gl3_payload.timestamp.month;
    }
  case DATE:
    {
      return sgl1->gl3_payload.timestamp.date;
    }
  case DAY:
    {
      return sgl1->gl3_payload.timestamp.day;
    }
  case HOUR:
    {
      return sgl1->gl3_payload.timestamp.hour;
    }
  case MIN:
    {
      return sgl1->gl3_payload.timestamp.min;
    }
  case SECGL1:
    {
      return sgl1->gl3_payload.timestamp.sec;
    }

  case ALIGNMENT:
    {
      return sgl1->gl3_payload.alignment;
    }
  case CROSSCTR:
    {
      return sgl1->gl3_payload.bunch_crossing_counter;
    }
  case BEAMCTR0:
    {
      return sgl1->gl3_payload.beam_crossing_counter[0];
    }
  case BEAMCTR1:
    {
      return sgl1->gl3_payload.beam_crossing_counter[1];
    }
  case GACCEPT:
    {
      return sgl1->gl3_payload.granule_accept_vector;
    }
  case ACPTORINP:
    {
      return sgl1->gl3_payload.accept_or_input;
    }
  case ACPTCTR:
    {
      return sgl1->gl3_payload.gl1_accept_counter;
    }
  case GRANCTR:
    {
      if( (ich<0) || (ich>31) ) {
	return 0;
      }
      else{
	return sgl1->gl3_payload.granule_accept[ich];
      }
    }
  case GDISABLE:
    {
      return sgl1->gl3_payload.granule_disables;
    }
 case FACCEPT:
    {
      return sgl1->gl3_payload.forced_accepts;
    }
  // GL1-2 
  case HEADER2:
    {
      return sgl1->gl2_payload.gl2_header;
    }
  case PACCEPT:
    {
      return sgl1->gl2_payload.partition_accept;
    }
  case MODEBITS:
    {
      return sgl1->gl2_payload.mode_bits;
    }
  case RBITS0:
    {
      return sgl1->gl2_payload.reduced_bits[0];
    }
  case RBITS1:
    {
      return sgl1->gl2_payload.reduced_bits[1];
    }
  case DCMFULL:
    {
      return sgl1->gl2_payload.dcm_full_fem_busy;
    }
  case FEMUNREL:
    {
      return sgl1->gl2_payload.fem_unreliable;
    }
  case GBUSY:
    {
      return sgl1->gl2_payload.granule_busy;
    }
  case PXBAR:
    {
      return sgl1->gl2_payload.part_busy_xbar_out;
    }
  case PBUSY:
    {
      return sgl1->gl2_payload.part_busy_bus;
    }

  // GL1-1
  case HEADER1:
    {
      if( (ich<0) || (ich>=sgl1->gl1_boards) ) {
	return 0;
      }
      else{
	return sgl1->gl1_payload[ich].gl1_header;
      }
    }
  case LUTINPUT:
    {// it is supposed that lut input is defined as: ich= i*8+lutinput,
      // were i is GL1-1 board number
      if( (ich<0) || (ich>=8*sgl1->gl1_boards) ) {
	return 0;
      }
      else{
        return sgl1->gl1_payload[ich/8].lut_input[ich%8];
      }
    }
  case RAWTRIG:
    {
      if( (ich<0) || (ich>=sgl1->gl1_boards) ) {
	return 0;
      }
      else{
        return sgl1->gl1_payload[ich].lut_output;
      }
    }
  case TRIGBUSY:
    {
      if( (ich<0) || (ich>=sgl1->gl1_boards) ) {
	return 0;
      }
      else{
        return sgl1->gl1_payload[ich].trigger_busy;
      }
    }
  case LIVETRIG:
    {
      if( (ich<0) || (ich>=sgl1->gl1_boards) ) {
	return 0;
      }
      else{
        return sgl1->gl1_payload[ich].live_trig_out;
      }
    }
  case SCALEDTRIG:
    {
      if( (ich<0) || (ich>=sgl1->gl1_boards) ) {
	return 0;
      }
      else{
        return sgl1->gl1_payload[ich].scaled_trig_out;
      }
    }
  case TRIGPARXBAR:
    {
      if( (ich<0) || (ich>=sgl1->gl1_boards) ) {
	return 0;
      }
      else{
        return sgl1->gl1_payload[ich].trig_part_xbar_out;
      }
    }
  default:
    {
      return 0;
    }
  }

  return 0;

}

void Packet_gl1::dump ( OSTREAM &os)
{
  int i,j,l,m,dlength;
  dlength = getDataLength();

	this->identify(os); 
	if(!sgl1) demangle();
	if(!sgl1){
	  os<<"Failed to fill sgl1. Exit"<<std::endl;
	  return;
	}

	// output some useful information:
	const char* days[]={
	  "Non",
	  "Sun",
	  "Mon",
	  "Tue",
	  "Wed",
	  "Thu",
	  "Fri",
	  "Sat"
	};
	std::ios::fmtflags oldFlags = os.flags();
	char oldFill;
	oldFill=os.fill('0');
	for(m=0;m<54;m++) os<<"=";
	os << std::endl;
	os << "GL1 data packet:" << std::endl;
	os << "Detected " << std::dec << sgl1->gl1_boards << " GL1-1 boards in the data packet." <<std::endl;
        os << "--> GL1-3 <--"<<std::endl;
	os << "GL1-3 header word    = 0x" << std::hex << SETW(4) << sgl1->gl3_payload.gl3_header << std::endl;
	os << "Time: " <<std::dec<<SETW(2)<<(sgl1->gl3_payload.timestamp.hour)<<":"<<SETW(2)<<(sgl1->gl3_payload.timestamp.min)<<":"
	   <<SETW(2)<<(sgl1->gl3_payload.timestamp.sec)<<" "<<SETW(2)<<days[(sgl1->gl3_payload.timestamp.day)]
	   <<", "<<(sgl1->gl3_payload.timestamp.month)<<"/"<<SETW(2)<<(sgl1->gl3_payload.timestamp.date)<<"/"
	   <<SETW(4)<<(sgl1->gl3_payload.timestamp.year)<<std::endl;
	os <<"Alignment = 0x"<<std::hex<<SETW(4)<<sgl1->gl3_payload.alignment<<std::endl;
        os << "Beam Cross. Counter  = 0x" << SETW(8) 
	   << sgl1->gl3_payload.beam_crossing_counter[1];
	os << SETW(8)
	   << sgl1->gl3_payload.beam_crossing_counter[0] << std::endl;
	os << "Bunch Cross. Counter = 0x" << SETW(4) << sgl1->gl3_payload.bunch_crossing_counter << std::endl;
	os << "Granule Accept Vector    = 0x" << SETW(8) << sgl1->gl3_payload.granule_accept_vector << std::endl;
	os << "Accept OR input          = 0x" << SETW(8) << sgl1->gl3_payload.accept_or_input << std::endl;
// 	os.flags(std::ios::right|(oldFlags^std::ios::left));
// 	os.fill(' ');
	os << "Accept Counter       = " << sgl1->gl3_payload.gl1_accept_counter << std::endl;
	os << "Granule Accept Counters:"<< std::endl;

	for(i=3;i>=0;i--)
	  {
	    for(j=7;j>=0;j--)
	      {
		l=i*8+j;
		if(l<10) os <<" ["<< std::dec<<SETW(1)<<l<<"]=0x"<<std::hex<<SETW(4)<<sgl1->gl3_payload.granule_accept[l]<<" ";
		else os <<"["<< std::dec<<SETW(2)<<l<<"]=0x"<<std::hex<<SETW(4)<<sgl1->gl3_payload.granule_accept[l]<<" ";
	      }
	    os<<std::endl;
	  }
	os.flags(oldFlags);
	os.fill('0');
	os << "Granule Disables         = 0x" << std::hex << SETW(8) << sgl1->gl3_payload.granule_disables << std::endl;
	os << "Forced Accepts           = 0x" << SETW(8) << sgl1->gl3_payload.forced_accepts << std::endl;
	os << "--> GL1-2 <--"<<std::endl;
	os << "GL1-2 header word    = 0x" << SETW(4) << sgl1->gl2_payload.gl2_header << std::endl;
	os << "Partition Accept Vector  = 0x" << SETW(8) << sgl1->gl2_payload.partition_accept << std::endl;
	os << "Mode bits = 0x"<< SETW(4) << sgl1->gl2_payload.mode_bits << std::endl;
	os << "Reduced bits             = 0x" << SETW(8)<<sgl1->gl2_payload.reduced_bits[1]
	   << sgl1->gl2_payload.reduced_bits[0] << std::endl;
        os << "DCM Full/FEM Busy Vector = 0x" << SETW(8) << sgl1->gl2_payload.dcm_full_fem_busy << std::endl;
        os << "FEM Unreliable Vector    = 0x" << SETW(8) << sgl1->gl2_payload.fem_unreliable << std::endl;
	os << "Granule busy             = 0x" << SETW(8) << sgl1->gl2_payload.granule_busy << std::endl;
	os << "Partition busy Xbar Out  = 0x" << SETW(8) << sgl1->gl2_payload.part_busy_xbar_out << std::endl;
	os << "Partition busy bus       = 0x" << SETW(8) << sgl1->gl2_payload.part_busy_bus << std::endl;
	for(i=0;i<sgl1->gl1_boards;i++)
	  {
	    os << std::dec << SETW(1) << "--> GL1-1[" << i << "] <--" << std::endl;
	    os << SETW(1) << "GL1-1[" << i << "] header word = 0x" << std::hex << SETW(4) << sgl1->gl1_payload[i].gl1_header << std::endl;
	    os << "LUT inputs ";
	    for(j=7;j>=0;j--)
	      os << "["<< std::dec << SETW(1) << j <<"]<-0x"<< std::hex << SETW(5) <<sgl1->gl1_payload[i].lut_input[j]<<" ";
	    os <<std::endl;
	    os << "LUT Out (Raw Triggers)   = 0x" << std::hex << SETW(8) << sgl1->gl1_payload[i].lut_output << std::endl;
	    os << "Trigger busy input       = 0x" << SETW(8) << sgl1->gl1_payload[i].trigger_busy << std::endl;
	    os << "Live Triggers            = 0x" << SETW(8) << sgl1->gl1_payload[i].live_trig_out << std::endl;
	    os << "Scaled Triggers          = 0x" << SETW(8) << sgl1->gl1_payload[i].scaled_trig_out << std::endl;
	    os << "Trig->Part Xbar Out      = 0x" << SETW(8) << sgl1->gl1_payload[i].trig_part_xbar_out << std::endl;
	    os << std::dec << std::endl;
	  }

	// Finally, a generic HEX dump:

        j = 0;
	int* k=( int *) findPacketDataStart(packet);
        if ( k ==0 ) return;

	while (1)
	  {
	    os << std::endl << std::dec << SETW(5) << j << " |  ";
	    for (l=0;l<4;l++)
	      {
		os << std::hex << SETW(8) << k[j++] << " ";
		if (j>=dlength) break;
	      }
	    if (j>=dlength) break;
	  }	
        os << std::endl;
	for(m=0;m<54;m++) os<<"=";
        os << std::dec << std::endl;
	os.fill(oldFill);
}

void Packet_gl1::demangle()
{

  /* 
     The job of this routine is to demangle the data into something
     that can be more properly addressed - the data is copied
     into a GL1 data structure that can be addressed properly.
  */
  int  check_length;
  if(sgl1) delete sgl1;
  sgl1=new GL1_EVENT_DATA;
  if(!sgl1){
    COUT<<"can't allocate memory for GL1_EVENT_DATA structure "<<std::endl;
    return;
  }
  int dlength = getDataLength();

  // First board will have 21 16-bit words in readout to keep the total number even at 46.
  // If there is an even number of GL1-1 boards total (or an odd number of additional 
  // GL1-1 boards) then the last GL1-1 board will have to read out 22
  // 16-bit words to keep the total an intergral number of 16-bit words.
	
  sgl1->gl1_boards = 1;
  check_length = 2*(dlength-46);  // remaining 16-bit words
  while( check_length>=21 ){
	  (sgl1->gl1_boards)++;
	  check_length-=21;
  }
 
  if(sgl1->gl1_boards>NUM_GL1_BOARDS)
  {
    sgl1->gl1_boards = NUM_GL1_BOARDS;
  }
  unsigned int* buf = (unsigned int *) findPacketDataStart(packet);
  if (buf == 0) return;

  // unpack the data into a GL1 structure:


  unsigned int idx;
  /* GL1-3 data */

  sgl1->gl3_payload.gl3_header = buf[0]&0xFFFF;
  sgl1->gl3_payload.timestamp.year  = 2000 + (10*((buf[0]>>20)&0xF) + ((buf[0]>>16)&0xF));
  sgl1->gl3_payload.timestamp.month = 10*((buf[1]>>12)&0xF) + ((buf[1]>>8)&0xF);
  sgl1->gl3_payload.timestamp.date  = 10*((buf[1]>>4)&0xF) + (buf[1]&0xF);
  sgl1->gl3_payload.timestamp.day   = (buf[1]>>24)&0x7;
  sgl1->gl3_payload.timestamp.hour  = 10*((buf[1]>>20)&0xF) + ((buf[1]>>16)&0xF);
  sgl1->gl3_payload.timestamp.min   = 10*((buf[2]>>12)&0xF) + ((buf[2]>>8)&0xF) ;
  sgl1->gl3_payload.timestamp.sec   = 10*((buf[2]>>4)&0xF) + (buf[2]&0xF);

  sgl1->gl3_payload.alignment = (buf[2]&0xFFFF0000)>>16;

  sgl1->gl3_payload.beam_crossing_counter[1] = 
    ((buf[3]&0xFFFF)<<16) + ((buf[3]&0xFFFF0000)>>16); 
  sgl1->gl3_payload.beam_crossing_counter[0] = 
    ((buf[4]&0xFFFF)<<16) + ((buf[4]&0xFFFF0000)>>16);

  sgl1->gl3_payload.bunch_crossing_counter = (buf[5]&0xFFFF);
  sgl1->gl3_payload.granule_accept_vector = (buf[5]&0xFFFF0000) + (buf[6]&0xFFFF);
  sgl1->gl3_payload.accept_or_input = (buf[6]&0xFFFF0000)+(buf[7]&0xFFFF);
  sgl1->gl3_payload.gl1_accept_counter = (buf[7]&0xFFFF0000)+(buf[8]&0xFFFF);
  for(unsigned int j=0;j<32;j+=2){
    sgl1->gl3_payload.granule_accept[31-j] = (buf[8+(j/2)]&0xFFFF0000)>>16;
    sgl1->gl3_payload.granule_accept[31-(j+1)] = buf[8+(j/2)+1]&0xFFFF;
  } 
  sgl1->gl3_payload.granule_disables = (buf[24]&0xFFFF0000)+(buf[25]&0xFFFF);
  sgl1->gl3_payload.forced_accepts = (buf[25]&0xFFFF0000)+(buf[26]&0xFFFF);

  /* GL1-2 */

  sgl1->gl2_payload.gl2_header = (buf[26]&0xFFFF0000)>>16;
  sgl1->gl2_payload.partition_accept = ((buf[27]&0xFFFF)<<16) + ((buf[27]&0xFFFF0000)>>16);
  sgl1->gl2_payload.mode_bits = (buf[28]&0xFFFF);

  sgl1->gl2_payload.reduced_bits[1] =
    (buf[28]&0xFFFF0000)+(buf[29]&0xFFFF); 
  sgl1->gl2_payload.reduced_bits[0] =
    (buf[29]&0xFFFF0000)+(buf[30]&0xFFFF); 

  sgl1->gl2_payload.dcm_full_fem_busy = (buf[30]&0xFFFF0000) + (buf[31]&0xFFFF);
  sgl1->gl2_payload.fem_unreliable = (buf[31]&0xFFFF0000) + (buf[32]&0xFFFF);
  sgl1->gl2_payload.granule_busy = (buf[32]&0xFFFF0000) + (buf[33]&0xFFFF);
  sgl1->gl2_payload.part_busy_xbar_out = (buf[33]&0xFFFF0000) + (buf[34]&0xFFFF);
  sgl1->gl2_payload.part_busy_bus = (buf[34]&0xFFFF0000) + (buf[35]&0xFFFF);
 
  /* GL1 data */

  idx = 35;
  for(int i=0; i<sgl1->gl1_boards;i++){

    if( (i%2)==0 ) {   // even numbered boards

      sgl1->gl1_payload[i].gl1_header = (buf[idx]&0xFFFF0000)>>16;
      idx++;

      /* Note the complicated way in which the 20-bit LUT inputs are packed */

      sgl1->gl1_payload[i].lut_input[7] = ((buf[idx]&0xFFFF)<<4)+((buf[idx]&0xF0000000)>>28);
      sgl1->gl1_payload[i].lut_input[6] = ((buf[idx]&0x0FFF0000)>>8)+((buf[idx+1]&0xFF00)>>8);
      sgl1->gl1_payload[i].lut_input[5] = ((buf[idx+1]&0xFF)<<12)+((buf[idx+1]&0xFFF00000)>>20);
      sgl1->gl1_payload[i].lut_input[4] = ((buf[idx+1]&0xF0000)<<16)+(buf[idx+2]&0xFFFF);
      sgl1->gl1_payload[i].lut_input[3] = ((buf[idx+2]&0xFFFF0000)>>12)+((buf[idx+3]&0xF000)>>12);
      sgl1->gl1_payload[i].lut_input[2] = ((buf[idx+3]&0xFFF)<<8)+((buf[idx+3]&0xFF000000)>>24);
      sgl1->gl1_payload[i].lut_input[1] = ((buf[idx+3]&0xFF0000)>>4)+((buf[idx+4]&0xFFF0)>>4);
      sgl1->gl1_payload[i].lut_input[0] = ((buf[idx+4]&0xF)<<16)+((buf[idx+4]&0xFFFF0000)>>16);
      idx+=5;

      sgl1->gl1_payload[i].lut_output = ((buf[idx]&0xFFFF0000)>>16) + ((buf[idx]&0xFFFF)<<16);
      idx++;
      sgl1->gl1_payload[i].trigger_busy = ((buf[idx]&0xFFFF0000)>>16) + ((buf[idx]&0xFFFF)<<16);
      idx++;
      sgl1->gl1_payload[i].live_trig_out = ((buf[idx]&0xFFFF0000)>>16) + ((buf[idx]&0xFFFF)<<16);
      idx++;
      sgl1->gl1_payload[i].scaled_trig_out = ((buf[idx]&0xFFFF0000)>>16) + ((buf[idx]&0xFFFF)<<16);
      idx++;
      sgl1->gl1_payload[i].trig_part_xbar_out = ((buf[idx]&0xFFFF0000)>>16) + ((buf[idx]&0xFFFF)<<16);
      idx++;

    }
    else   // odd numbered boards
    {
 
      sgl1->gl1_payload[i].gl1_header = buf[idx]&0xFFFF;

      /* Note the complicated way in which the 20-bit LUT inputs are packed */

      sgl1->gl1_payload[i].lut_input[7] = ((buf[idx]&0xFFFF0000)>>12)+((buf[idx+1]&0xF000)>>12);
      sgl1->gl1_payload[i].lut_input[6] = ((buf[idx+1]&0x0FFF)<<8) + ((buf[idx+1]&0xFF000000)>>24);
      sgl1->gl1_payload[i].lut_input[5] = ((buf[idx+1]&0xFF0000)>>4) + ((buf[idx+2]&0xFFF0)>>4);
      sgl1->gl1_payload[i].lut_input[4] = ((buf[idx+2]&0xF)<<16) + ((buf[idx+2]&0xFFFF0000)>>16);
      sgl1->gl1_payload[i].lut_input[3] = ((buf[idx+3]&0xFFFF)<<4) + ((buf[idx+3]&0xF0000000)>>28);
      sgl1->gl1_payload[i].lut_input[2] = ((buf[idx+3]&0xFFF0000)>>8) + ((buf[idx+4]&0xFF00)>>8);
      sgl1->gl1_payload[i].lut_input[1] = ((buf[idx+4]&0xFF)<<12) + ((buf[idx+4]&0xFFF00000)>>20);
      sgl1->gl1_payload[i].lut_input[0] = ((buf[idx+4]&0xF0000)) + ((buf[idx+5]&0xFFFF));
      idx+=5;

      sgl1->gl1_payload[i].lut_output = ((buf[idx]&0xFFFF0000)) + ((buf[idx+1]&0xFFFF));
      idx++;
      sgl1->gl1_payload[i].trigger_busy = ((buf[idx]&0xFFFF0000)) + ((buf[idx+1]&0xFFFF));
      idx++;
      sgl1->gl1_payload[i].live_trig_out = ((buf[idx]&0xFFFF0000)) + ((buf[idx+1]&0xFFFF));
      idx++;
      sgl1->gl1_payload[i].scaled_trig_out = ((buf[idx]&0xFFFF0000)) + ((buf[idx+1]&0xFFFF));
      idx++;
      sgl1->gl1_payload[i].trig_part_xbar_out = ((buf[idx]&0xFFFF0000)) + ((buf[idx+1]&0xFFFF));
      idx++;
      
    }

  }

}







