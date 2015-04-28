#include "oncsSubevent.h"
#include <string.h>


const char *oncs_get_mnemonic (const int structure,const int format);
const char *get_type_mnemonic (const int id);
const char *get_evt_mnemonic(const int id);


oncsSubevent::oncsSubevent(subevtdata_ptr sevt_ptr)
{
  SubeventHdr = sevt_ptr;

  data1_length = 0;
  data2_length = 0;
  data3_length = 0;
  data4_length = 0;

  decoded_data1 = NULL;
  decoded_data2 = NULL;
  decoded_data3 = NULL;
  decoded_data4 = NULL;

  is_data_type = 0;  //assume "p" type first.

}

oncsSubevent::~oncsSubevent()
{
  if (decoded_data1 != NULL) delete [] decoded_data1;
  if (decoded_data2 != NULL) delete [] decoded_data2;
  if (decoded_data3 != NULL) delete [] decoded_data3;
  if (decoded_data4 != NULL) delete [] decoded_data4;

  if (is_data_type) delete [] (int *) SubeventHdr;
}



// ----------------------------------------------

int oncsSubevent::getLength() const
{
  return SubeventHdr->sub_length;
}

// ----------------------------------------------

int oncsSubevent::getDataLength() const
{
  return SubeventHdr->sub_length -4;
}

// ----------------------------------------------

int oncsSubevent::getIdentifier() const
{
  return SubeventHdr->sub_id;
}

// ----------------------------------------------
//int oncsSubevent::get_type()
//{
//  return SubeventHdr->sub_type;
//}

// -----------------------------------------------
int oncsSubevent::getHitFormat() const
{
  return SubeventHdr->sub_decoding;
}

// -----------------------------------------------
int oncsSubevent::getPadding() const
{
  return SubeventHdr->sub_padding;
}



void  oncsSubevent::identify( OSTREAM& out ) const
{
  out << std::dec
      << "Packet " << SETW(5) <<  getIdentifier() 
      << " " << SETW(5) << getLength() 
      << " -1"  << " (ONCS Packet)";

  out << SETW(3) << getHitFormat() 
      << " (" << oncs_get_mnemonic( getStructure(), getHitFormat()) << ")" << std::endl;
}


int oncsSubevent::is_pointer_type() const
{
  if (is_data_type) return 0;
  else return 1;
}

int oncsSubevent::convert()
{
  if (is_data_type) return -1;

  int *tmp;

  tmp = new int[getLength()];
  int *from=  (int *) SubeventHdr;
  int *to=tmp;
  for (int k=0; k< getLength(); k++) 
    { 
      *to++ = *from++;
    }
  SubeventHdr = ( subevtdata_ptr )tmp;
  is_data_type = 1;
  return 0;

}


// ------------------------------------------------------

int   oncsSubevent::iValue(const int ich)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
        return 0;
    }

  // see if our array is long enough
  if (ich > data1_length) return 0;

  return decoded_data1[ich];
}

// ------------------------------------------------------

int   oncsSubevent::iValue(const int ich, const char *what)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
        return 0;
    }

  // see if our array is long enough
  if (ich > data1_length) return 0;

  return decoded_data1[ich];
}

// ------------------------------------------------------

int   oncsSubevent::iValue(const int ich, const int iy)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
        return 0;
    }

  // see if our array is long enough
  if (ich > data1_length) return 0;

  return decoded_data1[ich];
}

// ------------------------------------------------------

float oncsSubevent::rValue(const int ich)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
        return 0;
    }

  // see if our array is long enough
  if (ich > data1_length) return 0;

  return float(decoded_data1[ich]);
}

// ------------------------------------------------------

float oncsSubevent::rValue(const int ich, const char *what)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
        return 0;
    }

  // see if our array is long enough
  if (ich > data1_length) return 0;

  return float(decoded_data1[ich]);
}

// ------------------------------------------------------

float oncsSubevent::rValue(const int ich, const int iy)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
        return 0;
    }

  // see if our array is long enough
  if (ich > data1_length) return 0;

  return float(decoded_data1[ich]);
}

// ----------------------------------------------

int   
oncsSubevent::getArraylength(const char *what)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return -1;
    }

  return data1_length;
}

// ----------------------------------------------

int
oncsSubevent::fillIntArray (int iarr[],
			   const int nlen, int *nwout,
			   const char *what)
{

  if (strcmp(what,"") ==0)
    {
      // now let's derefence the proxy array. If we didn't decode
      // the data until now, we do it now
      if (decoded_data1 == NULL )
	{
	  if ( (decoded_data1 = decode(&data1_length))==NULL)
	    {
	      *nwout=0;
	      return -1;
	    }
	}

      // see if our array is long enough
      if (nlen < data1_length) return -2;

      // and copy the data to the output array
      int *from = decoded_data1;
      for (int i=0; i<data1_length; i++) *iarr++ = *from++;

      // tell how much we copied
      *nwout = data1_length;
    }
  
  if (strcmp(what,"RAW") ==0)
    {
      if (nlen < getLength() ) return -2;
      int *from = (int *) SubeventHdr;
      for (int i=0; i<getLength(); i++) *iarr++ = *from++;

      *nwout = getLength();
    }
  return 0;
}

// ------------------------------------------------------

int
oncsSubevent::fillFloatArray (float rarr[],
			     const int nlen, int *nwout,
			     const char *what)
{

  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	{
	  *nwout=0;
	  return -1;
	}
    }

  // see if our array is long enough
  if (nlen < data1_length) return -2;

  // and copy the data to the output array
  int *from = decoded_data1;
  for (int i=0; i<data1_length; i++) *rarr++ = float (*from++);

  // tell how much we copied
  *nwout = data1_length;
  return 0;
}

// ------------------------------------------------------
int*   
oncsSubevent::getIntArray (int *nwout, const char *what)
{

  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	{
	  *nwout=0;
	  return NULL;
	}
    }

  int *temp=new int[data1_length];
  int is = fillIntArray(temp, data1_length, nwout);
  if (is)
    {
      *nwout=0;
      return NULL;
    }
  return temp;
}

// ------------------------------------------------------
float*   oncsSubevent::getFloatArray (int *nwout, const char *what)
{
  // now let's derefence the proxy array. If we didn't decode
  // the data until now, we do it now
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	{
	  *nwout=0;
	  return NULL;
	}
    }

  float *temp=new float[data1_length];
  int is = fillFloatArray(temp, data1_length, nwout);
  if (is)
    {
      *nwout=0;
      return NULL;
    }
  return temp;
}

// ----------------------------------------------



// now the member functions common to all
// the individial p1, p2, and p4 classes.

// ------------Constructors first -----------------

oncsSubevent_w1::oncsSubevent_w1(subevtdata_ptr sevt_ptr)
  : oncsSubevent(sevt_ptr)
{}

oncsSubevent_w2::oncsSubevent_w2(subevtdata_ptr sevt_ptr)
  : oncsSubevent(sevt_ptr)
{}

oncsSubevent_w4::oncsSubevent_w4(subevtdata_ptr sevt_ptr)
  : oncsSubevent(sevt_ptr)
{}

// ---- and the dump routines, which just call the 
// ---- generic dump routines, but may be overloaded
// ---- by the subclasses.

void oncsSubevent_w1::dump(OSTREAM& out)  { gdump(2,out);}
void oncsSubevent_w2::dump(OSTREAM& out)  { gdump(2,out);}
void oncsSubevent_w4::dump(OSTREAM& out)  { gdump(2,out);}

void oncsSubevent_w4::gdump(const int i, OSTREAM& out) const
{

  int *SubeventData = &SubeventHdr->data;
  int j,l;
  identify(out);
  
  switch (i)
    {
    case (EVT_HEXADECIMAL):
      j = 0;
      while (1)
	{
	  out << std::endl << SETW(5) << j << " |  ";
	  for (l=0;l<4;l++)
	    {
	      out << std::hex << SETW(8) << SubeventData[j++] << " ";
	      if (j>=SubeventHdr->sub_length-SEVTHEADERLENGTH) break;
	    }
	  if (j>=SubeventHdr->sub_length-SEVTHEADERLENGTH) break;
	}
      break;

    case (EVT_DECIMAL):
      j = 0;
      while (1)
	{
	  out << std::dec << std::endl << SETW(5) << j << " |  ";

	  for (l=0;l<6;l++)
	    {
	      out << SETW(10) << SubeventData[j++] << " ";
	      if (j>=SubeventHdr->sub_length-SEVTHEADERLENGTH) break;
	    }
	  if (j>=SubeventHdr->sub_length-SEVTHEADERLENGTH) break;
	}
      break;

    default: 
      break;
    }
  out << std::endl;

}

// ---------------------------------------------------------------------

void oncsSubevent_w2::gdump(const int i, OSTREAM& out) const
{
  short *SubeventData = (short *) &SubeventHdr->data;
  int j,l;
  identify(out);

  switch (i)
    {
    case (EVT_HEXADECIMAL):
      j = 0;
      while (1)
	{
	  out << std::dec<< std::endl << SETW(5) << j << " |  ";
	  for (l=0;l<8;l++)
	    {
	      out << std::hex << SETW(4) << SubeventData[j++] << " ";
	      if (j>=2*(SubeventHdr->sub_length-SEVTHEADERLENGTH) ) break;
	    }
	  if (j>=2*(SubeventHdr->sub_length-SEVTHEADERLENGTH) ) break;
	}
      break;

    case (EVT_DECIMAL):
      j = 0;
      while (1)
	{
	  out << std::dec << std::endl << SETW(5) << j << " |  ";
	  for (l=0;l<8;l++)
	    {
	      out << std::dec << SETW(6) << SubeventData[j++] << " ";
	      if (j>=2*(SubeventHdr->sub_length-SEVTHEADERLENGTH) ) break;
	    }
	  if (j>=2*(SubeventHdr->sub_length-SEVTHEADERLENGTH) ) break;
	}
      break;

    default: break;
    }
  out << std::endl;
}

// ---------------------------------------------------------------------

void oncsSubevent_w1::gdump(const int i, OSTREAM& out) const
{
  char *SubeventData = (char *) &SubeventHdr->data;
  int j,l;
  char cstring[20];
  char *c;
  identify(out);

  j = 0;
  switch (i)
    {
    case (EVT_HEXADECIMAL):
      while (1)
	{
	  c = cstring;
	  out << std::dec << std::endl << SETW(5) << j << " |  ";
	  for (l=0;l<16;l++)
	    {
	      if (j < 4*(SubeventHdr->sub_length-SEVTHEADERLENGTH) ) 
		{
		  *c++ = SubeventData[j];
		  out << std::hex << SETW(2) << SubeventData[j++] << " ";
		}
	      else
		{
		  *c++ = 0x20;
		  out << "   ";
		}
	    }
	  *c = 0;
	  out << "  | " << cstring;
	  if (j >= 4*(SubeventHdr->sub_length-SEVTHEADERLENGTH) ) break;
	}
      break;

    case (EVT_DECIMAL):
      while (1)
	{
	  c = cstring;
	  out << std::dec << std::endl << SETW(5) << j << " |  ";
	  for (l=0;l<12;l++)
	    {
	      if (j < 4*(SubeventHdr->sub_length-SEVTHEADERLENGTH) ) 
		{
		  *c++ = SubeventData[j];
		  out << std::hex << SETW(4) << SubeventData[j++] << " ";
		}
	      else
		{
		  *c++ = 0x20;
		  out << "   ";
		}
	    }
	  *c = 0;
	  out << "  | " << cstring;
	  if (j >= 4*(SubeventHdr->sub_length-SEVTHEADERLENGTH) ) break;
	}
      break;

    default: break;
    }
  out << std::endl;
}

