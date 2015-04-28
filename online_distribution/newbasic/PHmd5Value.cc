#include <PHmd5Value.h>

#include <PHmd5Utils.h>
#include <iomanip>

PHmd5Value::PHmd5Value()
{
  int i;
  for (i = 0; i < PHMD5DIGESTLENGTH; i++)
    theDigest[i] = 0;
  isNotValid = 1;  // no value set yet.

}

PHmd5Value::PHmd5Value(const PHmd5Value &other)
{
  if ( other.getMD5(theDigest) ) isNotValid =1;
  else isNotValid = 0;
}

PHmd5Value::PHmd5Value(FILE *fp)
{
  if ( PHmd5Stream(fp, theDigest, &filesize) ) isNotValid =1;
  else isNotValid = 0;
}

PHmd5Value::PHmd5Value(const char *filename)
{
  if ( PHmd5File(filename, theDigest, &filesize) ) isNotValid =1;
  else isNotValid = 0;
}

PHmd5Value::~PHmd5Value()
{}

int PHmd5Value::Status() const
{
  return isNotValid;
}

int PHmd5Value::FileSize() const
{
  return filesize;
}


int PHmd5Value::operator== (const PHmd5Value &other) const
{
  //if either object has a bad status, we declare a non-match
  if (  isNotValid || other.Status() ) return 0;

  unsigned char otherDigest[PHMD5DIGESTLENGTH];
  other.getMD5 (otherDigest);
  int i;
  for (i = 0; i < PHMD5DIGESTLENGTH; i++)
    if (theDigest[i] != otherDigest[i]) return 0;
 
  return 1;
}

int PHmd5Value::getMD5(unsigned char * digest) const
{

  // if we weren't constructed right, don't hand out the
  //digest
  if ( isNotValid ) return 1;

  int i;
  for (i = 0; i < PHMD5DIGESTLENGTH; i++)
    digest[i] = theDigest[i];

  return 0;

}

int PHmd5Value::setMD5(const unsigned char * digest)
{
  
  int i;
  for (i = 0; i < PHMD5DIGESTLENGTH; i++)
    theDigest[i] = digest[i];
  isNotValid = 1;
  return 0;

}


int PHmd5Value::setFileMD5(const char * filename)
{
  if ( PHmd5File(filename, theDigest, &filesize) ) isNotValid =1;
  else isNotValid = 0;
  return isNotValid;
}

OSTREAM & operator<< ( OSTREAM &os, const PHmd5Value &md )
{
  int i;

  //  char v[2];

  if ( !md.Status() )
    {
      for (i = 0; i < PHMD5DIGESTLENGTH; i++)
	{
	  //  sprintf(v, "%02x", md.theDigest[i]);
	  os << std::setw(2) << std::hex << md.theDigest[i] << std::dec;
	}
    }
  else
    {
      for (i = 0; i < PHMD5DIGESTLENGTH; i++)
	os << "00";
    }
    
  return os;
}
