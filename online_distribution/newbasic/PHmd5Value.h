#ifndef __PHMD5VALUE_H__
#define __PHMD5VALUE_H__

/** 
This is a simple class maintaining a md5 checksum value for us. 
It is geared towards the md5 value from files, and accepts 
filenames, a FILE * in its constructor. 

It is made so that it gives the same result as the UNIX md5sum command.

I used and modified the md5_stream routine in the GNU textutils package
(which has a larger scope and doesn't compile easily with non-GNU
compilers).

I used the "MD5 message digest library" written by L. Peter Deutsch
(ghost@aladdin.com) to actually compute the MD5 sums. This package
compiles happily on all platforms, including Sun/Solaris w/ native
compiler.

*/

#include "PHmd5Utils.h"

#define PHMD5DIGESTLENGTH 16
#include "event_io.h"


#define WINDOWSEXPORT


#ifndef __CINT__
class WINDOWSEXPORT  PHmd5Value {
#else
class   PHmd5Value {
#endif
 public:
  /// construct an "empty" object
  PHmd5Value();

  /// the copy constructor
  PHmd5Value(const PHmd5Value &);

  /// give an open file pointer to gte the checksum of that file
  PHmd5Value(FILE *fp);

  /// give a filename to get the MD5 sum of that file
  PHmd5Value(const char *filename);

  // and the destructor
  ~PHmd5Value();

  /// The Status() is 0 if the md5 could be computed alright.
  int Status() const;

  /// overload the == operator to tell if two MD5 sums are equal
  int operator== (const PHmd5Value &) const;

  /// copy the digest to the "digest" parameter
  int getMD5(unsigned char * digest) const;

  /// set a new digest value from the argument
  int setMD5(const unsigned char * digest);

  /// if we change the file, this allows to update the MD5 sum
  int setFileMD5(const char * filename);

  /// this is for convenience -- once we digest the file, we might
    /// as well remember its size to get back at it later.
  int FileSize() const;

  /// define the ostream >> operator
  friend OSTREAM & operator<< ( OSTREAM &, const PHmd5Value &);

 protected:
  unsigned char theDigest[PHMD5DIGESTLENGTH];
  int isNotValid; // here we indicate that an actual value has been set alright (0);
  int filesize;
};



#endif /* __PHMD5VALUE_H__ */
