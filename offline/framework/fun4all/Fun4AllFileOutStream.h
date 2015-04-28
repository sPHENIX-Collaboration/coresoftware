#ifndef __FUN4ALLFILEOUTSTREAM_H__
#define __FUN4ALLFILEOUTSTREAM_H__

#include "Fun4AllEventOutStream.h"

#include <Event/phenixTypes.h>

#include <string>
#include <iostream>

class Event;
class oBuffer;

static const unsigned int LENGTH = (4*1024*1024);

class Fun4AllFileOutStream: public Fun4AllEventOutStream
{
 public:
  Fun4AllFileOutStream(const std::string &frule = "OUTDATA-%010d-%04d.PRDFF", const std::string &name = "FILEOUTSTREAM");
  virtual ~Fun4AllFileOutStream();
  int WriteEventOut(Event *evt);
  int CloseOutStream();
  void identify(std::ostream &os = std::cout) const;

 protected:
  std::string filerule;
  oBuffer *ob;
  int iseq;
  PHDWORD xb[LENGTH];
  int outfile_desc;
  unsigned long long byteswritten;
  unsigned long long MAXSIZE;
};

#endif /* __FUN4ALLFILEOUTSTREAM_H__ */
