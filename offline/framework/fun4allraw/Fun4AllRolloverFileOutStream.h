// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLRAW_FUN4ALLROLLOVERFILEOUTSTREAM_H
#define FUN4ALLRAW_FUN4ALLROLLOVERFILEOUTSTREAM_H

#include "Fun4AllFileOutStream.h"

#include <iostream>  // for cout, ostream
#include <string>    // for string

class Event;

class Fun4AllRolloverFileOutStream : public Fun4AllFileOutStream
{
 public:
  Fun4AllRolloverFileOutStream(const std::string &frule = "OUTDATA-%010d-%04d.PRDFF",
                               const unsigned int sizeInMB = 0,
                               const int offset = 0,
                               const int increment = 1,
                               const std::string &name = "Fun4AllRolloverFileOutStream");
  virtual ~Fun4AllRolloverFileOutStream() {}
  int WriteEventOut(Event *evt);
  void identify(std::ostream &os = std::cout) const;

 private:
  unsigned long long m_MaxFileFize;
  int m_CurrentSequence;
  int m_Offset;
  int m_Increment;
};

#endif /* __FUN4ALLFILEOUTSTREAM_H__ */
