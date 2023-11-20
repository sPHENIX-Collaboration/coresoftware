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
                               const unsigned int nEvents = 0,
                               const unsigned int sizeInMB = 0,
                               const int offset = 0,
                               const int increment = 1,
                               const std::string &name = "Fun4AllRolloverFileOutStream");

  virtual ~Fun4AllRolloverFileOutStream()=default;
  int WriteEventOut(Event *evt) override;
  void identify(std::ostream &os = std::cout) const;

 private:
  void open_new_file();
  uint64_t m_MaxFileFize {0};
  unsigned int m_MaxNEvents {0};
  int m_CurrentSequence {0};
  int m_Offset {0};
  int m_Increment {1};
};

#endif /* FUN4ALLRAW_FUN4ALLROLLOVERFILEOUTSTREAM_H */
