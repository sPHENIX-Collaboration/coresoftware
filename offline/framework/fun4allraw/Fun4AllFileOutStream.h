// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLRAW_FUN4ALLFILEOUTSTREAM_H
#define FUN4ALLRAW_FUN4ALLFILEOUTSTREAM_H

#include "Fun4AllEventOutStream.h"

#include <Event/phenixTypes.h>

#include <iostream>
#include <string>

class Event;
class oBuffer;

class Fun4AllFileOutStream : public Fun4AllEventOutStream
{
 public:
  static const unsigned int LENGTH = (4 * 1024 * 1024);
  Fun4AllFileOutStream(const std::string &frule = "OUTDATA-%010d-%04d.PRDFF", const std::string &name = "FILEOUTSTREAM");
  virtual ~Fun4AllFileOutStream();
  int WriteEventOut(Event *evt);
  int CloseOutStream();
  void identify(std::ostream &os = std::cout) const;
  oBuffer *GetoBuffer() { return m_ob; }
  void SetoBuffer(oBuffer *bf) { m_ob = bf; }
  unsigned long long MaxSize() const { return m_MaxSize; }
  void DeleteoBuffer();
  std::string FileRule() const { return m_FileRule; }
  int iSeq() const { return m_iSeq; }
  void iSeq(const int i) { m_iSeq = i; }
  unsigned long long BytesWritten() const { return m_BytesWritten; }
  void BytesWritten(const unsigned long long i) { m_BytesWritten = i; }
  int OutFileDescriptor() const { return m_OutFileDesc; }
  void OutFileDescriptor(const int i) { m_OutFileDesc = i; }
  PHDWORD *xb() { return m_xb; }

 private:
  std::string m_FileRule;
  oBuffer *m_ob = nullptr;
  int m_iSeq = 0;
  PHDWORD m_xb[LENGTH];
  int m_OutFileDesc = -1;
  unsigned long long m_BytesWritten = 0;
  unsigned long long m_MaxSize = 100000000000LL;  // 100GB
};

#endif
