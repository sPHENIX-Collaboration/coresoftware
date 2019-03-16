#ifndef FUN4ALL_FUN4ALLFILEOUTSTREAM_H
#define FUN4ALL_FUN4ALLFILEOUTSTREAM_H

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
  oBuffer *m_ob;
  int m_iSeq;
  PHDWORD m_xb[LENGTH];
  int m_OutFileDesc;
  unsigned long long m_BytesWritten;
  unsigned long long m_MaxSize;
};

#endif
