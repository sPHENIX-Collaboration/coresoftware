#ifndef __PACKET_GL1PSUM_H__
#define __PACKET_GL1PSUM_H__

#include <packet_w124.h>

/**
   This is the packet which deals with data in GL1PSUM format.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_gl1psum : public Packet_w4 {
#else
class  Packet_gl1psum : public Packet_w4 {
#endif
public:
  Packet_gl1psum(PACKET_ptr);
  ~Packet_gl1psum();
 virtual int  iValue(const int channel, const char *what);
 virtual int  iValue(const int channel);
 virtual void dump ( OSTREAM& );
 
protected:
  virtual int *decode (int *);
};

#endif /* __PACKET_GL1PSUM_H__ */
