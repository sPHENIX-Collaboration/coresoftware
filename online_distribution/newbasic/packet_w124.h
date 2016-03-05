#ifndef __PACKET_W124_H__
#define __PACKET_W124_H__

#include "packet_A.h"
#include "packetHeaders.h"
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/**
   Based on the packet\_A definition, we build the
   other classes of packets, which are the base classes 
   of packets with wordsizes 1, 2, and 4.

   Note that this class can be instantiated and will be used if we 
   encounter a packet with an as yet unknown decoding method, or if we
   want to customize the decoding step. It can still handle all the 
   operations such as return raw data or envelope information, but none
   of those which require the data to be decoded.
*/

#ifndef __CINT__
class WINDOWSEXPORT Packet_w1 : public Packet_A {
#else
class  Packet_w1 : public Packet_A {
#endif
public:

  Packet_w1();
  Packet_w1(PACKET_ptr);

  void  dump ( OSTREAM& =COUT ) ;
  void  gdump (const int how = EVT_HEXADECIMAL, OSTREAM& = COUT) const;

protected:
  inline int *decode (int *) {return 0;};

};

// ----------------------------------------------------
/**
   Based on the packet\_A definition, we build the
   other classes of packets, which are the base classes 
   of packets with wordsizes 1, 2, and 4.

   Note that this class can be instantiated and will be used if we 
   encounter a packet with an as yet unknown decoding method, or if we
   want to customize the decoding step. It can still handle all the 
   operations such as return raw data or envelope information, but none
   of those which require the data to be decoded.
*/

#ifndef __CINT__
class WINDOWSEXPORT Packet_w2 : public Packet_A {
#else
class  Packet_w2 : public Packet_A {
#endif
public:
  Packet_w2();
  Packet_w2(PACKET_ptr);

  void  dump ( OSTREAM& ) ;
  void  gdump (const int how = EVT_HEXADECIMAL, OSTREAM& = COUT) const;

protected:
  inline int *decode (int *) {return 0;};
};


// ----------------------------------------------------

/**
   Based on the packet\_A definition, we build the
   other classes of packets, which are the base classes 
   of packets with wordsizes 1, 2, and 4.

   Note that this class can be instantiated and will be used if we 
   encounter a packet with an as yet unknown decoding method, or if we
   want to customize the decoding step. It can still handle all the 
   operations such as return raw data or envelope information, but none
   of those which require the data to be decoded.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_w4 : public Packet_A {
#else
class  Packet_w4 : public Packet_A {
#endif
public:

  Packet_w4();
  Packet_w4(PACKET_ptr);


  void  dump ( OSTREAM& ) ;
  void  gdump (const int how = EVT_HEXADECIMAL, OSTREAM& = COUT) const;

protected:
  inline int *decode (int *) {return 0;};
};


#endif /* __PACKET_W124_H__ */



