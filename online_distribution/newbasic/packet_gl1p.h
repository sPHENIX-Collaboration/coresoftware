#ifndef __PACKET_GL1P_H__
#define __PACKET_GL1P_H__

/*  
    GL1P Accepted Event Data - 

    This header file defines a standard structure that can be 
    used to map GL1P accepted event data into something a 
    bit more "useful".
*/

/* Each GL1P will have an output packet in the data stream */

#define GL1P_HEADER 0
#define GL1P_EVNUMBER 1
#define GL1P_MODEBIT 2
#define GL1P_CLOCK 3
#define GL1P_SCALER 4
/* A structure for the GL1P accepted event data */

typedef struct{
  unsigned char  header; //GL1P header word. Depend on GL1P number. First is 7a, second - 8a, ...
  unsigned char  ev_number; //1-byte event number counter. 
  unsigned short  modebit; //GL1 mode bit 
  unsigned char  clock[4]; //Beam Crossing Counters for all 4 scalers
  unsigned int   scaler[4]; //Values of scalers A, B, C, D.
} GL1P_DATA;


#include <packet_w124.h>

/**
   This is the packet which deals with data in GL1 format.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_gl1p : public Packet_w4 {
#else
class  Packet_gl1p : public Packet_w4 {
#endif
public:
  Packet_gl1p(PACKET_ptr);
  ~Packet_gl1p();
 virtual int  iValue(const int channel, const char *what);
 virtual int  iValue(const int channel, const int what);
 virtual void dump ( OSTREAM& );
  virtual int    fillIntArray (int destination[],    // the data go here 
			       const int length,      // space we have in destination
			       int * nw,              // words actually used
			       const char * what=""); // type of data (see above)


protected:
  virtual void demangle ();
  virtual int *decode (int *);
  GL1P_DATA* sgl1p;
  int nGL1Pboards;
};

#endif /* __PACKET_GL1P_H__ */
