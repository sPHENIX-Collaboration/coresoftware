#ifndef __PACKET_GL1_H__
#define __PACKET_GL1_H__

/*  
    GL1 Accepted Event Data - 

    This header file defines a standard structure that can be 
    used to map GL1 accepted event data into something a 
    bit more "useful".

    The documentation for the GL1 output data packets can be 
    found at:

    http://www.physics.iastate.edu/~npexp/papers.html
    (the article you want is "GL1 Update: Monitor and L1DD Data")

  HISTORY:

    6/08/99 J.Lajoie - first version
    6/09/99 J.Lajoie - updated, added multiple GL1, comments
    7/12/99 J.Lajoie - updated, structure made more reasonable to match data
    9/12/99 J.Lajoie - updated, removed zero words for GL1-1 and GL1-2

*/

/* Each GL1 will have an output packet in the data stream */

#define NUM_GL1_BOARDS  4

#define HEADER3         0
#define YEAR            1
#define MONTH           2
#define DATE            3
#define DAY             4
#define HOUR            5
#define MIN             6
#define SECGL1          7 // previous SEC collided with /usr/include/sys/time.h
#define ALIGNMENT       8
#define CROSSCTR        9
#define BEAMCTR0       10
#define BEAMCTR1       11
#define GACCEPT        12
#define ACPTORINP      13
#define ACPTCTR        14
#define GRANCTR        15
#define GDISABLE       16
#define FACCEPT        17
#define HEADER2        18
#define PACCEPT        19
#define MODEBITS       20
#define RBITS0         21
#define RBITS1         22
#define DCMFULL        23
#define FEMUNREL       24
#define GBUSY          25
#define PXBAR          26
#define PBUSY          27
#define HEADER1        28
#define LUTINPUT       29
#define RAWTRIG        30
#define TRIGBUSY       31
#define LIVETRIG       32
#define SCALEDTRIG     33
#define TRIGPARXBAR    34

/* A structure for the GL1-1 accepted event data */

typedef struct{

  /* GL1 data */
  unsigned short int    gl1_header;                 /* 0x1aXX, XX=lower 8 bits of accept counter */
  unsigned       int    lut_input[8];               /* packed LUT input data */
  unsigned       int    lut_output;                 /* packed LUT output data (raw triggers) */
  unsigned       int    trigger_busy;               /* bit-coded trigger busies */
  unsigned       int    live_trig_out;              /* bit-coded live trigger output */
  unsigned       int    scaled_trig_out;            /* bit-coded scaled trigger output */
  unsigned       int    trig_part_xbar_out;         /* bit coded output of trig->part XBAR */ 

} GL1_1_DATA;

/* A structure for the GL1-2 accepted event data */

typedef struct{

  /* GL2 data */
  unsigned short int    gl2_header;                 /* 0x2aXX, XX=lower 8 bits of accept counter */
  unsigned       int    partition_accept;           /* bit-coded partition accept vector */
  unsigned short int    mode_bits;                  /* GL1 mode bits for this event */
  unsigned       int    reduced_bits[2];            /* 64 bit reduced bit input for this event */  
  unsigned       int    dcm_full_fem_busy;          /* bit-coded dcm busy/FEM full vector */ 
  unsigned       int    fem_unreliable;             /* bit-coded FEM unreliable vector */
  unsigned       int    granule_busy;               /* bit-coded granule busy vector */
  unsigned       int    part_busy_xbar_out;         /* output of partition busy XBAR */
  unsigned       int    part_busy_bus;              /* partition busy bus for this event */

} GL1_2_DATA;

typedef struct{

/* Time Stamp Structure */
  unsigned short year; 
  unsigned short month; 
  unsigned short date; 
  unsigned short day; 
  unsigned short hour; 
  unsigned short min; 
  unsigned short sec;

}  GL1_TIME_STAMP ; 
/* A structure for the GL1-3 accepted event data */

typedef struct{

  /* GL3 data */
  unsigned short int   gl3_header;                  /* >0x3aXX, XX=lower 8 bits of accept counter */
  GL1_TIME_STAMP       timestamp;                   /* BCD coded timestamp */
  unsigned short int   alignment;                   /* system alignment bits */ 
  unsigned       int   beam_crossing_counter[2];    /* 64-bit beam crossing counter */
  unsigned short int   bunch_crossing_counter;      /* bunch crossing counter (reset by fiducial) */
  unsigned       int   granule_accept_vector;       /* bit-coded granule accept vector */
  unsigned       int   accept_or_input;             /* bit-coded input to accept OR */
  unsigned       int   gl1_accept_counter;          /* 32-bit GL1 accept counter */ 
  unsigned short int   granule_accept[32];          /* 16-bit granule accept counters */
  unsigned       int   granule_disables;            /* bit-coded granule disables */
  unsigned       int   forced_accepts;              /* bit-coded forced accepts */

} GL1_3_DATA;

typedef struct{

  GL1_3_DATA   gl3_payload;
  GL1_2_DATA   gl2_payload;
  int          gl1_boards;
  GL1_1_DATA   gl1_payload[NUM_GL1_BOARDS];

} GL1_EVENT_DATA;


#include <packet_w124.h>

/**
   This is the packet which deals with data in GL1 format.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_gl1 : public Packet_w4 {
#else
class  Packet_gl1 : public Packet_w4 {
#endif
public:
  Packet_gl1(PACKET_ptr);
  ~Packet_gl1();
 virtual int  iValue(const int channel, const char *what);
 virtual int  iValue(const int channel, const int what);
 virtual void dump ( OSTREAM& );

protected:
  virtual void demangle ();
  virtual int *decode (int *);
  GL1_EVENT_DATA* sgl1;
};

#endif /* __PACKET_GL1_H__ */

















