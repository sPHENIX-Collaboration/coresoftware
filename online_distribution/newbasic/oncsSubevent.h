#ifndef __ONCSSUBEVENT_H__
#define __ONCSSUBEVENT_H__
// -*- c++ -*-


#include "packet.h"
#include "decoding_routines.h"
#include "oncsStructures.h"
#include <stddef.h>

#include <oncsSubConstants.h>

#ifndef __CINT__
class WINDOWSEXPORT oncsSubevent : public Packet {
#else
class  oncsSubevent : public Packet {
#endif

public:

  oncsSubevent(subevtdata_ptr data);
  virtual ~oncsSubevent();

  // access to information 
  virtual int   getLength() const;
  virtual int   getIdentifier() const;
  virtual int   getHitFormat() const;
  virtual int   getPadding() const;

  virtual int	getDataLength() const;
  virtual int	getStructure() const { return 0;};
  virtual int	getDebugLength() const { return 0;};
  virtual int	getErrorLength() const { return 0;};

  // debugging-type information
  virtual void  identify( OSTREAM& =COUT) const;

  // getting decoded values
  int    iValue(const int);
  int    iValue(const int,const char *);
  int    iValue(const int,const int);

  int    iValue(const int,const int, const char *){return 0;};
  int    iValue(const int channel,const int iy, const int iz) {return 0;};

  int    iValue(const int channel,const int iy, const int iz, const char *what) {return 0;};

  virtual float  rValue(const int);
  virtual float  rValue(const int,const char *);
  virtual float  rValue(const int,const int);

  virtual int    getArraylength(const char *);
  virtual int    fillIntArray (int [], const int, int *,const char *what="");
  virtual int    fillFloatArray (float [], const int, int *,const char *what="");
  virtual int*   getIntArray (int *,const char *what="");
  virtual float* getFloatArray (int *,const char *what="");

  // pointer or data based handling
  virtual int is_pointer_type() const;
  virtual int convert();


protected:

  virtual int *decode(int *) =0;


  subevtdata_ptr SubeventHdr;
  int data1_length;
  int data2_length;
  int data3_length;
  int data4_length;

  int *decoded_data1;
  int *decoded_data2;
  int *decoded_data3;
  int *decoded_data4;
  int is_data_type;
};

// --------------------------------------------------------------------------
#ifndef __CINT__
class WINDOWSEXPORT oncsSubevent_w1 : public oncsSubevent {
#else
class  oncsSubevent_w1 : public oncsSubevent {
#endif
public:
  oncsSubevent_w1(subevtdata_ptr);

  virtual void  dump ( OSTREAM& os = COUT) ;

  virtual void  gdump (const int how=EVT_HEXADECIMAL, OSTREAM& os = COUT) const;

protected:
  inline virtual int *decode (int *) {return 0;};

};

// ----------------------------------------------------

#ifndef __CINT__
class WINDOWSEXPORT oncsSubevent_w2 : public oncsSubevent {
#else
class  oncsSubevent_w2 : public oncsSubevent {
#endif
public:
  oncsSubevent_w2(subevtdata_ptr);

  virtual void  dump ( OSTREAM& os = COUT) ;

  virtual void  gdump (const int how=EVT_HEXADECIMAL, OSTREAM& os = COUT) const;

protected:
  inline virtual int *decode (int *) {return 0;};

};

// ----------------------------------------------------


#ifndef __CINT__
class WINDOWSEXPORT oncsSubevent_w4 : public oncsSubevent {
#else
class  oncsSubevent_w4 : public oncsSubevent {
#endif
public:
  oncsSubevent_w4(subevtdata_ptr);

  virtual void  dump ( OSTREAM& os = COUT) ;

  virtual void  gdump (const int how=EVT_HEXADECIMAL, OSTREAM& os = COUT) const;

protected:
  inline virtual int *decode (int *) {return 0;};

};


#endif /* __ONCSSUBEVENT_H__ */
