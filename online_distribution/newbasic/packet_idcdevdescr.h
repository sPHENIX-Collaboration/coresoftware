#ifndef __PACKET_IDCDEVDESCR_H__
#define __PACKET_IDCDEVDESCR_H__


#include "packet_w124.h"
#include <set>
#include <vector>

/**
   This is the packet which deals with data in IDCDEVDESCR format.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_idcdevdescr : public Packet_w1 {
#else
class  Packet_idcdevdescr : public Packet_w1 {
#endif

public:
  Packet_idcdevdescr(PACKET_ptr);
  ~Packet_idcdevdescr();

  double dValue(const int channel, const char *what);

  void  dump ( OSTREAM& os=COUT) ;

protected:

  int is_decoded;

  virtual int *decode ();
  int *decode(int*) { return decode(); }  // unhide Packet_w1::decode(int*)

  struct namedVector
  {
    char  name[200];
    std::vector<double> values;
  };



  struct valuesetPtrLess:
    public std::binary_function<const namedVector  *, const namedVector *, bool> 
  {
    bool operator() (const namedVector *lhs, const namedVector *rhs) const
    {
      return  ( strcmp ( lhs->name, rhs->name) < 0 );      
    }
  };
  

  typedef std::set<namedVector *, valuesetPtrLess> packetidcdevdescrSet;
  typedef packetidcdevdescrSet::iterator packetidcdevdescrSetiter;

  packetidcdevdescrSet nvset;


};

#endif /* __PACKET_IDCDEVDESCR_H__ */
