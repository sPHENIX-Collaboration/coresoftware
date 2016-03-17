// -*- c++ -*-
#ifndef __EVENT_H__
#define __EVENT_H__

//#include <stddef.h>

#include "packet.h"
#include <ctime>

// virtual base class for an "event"

const char *get_evt_mnemonic(const int);

/** 
This is the abstract Event class.

It defines all interfaces to deal with Event data. It is a pure
virtual class which does not define any implementation details (it has
only member function defintions, no data members). It is meant to be
subclassed by the classes which define the actual implementation.
*/

//class PHLeanTimeStamp;

#ifndef __CINT__
class WINDOWSEXPORT Event {
#else
class  Event {
#endif

public:
  inline virtual ~Event(){};

  // **** info & debug utils ******
  /** getEvtLength() returns the length of the event raw data
   in longwords (32bit). If you want to copy the event somewhere,
   you can find out whether or not you have enough space left.
  */
  virtual unsigned int getEvtLength() =0;

  /**
      getEvtType() returns the type of the event. Most of the
   events have the type DATAEVENT, but there are also
   BEGRUNEVENT, SPILLONEVENT, SPILLOFFEVENT, ENDRUNEVENT,
   and so on.
  */
  virtual int getEvtType() =0;

  /**
   getEvtSequence() returns the number of the event in a 
   particular run. This number is a property of the event. Its run number and the 
   sequence number uniquely indentify an event. It has nothing to do with the position
   of the event in any given data file.
  */
  virtual int getEvtSequence()  =0;
  
  /**
   getRunNumber() returns the number of the run to which this
   event belongs.
  */
  virtual int getRunNumber() =0;

  /**
   identify will write a short identification message to 
   standard output or to the ostream provided. 
	*/


  virtual void identify(std::ostream& os = std::cout) const = 0;

  /**
   listFrame will write a short descriptor of the frame header to the output stream 
   provided. If id == 0, it will list all frame headers. Otherwise it will list
   the header of the frame in which the packet with id resides. 
  */



  virtual void listFrame(  const int id=0, std::ostream& os=std::cout) const {};

  virtual void listHistory(  const int id=0, std::ostream& os=std::cout) const {};
  virtual void listError(  const int id=0, std::ostream& os=std::cout) const {};

  /**
   getFrameEntry will return a particular entry from the frame header which contains
   the packet with the specidied id. If id == 0, the first frame will be selected.

   keywords to return are:
   
   FRAMELENGTH
   FRAMEMARK
   FRAMEHDRVERSION
   FRAMEHDRLENGTH
   FRAMESTATUS
   FRAMESOURCEID
   FRAMEDATATYPE
   FRAMETYPE
   FRAMEALIGNLENGTH
   FRAMEALIGNMENTWORD[index]

   The last  FRAMEALIGNMENTWORD is an array, you need to call 
   getFrameEntry(id,"FRAMEALIGNMENTWORD",1) to get the index 1, and so on. 
   Asking for a word beyond the number of such words will return 0.
  */

  virtual unsigned int getFrameEntry(  const char *what, const int id=0, const int index =0) const { return 0; };


  // ***** packet handling *****
  /**
   getPacket creates and returns a pointer to the
   packet object with the specified packet id (think of it 
   as a house number). 
  */
  virtual Packet* getPacket(const int)=0;

  /**
   This interface allows to override the hitformat of the packet.
   For debugging purposes mostly.
    
  */
  virtual Packet* getPacket(const int, const int hitFormat)=0;

  /**
    getPacketList returns a simple-minded array of pointers
    to the packets in the given event. It returns the number 
    of packets returned in the array. The second parameter tells
    the packet how long the array is.
  */
  virtual int getPacketList(Packet*[], const int length) =0;

  /**
   existPacket returns 1 if such a packet exists in the
   event. 0 else.
  */
  virtual int existPacket (const int packetid)=0;

  // **** event copying *****
  /**
   the event's raw data are copied to destination. length is the
   size of destination, and nw is the number of words
   actually copied.

   if what = "RAW" then we just copy the payload data without the event header.
  */
  virtual int Copy ( int *destination, const unsigned int length, int *nw, const char * what ="" )=0;


  /**
   getErrorCode() returns a non-zero error code if something in the event
   structure is wrong. Test for 0 to see that the event passes some 
   consistency checks.
  */
  virtual int getErrorCode() {return 0;};


  /**
   getTagWord() returns the bit pattern that the event builder may have put into
   the header for fast event selection purposes. This is mainly used by the
   ET system to distribute events based on its properties.
  */
  virtual unsigned int getTagWord( const int i=0) const { return 0;};

  virtual int is_pointer_type() const =0;

  // **** convert event from pointer- to data-based *****
  /**
   converting the Event object means that the actual data it manages 
   are copied to an internal location and are thus safe from being overwritten. 
   For efficiency,  most Event objects initially maintain a pointer to the 
   external raw data only. They are referred to as being pointer-based. 
   If an Event object has already been converted, a second convert operation 
   has no effect.  
  */
  virtual int convert ()=0;

  virtual int getDate() = 0;
  virtual time_t getTime() const = 0;
};

#endif /* __EVENT_H__ */
