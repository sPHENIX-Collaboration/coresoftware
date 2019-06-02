// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_EVENTHEADER_H
#define FFAOBJECTS_EVENTHEADER_H

#include <phool/PHObject.h>

#include <ctime>
#include <iostream>          // for cout, ostream

//! base class for EventHeaders
class EventHeader: public PHObject
{
 public:

  /// dtor
  virtual ~EventHeader() {}

  /// Clear Event
  virtual void Reset();

  /** identify Function from PHObject
      @param os Output Stream 
   */
  virtual void identify(std::ostream& os = std::cout) const;

  /// isValid returns non zero if object contains valid data
  virtual int isValid() const;

  /// get Run Number
  virtual int get_RunNumber() const {return -9999;}
  /// set Run Number
  virtual void set_RunNumber(const int run) {return;}

  /// get Event Number
  virtual int get_EvtSequence() const {return -9999;}
  /// set Event Number
  virtual void set_EvtSequence(const int /*ival*/) {return;}

  /// get Event Type (Data,rejected,EOR,BOR,...)
  virtual int get_EvtType() const {return -9999;}
  /// set Event Type (Data,rejected,EOR,BOR,...)
  virtual void set_EvtType(const int /*ival*/) {return;}

  /// get ATP TimeStamp (unix time, convert with ctime()
  virtual time_t get_TimeStamp() const {return 0;}
  /// set TimeStamp
  virtual void set_TimeStamp(const time_t /*evttime*/) {return;}

 private: // prevent doc++ from showing ClassDef
  ClassDef(EventHeader,1)

};

#endif



