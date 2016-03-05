// -*- c++ -*-
#ifndef __FILEEVENTITERATOR_H__
#define __FILEEVENTITERATOR_H__

#include <stdio.h>

#include "Eventiterator.h"
#include "Event.h"

#include "gzbuffer.h"
#include "oncsBuffer.h"

/**
   The fileEventiterator reads the event data from a data file on disk. 
   It creates and returns pointers to Event objects. At the end of the file 
   it returns 0 when there are no events left.
*/
#ifndef __CINT__
class WINDOWSEXPORT fileEventiterator : public Eventiterator {
#else
class  fileEventiterator : public Eventiterator {
#endif
public:

  virtual ~fileEventiterator();

  /// This simple constructor just needs the file name of the data file.
  fileEventiterator(const char *filename);

  /**
  This constructor gives you a status so you can learn that the creation
  of the fileEventiterator object was successful. If the status is not 0,
  something went wrong and you should delete the object again.
  */
  fileEventiterator(const char *filename, int &status);

  const char * getIdTag() const;

  virtual void identify(std::ostream& os = std::cout) const;

  virtual const char * getCurrentFileName() const;

/**
   this member function returns a pointer to the Event object, or
   NULL if there are no events left.
*/   
  Event *getNextEvent();

  int  setVerbosity(const int v) 
  { 
    verbosity=v;
    return 0; 
  }; 

  int  getVerbosity() const 
  { 
    return verbosity; 
  };


private:
  int open_file(const char *filename);
  int read_next_buffer();

  char *thefilename;
  int fd;
  
  PHDWORD *bp;
  unsigned int allocatedsize;

  int current_index;
  int last_read_status;
  unsigned int buffer_size;
  buffer *bptr;

  int events_so_far;
  int verbosity;
  int _defunct;
};

#endif /* __FILEEVENTITERATOR_H__ */

