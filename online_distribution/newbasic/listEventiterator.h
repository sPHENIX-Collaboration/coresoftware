#ifndef __LISTEVENTITERATOR_H__
#define __LISTEVENTITERATOR_H__

#include <cstdio>

#ifndef __CINT__
#include <fstream>
#endif

#include "fileEventiterator.h"
#include "Event.h"

#define FEMAXFILENAMELENGTH 256

/**
   The listEventiterator reads the event data from a list of data files
   on disk.  It creates and returns pointers to Event objects. At the
   end of the file it returns 0 when there are no events left.  */

#ifndef __CINT__
class WINDOWSEXPORT listEventiterator : public Eventiterator {
#else
class  listEventiterator : public Eventiterator {
#endif
public:

  virtual ~listEventiterator();

  /// This simple constructor just needs the file name of the file 
  /// that contains, line by line, teh filenames.

  listEventiterator(const char *filename);

  /**
  This constructor gives you a status so you can learn that the creation
  of the fileEventiterator object was successful. If the status is not 0,
  something went wrong and you should delete the object again.
  */
  listEventiterator(const char *filename, int &status);

  const char * getIdTag() const;

  virtual void identify(std::ostream& os = std::cout) const;

  virtual const char * getCurrentFileName() const;

/**
   this member function returns a pointer to the Event object, or
   NULL if there are no events left.
*/   
  Event *getNextEvent();


private:

  Eventiterator * getNextIterator();

  char thefilename[FEMAXFILENAMELENGTH];
  char * listfilename;

#ifndef __CINT__
  std::ifstream *liststream;
#endif

  Eventiterator *it;
  int defunct;
  int finished;

};



#endif /* __LISTEVENTITERATOR_H__ */
