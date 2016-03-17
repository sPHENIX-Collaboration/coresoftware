// -*- c++ -*-
#ifndef __ONCSEVENTITERATOR_H__
#define __ONCSEVENTITERATOR_H__

#include "Eventiterator.h"
#include "oncsBuffer.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsEventiterator : public Eventiterator {
#else
class  oncsEventiterator : public Eventiterator {
#endif
public:

  virtual ~oncsEventiterator();
  oncsEventiterator(const char *filename);
  oncsEventiterator(const char *filename, int &status);

  const char * getIdTag() const;
  virtual void identify(std::ostream& os = std::cout) const;


  Event *getNextEvent();

private:
  int read_next_buffer();
  
  char * thefilename;
  int fd;
  int initialbuffer[BUFFERSIZE];
  int *bp;
  int allocatedsize;

  int current_index;
  int last_read_status;
  int buffer_size;
  oncsBuffer *bptr;

};

#endif /* __ONCSEVENTITERATOR_H__ */

