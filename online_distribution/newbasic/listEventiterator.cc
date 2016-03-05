//
// listEventIterator   mlp 4/19/1997
//
// this iterator reads events froma data file. 


#include "listEventiterator.h"
#include <stddef.h>
#include <string.h>

// there are two similar constructors, one with just the
// filename, the other with an additional status value
// which is non-zero on return if anything goes wrong. 


listEventiterator::~listEventiterator()
{
     if (it) delete it;
     if (liststream) delete liststream;
     if ( listfilename ) delete [] listfilename;
}  


listEventiterator::listEventiterator(const char *filename)
{
  listfilename = new char[strlen(filename)+1];
  strcpy (listfilename, filename);

  liststream = new std::ifstream(filename);
  finished = 0;
  defunct = 0;
#ifdef OSF1
  if (liststream) 
#else
  if (!liststream->is_open() ) 
#endif
    {
      defunct = 1;
      finished = 1;
    }
  it = 0;
  strcpy(thefilename, " ");
}


listEventiterator::listEventiterator(const char *filename, int &status)
{
  listfilename = new char[strlen(filename)+1];
  strcpy (listfilename, filename);

  liststream = new std::ifstream(filename);
  finished = 0;
  defunct = 0;
  status = 0;
#ifdef OSF1
  if (liststream) 
#else
  if (!liststream->is_open() ) 
#endif
    {
      status = 1;
      defunct = 1;
      finished = 1;
    }
  it = 0;
  strcpy(thefilename, " ");
}



void listEventiterator::identify (OSTREAM &os) const
{ 
  os << "listEventiterator from list " << listfilename << ", current from " << thefilename << std::endl;

};


const char * listEventiterator::getCurrentFileName() const
{ 
  static char namestr[512];

  strcpy (namestr, thefilename);
  return namestr;
 
};


const char *  listEventiterator::getIdTag () const
{ 
  //  sprintf (idline, " -- listEventiterator reading from %s", thefilename);
  return "listEventiterator";
};



// and, finally, the only non-constructor member function to
// retrieve events from the iterator.

Eventiterator * listEventiterator::getNextIterator()
{

  Eventiterator *fit;
  int status = 0;

  if (finished) return 0;

  if (it) return it; //protect against overwriting our own.

  while ( liststream->getline(thefilename,FEMAXFILENAMELENGTH) )
    {
      fit = new fileEventiterator(thefilename,status);
      if (! status)
	{
	  std::cout << "opened file " << thefilename << std::endl;
	  return fit;
	}
    }
  finished = 1;
  delete liststream;
  liststream = 0;
  return 0;
}
	    

Event * listEventiterator::getNextEvent()
{
  while ( ! finished)   // we still have more file to process
    {
      if (it ==0 )      // if we don't have one open...
	{
	  it = getNextIterator();     // try to open it
	  if (! it)                   // no more there, done
	    {
	      return 0;
	    }
	}
      Event *e = it->getNextEvent();  // get the nexct event
      if (e) return e;                // and if we got it, we're done
      delete it;                      // if that prdf was exhausted, try the next one
      it = 0;
    }
  return 0;
}
