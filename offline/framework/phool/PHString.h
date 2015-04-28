//  Declaration of class PHString
//
//  Purpose: yet another string class (nice, though)
//
//  Author: Matthias Messer

#ifndef PHSTRING_H
#define PHSTRING_H

#include "phool.h"

#include <iosfwd>
#include <cstring>

PHString operator+ (const PHString&, const PHString&);

template <class T> class PHPointerList;

class PHString { 

public: 
   PHString(); 
   PHString(const PHString&); 
   PHString(const char *); 
   ~PHString(); 

public:
   PHBoolean find(const PHString&) const;
   PHBoolean operator== (const PHString&) const;
   PHBoolean operator!= (const PHString&) const;
   PHString& operator=  (const PHString&);
   PHString& operator+= (const PHString&);
   PHBoolean operator== (const char* s) const { return strcmp(string, s)==0; }
   friend PHString operator+ (const PHString&, const PHString&);
   size_t split(PHPointerList<PHString>&, const char*) const;
   size_t length() const { return string ? strlen(string) : 0; }
   char* getString() const { return string; }
   void setString(char * s);
   
private:
   void copy_(const char*);
   
private: 
   char *string;
}; 

std::ostream & operator << (std::ostream &, const PHString &);

PHString join(PHPointerList<PHString>&, const PHString&);

#endif /* PHSTRING_H */ 


