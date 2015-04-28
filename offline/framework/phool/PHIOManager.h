#ifndef __PHIOMANAGER_H__
#define __PHIOMANAGER_H__

//  Declaration of class PHIOManager
//  Purpose: Abstract base class for file IO
//  Author: Matthias Messer

#include "phool.h"
#include "PHString.h"

class PHCompositeNode;

class PHIOManager { 
public: 
   PHIOManager();
   virtual ~PHIOManager(){} 

public:
   PHString getFilename() const;
   size_t getEventNumber() const { return eventNumber; }
   void setEventNumber(const size_t evno) { eventNumber = evno; return;}
   virtual void closeFile() = 0;
   virtual PHBoolean write(PHCompositeNode *) = 0;
   virtual void print() const = 0;
   
protected: 
   PHString        filename;
   size_t          eventNumber;
}; 

#endif /* __PHIOMANAGER_H__ */
