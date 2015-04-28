#ifndef __PHRAWDATANODE_H__
#define __PHRAWDATANODE_H__

//  Declaration of class PHRawDataNode
//  Purpose: Node digested by the PHRawOManager
//  Author: Matthias Messer

#include "PHDataNode.h"

#include <Event/phenixTypes.h>

class PHRawDataNode : public PHDataNode<PHDWORD> 
{ 

public: 
   PHRawDataNode();
   PHRawDataNode(PHDWORD *, const PHString&, const int, const int, const int, const int);
   virtual ~PHRawDataNode();

public:
   virtual PHBoolean write(PHIOManager *, const PHString& = "");

   int getLength()     const { return length; }
   int getID()         const { return ID; }
   int getWordLength() const { return wordLength; }
   int getHitFormat()  const { return hitFormat; }

   void setLength(int val)     { length = val; }
   void setID(int val)         { ID = val; }
   void setWordLength(int val) { wordLength = val; }
   void setHitFormat(int val)  { hitFormat = val; }

private: 
   int length;
   int ID;
   int wordLength;
   int hitFormat;
}; 

#endif /* __PHRAWDATANODE_H__ */
