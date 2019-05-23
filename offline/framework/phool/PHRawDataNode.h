#ifndef PHOOL_PHRAWDATANODE_H
#define PHOOL_PHRAWDATANODE_H

//  Declaration of class PHRawDataNode
//  Purpose: Node digested by the PHRawOManager
//  Author: Matthias Messer

#include "PHDataNode.h"

#include <Event/phenixTypes.h>

#include <string>

class PHIOManager;

class PHRawDataNode : public PHDataNode<PHDWORD>
{
 public:
  PHRawDataNode(PHDWORD *, const std::string &, const int, const int, const int, const int);
  virtual ~PHRawDataNode();

 public:
  virtual bool write(PHIOManager *, const std::string & = "");

  int getLength() const { return length; }
  int getID() const { return ID; }
  int getWordLength() const { return wordLength; }
  int getHitFormat() const { return hitFormat; }
  void setLength(const int val) { length = val; }
  void setID(const int val) { ID = val; }
  void setWordLength(const int val) { wordLength = val; }
  void setHitFormat(const int val) { hitFormat = val; }

 private:
  PHRawDataNode() = delete;
  int length;
  int ID;
  int wordLength;
  int hitFormat;
};

#endif
