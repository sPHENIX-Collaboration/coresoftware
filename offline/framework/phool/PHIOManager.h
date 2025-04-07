#ifndef PHOOL_PHIOMANAGER_H
#define PHOOL_PHIOMANAGER_H

//  Declaration of class PHIOManager
//  Purpose: Abstract base class for file IO
//  Author: Matthias Messer

#include <cstddef>
#include <string>

class PHCompositeNode;

class PHIOManager
{
 public:
  virtual ~PHIOManager() {}

 public:
  const std::string &getFilename() const { return filename; }
  size_t getEventNumber() const { return eventNumber; }
  void setEventNumber(const size_t evno)
  {
    eventNumber = evno;
    return;
  }
  virtual void closeFile() = 0;
  virtual bool write(PHCompositeNode *) = 0;
  virtual void print() const = 0;

 protected:
  PHIOManager() = default;
  size_t eventNumber{0};
  std::string filename;
};

#endif
