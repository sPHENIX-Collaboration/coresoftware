#ifndef PHIOMANAGER_H__
#define PHIOMANAGER_H__

//  Declaration of class PHIOManager
//  Purpose: Abstract base class for file IO
//  Author: Matthias Messer

#include <string>

class PHCompositeNode;

class PHIOManager
{
 public:
  virtual ~PHIOManager() {}
 public:
  std::string getFilename() const { return filename; }
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
  PHIOManager();
  std::string filename;
  size_t eventNumber;
};

#endif /* PHIOMANAGER_H__ */
