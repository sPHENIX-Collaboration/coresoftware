#ifndef PHSARTRE_PHSARTREGENTRIGGER_H
#define PHSARTRE_PHSARTREGENTRIGGER_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

class Event;

class PHSartreGenTrigger {

 protected:  
  //! constructor
  PHSartreGenTrigger(const std::string &name = "PHSartreGenTrigger");

 public:
  virtual ~PHSartreGenTrigger();

  #ifndef __CINT__
  virtual bool Apply(Event *event) {
    std::cout << "PHSartreGenTrigger::Apply - in virtual function" << std::endl;
    return false;
  }
  #endif

  virtual std::string GetName() { return _name; }
  
  std::vector<int> convertToInts(std::string s);

  void Verbosity(int v) { _verbosity = v; }

protected:
  int _verbosity;  
  
private:
  std::string _name;
};

#endif	

