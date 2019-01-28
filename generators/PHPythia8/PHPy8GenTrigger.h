#ifndef __PHPY8GENTRIGGER_H__
#define __PHPY8GENTRIGGER_H__

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

namespace Pythia8 {
  class Pythia;
};

class PHPy8GenTrigger {

 protected:  
  //! constructor
  PHPy8GenTrigger(const std::string &name = "PHPy8GenTrigger");

 public:
  virtual ~PHPy8GenTrigger();

  virtual bool Apply(Pythia8::Pythia *pythia) {
    std::cout << "PHPy8GenTrigger::Apply - in virtual function" << std::endl;
    return false;
  }

  virtual std::string GetName() { return _name; }
  
  std::vector<int> convertToInts(std::string s);

  void Verbosity(int v) { _verbosity = v; }

protected:
  int _verbosity;  
  
private:
  std::string _name;
};

#endif	

