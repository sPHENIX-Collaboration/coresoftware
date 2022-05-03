#ifndef PHPYTHIA8_PHPY8GENTRIGGER_H
#define PHPYTHIA8_PHPY8GENTRIGGER_H

#include <iostream>
#include <string>
#include <vector>

namespace Pythia8
{
  class Pythia;
}

class PHPy8GenTrigger
{
 protected:
  //! constructor
  PHPy8GenTrigger(const std::string &name = "PHPy8GenTrigger");

 public:
  virtual ~PHPy8GenTrigger(){}

  virtual bool Apply(Pythia8::Pythia */*pythia*/)
  {
    std::cout << "PHPy8GenTrigger::Apply - in virtual function" << std::endl;
    return false;
  }

  virtual std::string GetName() { return m_Name; }

  std::vector<int> convertToInts(std::string s);
  int Verbosity() const { return m_Verbosity; }
  void Verbosity(int v) { m_Verbosity = v; }

 private:
  int m_Verbosity;
  std::string m_Name;
};

#endif
