#ifndef PHPYTHIA6_PHPY6GENTRIGGER_H
#define PHPYTHIA6_PHPY6GENTRIGGER_H

#include <iostream>
#include <string>
#include <vector>

namespace HepMC
{
  class GenEvent;
}

class PHPy6GenTrigger
{
 protected:
  //! constructor
  PHPy6GenTrigger(const std::string& name = "PHPy6GenTrigger");

 public:
  virtual ~PHPy6GenTrigger();

  virtual bool Apply(const HepMC::GenEvent* /*evt*/)
  {
    std::cout << "PHPy6GenTrigger::Apply - in virtual function" << std::endl;
    return false;
  }

  virtual std::string GetName()
  {
    return m_Name;
  }

  std::vector<int> convertToInts(std::string s);

  void Verbosity(int v) { m_Verbosity = v; }
  int Verbosity() const {return m_Verbosity;}

 private:
  int m_Verbosity = 0;
  std::string m_Name;
};

#endif
