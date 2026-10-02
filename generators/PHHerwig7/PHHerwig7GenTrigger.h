#ifndef PHHERWIG7_PHHERWIG7GENTRIGGER_H
#define PHHERWIG7_PHHERWIG7GENTRIGGER_H

#include <iostream>
#include <string>
#include <vector>

namespace HepMC
{
  class GenEvent;
}

class PHHerwig7GenTrigger
{
 protected:
  //! constructor
  PHHerwig7GenTrigger(const std::string &name = "PHHerwig7GenTrigger");

 public:
  virtual ~PHHerwig7GenTrigger() {}

  virtual bool Apply(HepMC::GenEvent * /*herwig*/)
  {
    std::cout << "PHHerwig7GenTrigger::Apply - in virtual function" << std::endl;
    return false;
  }

  virtual std::string GetName() { return m_Name; }

  std::vector<int> convertToInts(const std::string &s);
  int Verbosity() const { return m_Verbosity; }
  void Verbosity(int v) { m_Verbosity = v; }

 private:
  int m_Verbosity;
  std::string m_Name;
};

#endif
