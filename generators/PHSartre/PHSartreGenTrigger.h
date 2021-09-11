#ifndef PHSARTRE_PHSARTREGENTRIGGER_H
#define PHSARTRE_PHSARTREGENTRIGGER_H

#include <iostream>
#include <string>
#include <vector>

class Event;

class PHSartreGenTrigger
{
 protected:
  //! constructor
  PHSartreGenTrigger(const std::string &name = "PHSartreGenTrigger");

 public:
  virtual ~PHSartreGenTrigger();

  virtual bool Apply(Event */*event*/)
  {
    std::cout << "PHSartreGenTrigger::Apply - in virtual function" << std::endl;
    return false;
  }

  virtual std::string GetName() { return m_Name; }

  std::vector<int> convertToInts(std::string s);
  int Verbosity() const { return m_Verbosity; }
  void Verbosity(int v) { m_Verbosity = v; }

 protected:
 private:
  int m_Verbosity;
  std::string m_Name;
};

#endif
