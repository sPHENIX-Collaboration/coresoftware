// PHJewelGenTrigger -- base class of the truth triggers PHJewel applies inside
// its generation loop (PHHerwig7GenTrigger / PHPy8GenTrigger pattern).
// Register with PHJewel::register_trigger(); PHJewel owns nothing and
// combines the results with set_trigger_OR() (default) or set_trigger_AND().
#ifndef PHJEWEL_PHJEWELGENTRIGGER_H
#define PHJEWEL_PHJEWELGENTRIGGER_H

#include <iostream>
#include <string>

namespace HepMC
{
  class GenEvent;
}

class PHJewelGenTrigger
{
 protected:
  explicit PHJewelGenTrigger(const std::string &name = "PHJewelGenTrigger")
    : m_Name(name)
  {
  }

 public:
  virtual ~PHJewelGenTrigger() = default;

  virtual bool Apply(const HepMC::GenEvent * /*evt*/)
  {
    std::cout << "PHJewelGenTrigger::Apply - in virtual function" << std::endl;
    return false;
  }

  virtual void PrintConfig() const {}

  const std::string &GetName() const { return m_Name; }
  int Verbosity() const { return m_Verbosity; }
  void Verbosity(int v) { m_Verbosity = v; }

 private:
  int m_Verbosity{0};
  std::string m_Name;
};

#endif
