// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_RUNHEADERV1_H
#define FFAOBJECTS_RUNHEADERV1_H

#include "RunHeader.h"

#include <iostream>
#include <map>
#include <string>

class RunHeaderv1 : public RunHeader
{
 public:
  RunHeaderv1() = default;
  virtual ~RunHeaderv1() = default;

  void Reset() override;
  void identify(std::ostream &os = std::cout) const override;
  int isValid() const override;

  int get_RunNumber() const override { return RunNumber; }
  void set_RunNumber(const int run) override
  {
    RunNumber = run;
    return;
  }

  void set_floatval(const std::string &name, const float fval) override;
  float get_floatval(const std::string &name) const override;

  void set_intval(const std::string &name, const int ival) override;
  int get_intval(const std::string &name) const override;

 private:
  int RunNumber = 0;
  std::map<std::string, int> m_IntRunProperties;
  std::map<std::string, float> m_FloatRunProperties;

// rootcling and clang complain about inconsistent overrides in the ClassDef
// this can be supressed with ignoring -Winconsistent-missing-override
// this pragma is not known to gcc, so we need an #ifdef __clang__ here
#pragma GCC diagnostic push
#if defined(__clang__)
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#endif
  ClassDef(RunHeaderv1, 1)
#pragma GCC diagnostic pop
};

#endif /* __RUNHEADERV1_H */
