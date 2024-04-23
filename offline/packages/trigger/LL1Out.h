#ifndef TRIGGER_LL1OUT_H__
#define TRIGGER_LL1OUT_H__

#include "TriggerDefs.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <ostream>
///
class LL1Out : public PHObject
{
 public:
  typedef std::map<unsigned int, std::vector<unsigned int>*> Map;
  typedef Map::const_iterator ConstIter;
  typedef Map::iterator Iter;
  typedef std::pair<Iter, Iter> Range;
  typedef std::pair<ConstIter, ConstIter> ConstRange;

  ///
  LL1Out() = default;
  ///
  virtual ~LL1Out() override = default;
  /// Clear Event from memory
  virtual void Reset() override { return; }

  void identify(std::ostream& os = std::cout) const override;
  int isValid() const override { return 1; }

  virtual std::string getLL1Type() const { return "none"; }

  virtual void setLL1Type(std::string& /*ll1type*/) {}
  virtual void setTriggerType(std::string& /*triggertype*/) {}

  virtual TriggerDefs::TriggerKey getTriggerKey() const { return 0; }
  virtual void setTriggerKey(TriggerDefs::TriggerKey /*key*/) {}

  virtual std::vector<unsigned int>* GetTriggerBits() { return nullptr; }

  virtual bool passesTrigger() { return 0; }

  virtual void add_word(int /*key*/, std::vector<unsigned int>* /*trigger_words*/) {}

  virtual void set_event_number(unsigned int /*evt*/) {}
  virtual unsigned int get_event_number() { return 0; }

  virtual void set_clock_number(unsigned int /*clk*/) {}
  virtual unsigned int get_clock_number() { return 0; }

  virtual ConstRange getTriggerWords() const;
  virtual Range getTriggerWords();

 private:  // so the ClassDef does not show up with doc++
  ClassDefOverride(LL1Out, 1);
};

#endif
