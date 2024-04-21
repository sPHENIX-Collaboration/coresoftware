#ifndef TRIGGER_LL1OUTV1_H
#define TRIGGER_LL1OUTV1_H

#include "LL1Out.h"
#include "TriggerDefs.h"

#include <phool/PHObject.h>

#include <ostream>
#include <string>
#include <vector>

///
class LL1Outv1 : public LL1Out
{
 public:
  typedef std::map<unsigned int, std::vector<unsigned int>*> Map;
  typedef Map::const_iterator ConstIter;
  typedef Map::iterator Iter;
  typedef std::pair<Iter, Iter> Range;
  typedef std::pair<ConstIter, ConstIter> ConstRange;

  ///
  LL1Outv1();

  LL1Outv1(const std::string& triggertype, const std::string& ll1type);
  ///
  ~LL1Outv1() override;

  /// Clear Event from memory
  void Reset() override;

  void identify(std::ostream& os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  // Get Trigger Type
  std::string getLL1Type() const override { return _ll1_type; }

  // Set Trigger Type
  void setLL1Type(std::string& ll1type) override { _ll1_type = ll1type; }
  void setTriggerType(std::string& triggertype) override { _trigger_type = triggertype; }

  TriggerDefs::TriggerKey getTriggerKey() const override { return _trigger_key; }
  void setTriggerKey(TriggerDefs::TriggerKey key) override { _trigger_key = key; }

  std::vector<unsigned int>* GetTriggerBits() override { return _trigger_bits; }
  bool passesTrigger() override;

  void add_word(int key, std::vector<unsigned int>* trigger_words) override { _trigger_words[key] = trigger_words; }

  void set_event_number(unsigned int evt) override { _event_number = evt; }
  unsigned int get_event_number() override { return _event_number; }
  void set_clock_number(unsigned int clk) override { _clock_number = clk; }
  unsigned int get_clock_number() override { return _clock_number; }

  ConstRange getTriggerWords() const override;
  Range getTriggerWords() override;

 protected:
  int idx{0};
  unsigned int _event_number{0};
  unsigned int _clock_number{0};
  unsigned int _thresholds[10]{0};

  TriggerDefs::TriggerKey _trigger_key = TriggerDefs::TRIGGERKEYMAX;

  std::string _ll1_type{"NONE"};
  std::string _trigger_type{"NONE"};
  std::vector<unsigned int>* _trigger_bits{nullptr};
  Map _trigger_words;

 private:  // so the ClassDef does not show up with doc++
  ClassDefOverride(LL1Outv1, 1);
};

#endif
