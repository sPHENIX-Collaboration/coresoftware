#ifndef __LL1OUTV1_H
#define __LL1OUTV1_H

#include <string>
#include <ostream>
#include <phool/PHObject.h>
#include "LL1Out.h"
#include <vector>

using namespace std;
///
class LL1Outv1 : public LL1Out
{
 public:

  ///
  LL1Outv1();
  ///
  virtual ~LL1Outv1() override;

  /// Clear Event from memory
  virtual void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream 
  */
  virtual void identify(std::ostream& os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  virtual int isValid() const override;

  // Get Trigger Type
  virtual std::string get_TriggerType() const {return _trigger_type;}

  // Set Trigger Type  
  virtual void set_TriggerType(std::string &triggertype) { _trigger_type = triggertype;}

  virtual void AddTriggerBits(unsigned int t_bits);
  virtual unsigned int GetTriggerBits(unsigned int i_sample);

  virtual void AddTriggerWords(vector<unsigned int> t_words);
  virtual unsigned int GetTriggerWord(unsigned int i_word, unsigned int i_sample);

  virtual void AddTriggerPrimitives(vector<unsigned int> t_prims);
  virtual unsigned int GetTriggerPrimitive(unsigned int i_adc, unsigned int i_sample);

 protected:
  
  virtual void Init() override;
  
 private:
  
  std::string _trigger_type;

  int idx;

  vector<unsigned int> _trigger_bits;
  vector<vector<unsigned int>> _trigger_prims;
  vector<vector<unsigned int>> _trigger_words;
  unsigned int _thresholds[10];

 private: // so the ClassDef does not show up with doc++
  ClassDefOverride(LL1Outv1,1);
};

#endif
