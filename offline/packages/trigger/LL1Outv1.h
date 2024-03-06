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
  ~LL1Outv1() override;

  /// Clear Event from memory
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream 
  */
  void identify(std::ostream& os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  // Get Trigger Type
  std::string get_TriggerType() const {return _trigger_type;}

  // Set Trigger Type  
  void set_TriggerType(std::string &triggertype) { _trigger_type = triggertype;}

  void AddTriggerBits(unsigned int t_bits);
  unsigned int GetTriggerBits(unsigned int i_sample);

  void AddTriggerWords(vector<unsigned int> t_words);
  unsigned int GetTriggerWord(unsigned int i_word, unsigned int i_sample);

  void AddTriggerPrimitives(vector<unsigned int> t_prims);
  unsigned int GetTriggerPrimitive(unsigned int i_adc, unsigned int i_sample);

 protected:
  
  void Init() override;
  
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
