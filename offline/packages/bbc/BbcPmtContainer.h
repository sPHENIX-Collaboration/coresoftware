//  virtual Bbc PMT Container class

#ifndef __BBCPMTCONTAINER_H__
#define __BBCPMTCONTAINER_H__

#include "phool/PHObject.h"
#include <iostream>

///
class BbcPmtContainer : public PHObject
{
 public:
  /// dtor
  virtual ~BbcPmtContainer() {}

  /** identify Function from PHObject
      @param os Output Stream 
   */
  virtual void identify(std::ostream& os = std::cout) const override;

  /// Clear Event
  virtual void Reset() override;

  /// isValid returns non zero if object contains valid data
  virtual int isValid() const override;

  /** set number of PMTs for Bbc
      @param ival Number of Bbc Pmt's
   */
  virtual void set_npmt(const Short_t ival);

  /// get Number of Bbc Pmt's
  virtual Short_t get_npmt() const;

  /** get Adc of Pmt iPmt in TClonesArray
      @param iPmt no of Pmt in TClonesArray
   */
  virtual Float_t get_adc(const int iPmt) const;

  /** get Tdc0 of Pmt iPmt in TClonesArray
      @param iPmt no of Pmt in TClonesArray
   */
  virtual Float_t get_tdc0(const int iPmt) const;

  /** get Tdc1 of Pmt iPmt in TClonesArray
      @param iPmt no of Pmt in TClonesArray
   */
  virtual Float_t get_tdc1(const int iPmt) const;

  /** Add Bbc Raw hit object to TCLonesArray
      @param ipmt Pmt id
      @param adc  Adc value
      @param tdc0 Tdc0 value
      @param tdc1 Tdc1 value
  */
  virtual void AddBbcPmt(const Short_t ipmt, const Float_t adc, const Float_t tdc0, const Float_t tdc1);

 private:
  void virtual_warning(const char *funcname) const;

  ClassDefOverride(BbcPmtContainer,1)
};

#endif
