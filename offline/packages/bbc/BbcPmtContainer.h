// Tell emacs that this is a C++ source
//  -*- C++ -*-.
//  virtual Bbc PMT Container class

#ifndef BBC_BBCPMTCONTAINER_H
#define BBC_BBCPMTCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>
#include <string>

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
  virtual void set_npmt(const short ival);

  /// get Number of Bbc Pmt's
  virtual short get_npmt() const;

  /** get id of Pmt iPmt in TClonesArray
      @param iPmt no of Pmt in TClonesArray
   */
  virtual short get_pmt(const int iPmt) const;

  /** get Adc of Pmt iPmt in TClonesArray
      @param iPmt no of Pmt in TClonesArray
   */
  virtual float get_adc(const int iPmt) const;

  /** get Tdc0 of Pmt iPmt in TClonesArray
      @param iPmt no of Pmt in TClonesArray
   */
  virtual float get_tdc0(const int iPmt) const;

  /** get Tdc1 of Pmt iPmt in TClonesArray
      @param iPmt no of Pmt in TClonesArray
   */
  virtual float get_tdc1(const int iPmt) const;

  /** Add Bbc Raw hit object to TCLonesArray
      @param ipmt Pmt id
      @param adc  Adc value
      @param tdc0 Tdc0 value
      @param tdc1 Tdc1 value
  */
  virtual void AddBbcPmt(const short ipmt, const float adc, const float tdc0, const float tdc1);

 private:
  void virtual_warning(const std::string& funcname) const;

  ClassDefOverride(BbcPmtContainer, 1)
};

#endif
