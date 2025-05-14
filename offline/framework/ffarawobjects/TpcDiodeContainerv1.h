#ifndef FUN4ALLRAW_TPCDIODECONTAINERV1_H
#define FUN4ALLRAW_TPCDIODECONTAINERV1_H

#include "TpcDiodeContainer.h"

class TpcDiode;
class TClonesArray;

class TpcDiodeContainerv1 : public TpcDiodeContainer
{
 public:
  TpcDiodeContainerv1();
  ~TpcDiodeContainerv1() override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  TpcDiode *AddDiode() override;
  TpcDiode *AddDiode(TpcDiode *tpcdiode) override;
  unsigned int get_ndiodes() override;
  TpcDiode *get_diode(unsigned int index) override;
  unsigned int get_Laser() override;
  std::vector<TpcDiode *> get_PO1() override;
  std::vector<TpcDiode *> get_PO2() override;
  std::vector<TpcDiode *> get_EGG() override;

 private:
  TClonesArray *TpcDiodesTCArray = nullptr;

  ClassDefOverride(TpcDiodeContainerv1, 1)
};

#endif
