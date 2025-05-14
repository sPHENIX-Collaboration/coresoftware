#ifndef FUN4ALLRAW_TPCDIODECONTAINER_H
#define FUN4ALLRAW_TPCDIODECONTAINER_H

#include <phool/PHObject.h>

class TpcDiode;

class TpcDiodeContainer : public PHObject
{
 public:
  TpcDiodeContainer() = default;
  virtual ~TpcDiodeContainer() = default;

  virtual TpcDiode *AddDiode() { return nullptr; }
  virtual TpcDiode *AddDiode(TpcDiode *) { return nullptr; }
  virtual unsigned int get_ndiodes() { return 0; }
  virtual TpcDiode *get_diode(unsigned int) { return nullptr; }
  virtual unsigned int get_Laser() { return 0; }
  virtual std::vector<TpcDiode *> get_PO1() { return std::vector<TpcDiode *>{nullptr}; }
  virtual std::vector<TpcDiode *> get_PO2() { return std::vector<TpcDiode *>{nullptr}; }
  virtual std::vector<TpcDiode *> get_EGG() { return std::vector<TpcDiode *>{nullptr}; }
  virtual void setStatus(const unsigned int) { return; }
  virtual unsigned int getStatus() const { return 0; }
  // virtual void setBco(const uint64_t) { return; }
  // virtual uint64_t getBco() const { return 0; }

 private:
  ClassDefOverride(TpcDiodeContainer, 0)
};

#endif
