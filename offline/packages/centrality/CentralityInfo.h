#ifndef CENTRALITY_IO_CENTRALITYINFO_H
#define CENTRALITY_IO_CENTRALITYINFO_H

//===========================================================
/// \file CentralityInfo.h
/// \brief Lightweight centrality information storage node
/// \author Dennis V. Perepelitsa
//===========================================================

#include <phool/PHObject.h>

class CentralityInfo : public PHObject
{
 public:
  ~CentralityInfo() override{};

  void identify(std::ostream &os = std::cout) const override { os << "CentralityInfo base class" << std::endl; };
  int isValid() const override { return 0; }

  enum PROP
  {

    //! Minimum Bias Detector (MBD) North-side charge sum
    mbd_N = 0,
    //! MBD South-side charge sum
    mbd_S = 1,
    //! MBD North+South charge sum
    mbd_NS = 2,

    //! sPHENIX Event Plane Detector (sEPD) North-side energy sum
    epd_N = 3,
    //! sEPD South-side energy sum
    epd_S = 4,
    //! sEPD North+South energy sum
    epd_NS = 5,

    //! Impact parameter (b) in HIJING event
    bimp = 6

  };

  virtual bool has_quantity(const PROP /*prop_id*/) const { return false; }
  virtual float get_quantity(const PROP /*prop_id*/) const { return -99; }
  virtual void set_quantity(const PROP /*prop_id*/, const float /*value*/) { return; }

  virtual bool has_centile(const PROP /*prop_id*/) const { return false; }
  virtual float get_centile(const PROP /*prop_id*/) const { return -99; }
  virtual void set_centile(const PROP /*prop_id*/, const float /*value*/) { return; }

 protected:
  CentralityInfo() {}

 private:
  ClassDefOverride(CentralityInfo, 1);
};

#endif
