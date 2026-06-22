/**
 * @file trackbase/TrkrClusterv6.h
 * @author Ishan Goel
 * @date May 2026
 * @brief Version 6 of TrkrCluster
 */
#ifndef TRACKBASE_TRKRCLUSTERV6_H
#define TRACKBASE_TRKRCLUSTERV6_H

#include "TrkrCluster.h"
#include "TrkrDefs.h"

#include <limits>
#include <iostream>

class PHObject;

/**
 * @brief Version 6 of TrkrCluster
 *
 * This version of TrkrCluster is blown up to contain a maximum of information
 */

class TrkrClusterv6 : public TrkrCluster
{
 public:
  //! ctor
  TrkrClusterv6() = default;

  //! dtor
  ~TrkrClusterv6() override = default;

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = TrkrClusterv6(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new TrkrClusterv6(*this); }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;

  //! copy content from base class
  void CopyFrom(const TrkrCluster&) override;

  //! copy content from base class
  void CopyFrom(TrkrCluster* source) override
  {
    if (!source)
    {
      return;
    }
    CopyFrom(*source);
  }

  //
  // cluster position
  //
  float getPosition(int coor) const override
  {
    return (coor >= 0 && coor < 2) ? m_local[coor] : std::numeric_limits<float>::quiet_NaN();
  }
  void setPosition(const int coor, const float xi) override
  {
    if (coor >= 0 && coor < 2)
    {
      m_local[coor] = xi;
    }
  }
  float getLocalX() const override { return m_local[0]; }
  void setLocalX(const float loc0) override { m_local[0] = loc0; }
  float getLocalY() const override { return m_local[1]; }
  void setLocalY(const float loc1) override { m_local[1] = loc1; }

  TrkrDefs::subsurfkey getSubSurfKey() const override { return m_subsurfkey; }
  void setSubSurfKey(const TrkrDefs::subsurfkey id) override { m_subsurfkey = id; }

  //
  // cluster info
  //
  unsigned int getAdc() const override { return m_adc; }
  void setAdc(const unsigned int adc) override { m_adc = adc; }

  unsigned int getMaxAdc() const override { return m_maxadc; }
  void setMaxAdc(const uint16_t maxadc) override { m_maxadc = maxadc; }

  unsigned int getCenAdc() const override { return m_cenadc; }
  void setCenAdc(const uint16_t cenadc) override { m_cenadc = cenadc; }

  float getPadCen() const override { return m_padcen; }
  void setPadCen(const float padcen) override { m_padcen = padcen; }

  float getTBinCen() const override { return m_tbincen; }
  void setTBinCen(const float tbincen) override { m_tbincen = tbincen; }

  int getPadMax() const override { return m_padmax; }
  void setPadMax(const int padmax) override { m_padmax = padmax; }

  int getTBinMax() const override { return m_tbinmax; }
  void setTBinMax(const int tbinmax) override { m_tbinmax = tbinmax; }

  //
  // convenience interface
  //
  float getRPhiError() const override { return m_phierr; }
  float getZError() const override { return m_zerr; }

  void setPhiError(const float phierror) override { m_phierr = phierror; }
  void setZError(const float zerror) override { m_zerr = zerror; }

  char getSize() const override { return m_phisize * m_zsize; }
  // void setSize(const char size) { m_size = size; }

  float getRSize() const override { return (float) m_rsize; }
  void setRSize(const char rsize) override { m_rsize = rsize; }

  float getPhiSize() const override { return (float) m_phisize; }
  void setPhiSize(const char phisize) override { m_phisize = phisize; }

  float getZSize() const override { return (float) m_zsize; }
  void setZSize(const char zsize) override { m_zsize = zsize; }

  char getOverlap() const override { return m_overlap; }
  void setOverlap(const char overlap) override { m_overlap = overlap; }

  char getEdge() const override { return m_edge; }
  void setEdge(const char edge) override { m_edge = edge; }

  char getSLEdge() const override { return m_sledge; }
  void setSLEdge(const char sledge) override { m_sledge = sledge; }

  char getSREdge() const override { return m_sredge; }
  void setSREdge(const char sredge) override { m_sredge = sredge; }

  char getTLEdge() const override { return m_tledge; }
  void setTLEdge(const char tledge) override { m_tledge = tledge; }

  char getTREdge() const override { return m_tredge; }
  void setTREdge(const char tredge) override { m_tredge = tredge; }

  char getDLEdge() const override { return m_dledge; }
  void setDLEdge(const char dledge) override { m_dledge = dledge; }

  char getDREdge() const override { return m_dredge; }
  void setDREdge(const char dredge) override { m_dredge = dredge; }

  char getHLEdge() const override { return m_hledge; }
  void setHLEdge(const char hledge) override { m_hledge = hledge; }

  char getHREdge() const override { return m_hredge; }
  void setHREdge(const char hredge) override { m_hredge = hredge; }

  int getSLMix() const override { return m_slmix; }
  void setSLMix(const char slmix) override { m_slmix = slmix; }

  int getSRMix() const override { return m_srmix; }
  void setSRMix(const char srmix) override { m_srmix = srmix; }

  int getTLMix() const override { return m_tlmix; }
  void setTLMix(const char tlmix) override { m_tlmix = tlmix; }

  int getTRMix() const override { return m_trmix; }
  void setTRMix(const char trmix) override { m_trmix = trmix; }

  char getPhiBinLo() const override { return m_phibinlo; }
  void setPhiBinLo(const char phibinlo) override { m_phibinlo = phibinlo; }

  char getPhiBinHi() const override { return m_phibinhi; }
  void setPhiBinHi(const char phibinhi) override { m_phibinhi = phibinhi; }

  char getTBinLo() const override { return m_tbinlo; }
  void setTBinLo(const char tbinlo) override { m_tbinlo = tbinlo; }

  char getTBinHi() const override { return m_tbinhi; }
  void setTBinHi(const char tbinhi) override { m_tbinhi = tbinhi; }

  float getPadPhase() const override { return m_padphase; }
  void setPadPhase(const float padphase) override { m_padphase = padphase; }

  float getTBinPhase() const override { return m_tbinphase; }
  void setTBinPhase(const float tbinphase) override { m_tbinphase = tbinphase; }

 private:
  float m_local[2]{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
  //< 2D local position [cm] 2 * 32 64bit  - cumul 1*64
  TrkrDefs::subsurfkey m_subsurfkey {TrkrDefs::SUBSURFKEYMAX};  //< unique identifier for hitsetkey-surface maps 16 bit
  float m_phierr{0};
  float m_zerr{0};
  unsigned short m_adc{0};     //< cluster sum adc 16
  unsigned short m_maxadc{0};  //< cluster max adc 16
  unsigned short m_cenadc{0};  //< cluster centroid adc 16
  float m_padcen{0};
  float m_tbincen{0};
  int m_padmax{0};
  int m_tbinmax{0};
  char m_rsize{0};        
  char m_phisize{0};               
  char m_zsize{0};                 
  char m_overlap{0};               
  char m_edge{0};                   
  char m_sledge{0};                
  char m_sredge{0};                
  char m_tledge{0};                
  char m_tredge{0};                
  char m_dledge{0};                
  char m_dredge{0};                
  char m_hledge{0};                
  char m_hredge{0};                
  char m_slmix{0};                 
  char m_srmix{0};                 
  char m_tlmix{0};                 
  char m_trmix{0};                 
  char m_phibinlo{0};
  char m_phibinhi{0};
  char m_tbinlo{0};
  char m_tbinhi{0};
  float m_padphase{0};
  float m_tbinphase{0};

  ClassDefOverride(TrkrClusterv6, 1)
};

#endif  // TRACKBASE_TRKRCLUSTERV6_H
