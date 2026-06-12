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
  void setPosition(int coor, float xi) override
  {
    if (coor >= 0 && coor < 2)
    {
      m_local[coor] = xi;
    }
  }
  float getLocalX() const override { return m_local[0]; }
  void setLocalX(float loc0) override { m_local[0] = loc0; }
  float getLocalY() const override { return m_local[1]; }
  void setLocalY(float loc1) override { m_local[1] = loc1; }

  TrkrDefs::subsurfkey getSubSurfKey() const override { return m_subsurfkey; }
  void setSubSurfKey(TrkrDefs::subsurfkey id) override { m_subsurfkey = id; }

  //
  // cluster info
  //
  unsigned int getAdc() const override { return m_adc; }
  void setAdc(unsigned int adc) override { m_adc = adc; }

  unsigned int getMaxAdc() const override { return m_maxadc; }
  void setMaxAdc(uint16_t maxadc) override { m_maxadc = maxadc; }

  unsigned int getCenAdc() const override { return m_cenadc; }
  void setCenAdc(uint16_t cenadc) { m_cenadc = cenadc; }

  float getPadCen() const override { return m_padcen; }
  void setPadCen(float padcen) { m_padcen = padcen; }

  float getTBinCen() const override { return m_tbincen; }
  void setTBinCen(float tbincen) { m_tbincen = tbincen; }

  float getPadMax() const override { return m_padmax; }
  void setPadMax(float padmax) { m_padmax = padmax; }

  float getTBinMax() const override { return m_tbinmax; }
  void setTBinMax(float tbinmax) { m_tbinmax = tbinmax; }

  //
  // convenience interface
  //
  float getRPhiError() const override { return m_phierr; }
  float getZError() const override { return m_zerr; }

  void setPhiError(float phierror) { m_phierr = phierror; }
  void setZError(float zerror) { m_zerr = zerror; }

  char getSize() const override { return m_phisize * m_zsize; }
  // void setSize(char size) { m_size = size; }

  float getRSize() const override { return (float) m_rsize; }
  void setRSize(unsigned char rsize) { m_rsize = rsize; }

  float getPhiSize() const override { return (float) m_phisize; }
  void setPhiSize(char phisize) { m_phisize = phisize; }

  float getZSize() const override { return (float) m_zsize; }
  void setZSize(char zsize) { m_zsize = zsize; }

  char getOverlap() const override { return m_overlap; }
  void setOverlap(char overlap) override { m_overlap = overlap; }

  char getEdge() const override { return m_edge; }
  void setEdge(char edge) override { m_edge = edge; }

  char getSLEdge() const override { return m_sledge; }
  void setSLEdge(char sledge) { m_sledge = sledge; }

  char getSREdge() const override { return m_sredge; }
  void setSREdge(char sredge) { m_sredge = sredge; }

  char getTLEdge() const override { return m_tledge; }
  void setTLEdge(char tledge) { m_tledge = tledge; }

  char getTREdge() const override { return m_tredge; }
  void setTREdge(char tredge) { m_tredge = tredge; }

  char getDLEdge() const override { return m_dledge; }
  void setDLEdge(char dledge) { m_dledge = dledge; }

  char getDREdge() const override { return m_dredge; }
  void setDREdge(char dredge) { m_dredge = dredge; }

  char getHLEdge() const override { return m_hledge; }
  void setHLEdge(char hledge) { m_hledge = hledge; }

  char getHREdge() const override { return m_hredge; }
  void setHREdge(char hredge) { m_hredge = hredge; }

  int getSLMix() const override { return m_slmix; }
  void setSLMix(char slmix) { m_slmix = slmix; }

  int getSRMix() const override { return m_srmix; }
  void setSRMix(char srmix) { m_srmix = srmix; }

  int getTLMix() const override { return m_tlmix; }
  void setTLMix(char tlmix) { m_tlmix = tlmix; }

  int getTRMix() const override { return m_trmix; }
  void setTRMix(char trmix) { m_trmix = trmix; }

  float getPhiBinLo() const override { return m_phibinlo; }
  void setPhiBinLo(float phibinlo) { m_phibinlo = phibinlo; }

  float getPhiBinHi() const override { return m_phibinhi; }
  void setPhiBinHi(float phibinhi) { m_phibinhi = phibinhi; }

  float getTBinLo() const override { return m_tbinlo; }
  void setTBinLo(float tbinlo) { m_tbinlo = tbinlo; }

  float getTBinHi() const override { return m_tbinhi; }
  void setTBinHi(float tbinhi) { m_tbinhi = tbinhi; }

  float getPadPhase() const override { return m_padphase; }
  void setPadPhase(float padphase) { m_padphase = padphase; }

  float getTBinPhase() const override { return m_tbinphase; }
  void setTBinPhase(float tbinphase){ m_tbinphase = tbinphase; }

 private:
  float m_local[2]{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
  //< 2D local position [cm] 2 * 32 64bit  - cumul 1*64
  TrkrDefs::subsurfkey m_subsurfkey {TrkrDefs::SUBSURFKEYMAX};  //< unique identifier for hitsetkey-surface maps 16 bit
  float m_phierr{0};
  float m_zerr{0};
  unsigned short int m_adc{0};     //< cluster sum adc 16
  unsigned short int m_maxadc{0};  //< cluster max adc 16
  unsigned short int m_cenadc{0};  //< cluster centroid adc 16
  float m_padcen{0};
  float m_tbincen{0};
  float m_padmax{0};
  float m_tbinmax{0};
  unsigned char m_rsize{0};        // 8bit
  char m_phisize{0};               // 8bit
  char m_zsize{0};                 // 8bit
  char m_overlap{0};               // 8bit
  char m_edge{0};                  // 8bit - cumul 2*64
  char m_sledge{0};                // 8bit
  char m_sredge{0};                // 8bit
  char m_tledge{0};                // 8bit
  char m_tredge{0};                // 8bit
  char m_dledge{0};                // 8bit
  char m_dredge{0};                // 8bit
  char m_hledge{0};                // 8bit
  char m_hredge{0};                // 8bit
  char m_slmix{0};                 // 8bit
  char m_srmix{0};                 // 8bit
  char m_tlmix{0};                 // 8bit
  char m_trmix{0};                 // 8bit
  float m_phibinlo{0};
  float m_phibinhi{0};
  float m_tbinlo{0};
  float m_tbinhi{0};
  float m_padphase{0};
  float m_tbinphase{0};

  ClassDefOverride(TrkrClusterv6, 1)
};

#endif  // TRACKBASE_TRKRCLUSTERV6_H
