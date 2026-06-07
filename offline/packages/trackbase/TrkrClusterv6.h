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

#include <climits>
#include <cstdint>
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
  uint16_t getAdc() const override { return m_adc; }
  void setAdc(uint16_t adc) override { m_adc = adc; }

  uint16_t getMaxAdc() const override { return m_maxadc; }
  void setMaxAdc(uint16_t maxadc) override { m_maxadc = maxadc; }

  uint16_t getCenAdc() const override { return m_cenadc; }
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

  /// deprecated global funtions with a warning
  float getX() const override
  {
    std::cout << "Deprecated getx trkrcluster function!" << std::endl;
    return std::numeric_limits<float>::quiet_NaN();
  }
  float getY() const override
  {
    std::cout << "Deprecated gety trkrcluster function!" << std::endl;
    return std::numeric_limits<float>::quiet_NaN();
  }
  float getZ() const override
  {
    std::cout << "Deprecated getz trkrcluster function!" << std::endl;
    return std::numeric_limits<float>::quiet_NaN();
  }
  void setX(float) override
  {
    std::cout << "Deprecated setx trkrcluster function!" << std::endl;
  }
  void setY(float) override
  {
    std::cout << "Deprecated sety trkrcluster function!" << std::endl;
  }
  void setZ(float) override
  {
    std::cout << "Deprecated setz trkrcluster function!" << std::endl;
  }
  float getSize(unsigned int, unsigned int) const override
  {
    std::cout << "Deprecated getsize trkrcluster function!" << std::endl;
    return std::numeric_limits<float>::quiet_NaN();
  }
  void setSize(unsigned int, unsigned int, float) override
  {
    std::cout << "Deprecated setsize trkrcluster function!" << std::endl;
  }
  float getError(unsigned int, unsigned int) const override
  {
    std::cout << "Deprecated geterr trkrcluster function!" << std::endl;
    return std::numeric_limits<float>::quiet_NaN();
  }
  void setError(unsigned int, unsigned int, float) override
  {
    std::cout << "Deprecated seterr trkrcluster function!" << std::endl;
  }

  uint8_t getSize() const override { return m_phisize * m_zsize; }
  // void setSize(uint8_t size) { m_size = size; }

  float getRSize() const override { return (float) m_rsize; }
  void setRSize(uint8_t rsize) { m_rsize = rsize; }

  float getPhiSize() const override { return (float) m_phisize; }
  void setPhiSize(uint8_t phisize) { m_phisize = phisize; }

  float getZSize() const override { return (float) m_zsize; }
  void setZSize(uint8_t zsize) { m_zsize = zsize; }

  uint8_t getOverlap() const override { return m_overlap; }
  void setOverlap(uint8_t overlap) override { m_overlap = overlap; }

  uint8_t getEdge() const override { return m_edge; }
  void setEdge(uint8_t edge) override { m_edge = edge; }

  uint8_t getSLEdge() const override { return m_sledge; }
  void setSLEdge(uint8_t sledge) { m_sledge = sledge; }

  uint8_t getSREdge() const override { return m_sredge; }
  void setSREdge(uint8_t sredge) { m_sredge = sredge; }

  uint8_t getTLEdge() const override { return m_tledge; }
  void setTLEdge(uint8_t tledge) { m_tledge = tledge; }

  uint8_t getTREdge() const override { return m_tredge; }
  void setTREdge(uint8_t tredge) { m_tredge = tredge; }

  uint8_t getDLEdge() const override { return m_dledge; }
  void setDLEdge(uint8_t dledge) { m_dledge = dledge; }

  uint8_t getDREdge() const override { return m_dredge; }
  void setDREdge(uint8_t dredge) { m_dredge = dredge; }

  uint8_t getHLEdge() const override { return m_hledge; }
  void setHLEdge(uint8_t hledge) { m_hledge = hledge; }

  uint8_t getHREdge() const override { return m_hredge; }
  void setHREdge(uint8_t hredge) { m_hredge = hredge; }

  int getSLMix() const override { return m_slmix; }
  void setSLMix(uint8_t slmix) { m_slmix = slmix; }

  int getSRMix() const override { return m_srmix; }
  void setSRMix(uint8_t srmix) { m_srmix = srmix; }

  int getTLMix() const override { return m_tlmix; }
  void setTLMix(uint8_t tlmix) { m_tlmix = tlmix; }

  int getTRMix() const override { return m_trmix; }
  void setTRMix(uint8_t trmix) { m_trmix = trmix; }

  int getPhiBinLo() const override { return m_phibinlo; }
  void setPhiBinLo(int phibinlo) { m_phibinlo = phibinlo; }

  int getPhiBinHi() const override { return m_phibinhi; }
  void setPhiBinHi(int phibinhi) { m_phibinhi = phibinhi; }

  int getTBinLo() const override { return m_tbinlo; }
  void setTBinLo(int tbinlo) { m_tbinlo = tbinlo; }

  int getTBinHi() const override { return m_tbinhi; }
  void setTBinHi(int tbinhi) { m_tbinhi = tbinhi; }

  float getPadPhase() const override { return m_padphase; }
  void setPadPhase(float padphase) { m_padphase = padphase; }

  float getTBinPhase() const override { return m_tbinphase; }
  void setTBinPhase(float tbinphase){ m_tbinphase = tbinphase; }

  // float getPhiSize() const override
  //{ std::cout << "Deprecated size function"<< std::endl; return std::numeric_limits<float>::quiet_NaN();}
  // float getZSize() const override
  //{std::cout << "Deprecated size function" << std::endl; return std::numeric_limits<float>::quiet_NaN();}
  // float getPhiError() const override
  //{ std::cout << "Deprecated getPhiError function"<< std::endl; return std::numeric_limits<float>::quiet_NaN();}

 private:
  float m_local[2]{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
  //< 2D local position [cm] 2 * 32 64bit  - cumul 1*64
  TrkrDefs::subsurfkey m_subsurfkey{TrkrDefs::SUBSURFKEYMAX};  //< unique identifier for hitsetkey-surface maps 16 bit
  float m_phierr{0};
  float m_zerr{0};
  uint16_t m_adc{0};               //< cluster sum adc 16
  uint16_t m_maxadc{0};            //< cluster max adc 16
  uint16_t m_cenadc{0};            //< cluster centroid adc 16
  float m_padcen{0};
  float m_tbincen{0};
  float m_padmax{0};
  float m_tbinmax{0};
  uint8_t m_rsize{0};              // 8bit
  uint8_t m_phisize{0};            // 8bit
  uint8_t m_zsize{0};              // 8bit
  uint8_t m_overlap{0};            // 8bit
  uint8_t m_edge{0};               // 8bit - cumul 2*64
  uint8_t m_sledge{0};             // 8bit
  uint8_t m_sredge{0};             // 8bit
  uint8_t m_tledge{0};             // 8bit
  uint8_t m_tredge{0};             // 8bit
  uint8_t m_dledge{0};             // 8bit
  uint8_t m_dredge{0};             // 8bit
  uint8_t m_hledge{0};             // 8bit
  uint8_t m_hredge{0};             // 8bit
  uint8_t m_slmix{0};              // 8bit
  uint8_t m_srmix{0};              // 8bit
  uint8_t m_tlmix{0};              // 8bit
  uint8_t m_trmix{0};              // 8bit
  int m_phibinlo{0};
  int m_phibinhi{0};
  int m_tbinlo{0};
  int m_tbinhi{0};
  float m_padphase{0};
  float m_tbinphase{0};

  ClassDefOverride(TrkrClusterv6, 1)
};

#endif  // TRACKBASE_TRKRCLUSTERV6_H
