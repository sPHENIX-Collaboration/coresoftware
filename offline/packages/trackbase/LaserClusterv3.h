/**
 * @file trackbase/LaserClusterv3.h
 * @author Ben Kimelman
 * @date July 2026
 * @brief Version 3 of CMFLashCluster
 */
#ifndef TRACKBASE_LASERCLUSTERV3_H
#define TRACKBASE_LASERCLUSTERV3_H

#include "LaserCluster.h"

#include <iostream>
#include <vector>

class PHObject;

/**
 * @brief Version 3 of LaserCluster
 *
 * Note - D. McGlinchey June 2018:
 *   CINT does not like "override", so ignore where CINT
 *   complains. Should be checked with ROOT 6 once
 *   migration occurs.
 */


class LaserClusterv3 : public LaserCluster
{
 public:
  //! ctor
  LaserClusterv3() = default;

  // PHObject virtual overloads
  void Reset() override {}
  int isValid() const override;
  PHObject* CloneMe() const override { return new LaserClusterv3(*this); }
 
  //! copy content from base class
  void CopyFrom( const LaserCluster& ) override;

  //! copy content from base class
  void CopyFrom( LaserCluster* source ) override
  { CopyFrom( *source ); }

  bool getFitMode() const override { return m_fitMode; }
  void setFitMode(bool fitMode) override { m_fitMode = fitMode; }

  unsigned int getLayerInt() const override { return m_posHardware[0]; }
  void setLayerInt(unsigned int layer) override { m_posHardware[0] = layer; }
  unsigned int getIPhiInt() const override { return m_posHardware[1]; }
  void setIPhiInt(unsigned int iphi) override { m_posHardware[1] = iphi; }
  unsigned int getITInt() const override { return m_posHardware[2]; }
  void setITInt(unsigned int it) override { m_posHardware[2] = it; }

  unsigned int getNhits() const override {return (unsigned int)m_hits.size();}

  //
  // cluster info
  //
  unsigned int getAdc() const override;

  void setNLayers(unsigned int nLayers) override { m_nLayers = nLayers; }
  unsigned int getNLayers() const override { return m_nLayers; }

  void setNIPhi(unsigned int nIPhi) override { m_nIPhi = nIPhi; }
  unsigned int getNIPhi() const override { return m_nIPhi; }

  void setNIT(unsigned int nIT) override { m_nIT = nIT; }
  unsigned int getNIT() const override { return m_nIT; }

  void setSDLayer(float SDLayer) override { m_SDLayer = SDLayer; }
  float getSDLayer() const override { return m_SDLayer; }

  void setSDIPhi(float SDIPhi) override { m_SDIPhi = SDIPhi; }
  float getSDIPhi() const override { return m_SDIPhi; }

  void setSDIT(float SDIT) override { m_SDIT = SDIT; }
  float getSDIT() const override { return m_SDIT; }

  void setSDWeightedLayer(float SDLayer) override { m_SDWeightedLayer = SDLayer; }
  float getSDWeightedLayer() const override { return m_SDWeightedLayer; }

  void setSDWeightedIPhi(float SDIPhi) override { m_SDWeightedIPhi = SDIPhi; }
  float getSDWeightedIPhi() const override { return m_SDWeightedIPhi; }

  void setSDWeightedIT(float SDIT) override { m_SDWeightedIT = SDIT; }
  float getSDWeightedIT() const override { return m_SDWeightedIT; }

  void addHit(TrkrDefs::hitsetkey hitsetkey, TrkrDefs::hitkey hitkey, uint16_t adc) override { m_hits.push_back(LaserClusterHitInfo(hitsetkey, hitkey, adc)); };
  LaserClusterHitInfo getHit(int hitIndex) const override { return m_hits[hitIndex]; };

  void identify(std::ostream& os = std::cout) const override;

 protected:

  std::vector<LaserClusterHitInfo> m_hits;
  
  unsigned int m_posHardware[3] = {std::numeric_limits<unsigned int>::max(), std::numeric_limits<unsigned int>::max(), std::numeric_limits<unsigned int>::max()};
  bool m_fitMode{false};

  /// number of TPC clusters used to create this central mebrane cluster
  unsigned int m_nhits = std::numeric_limits<unsigned int>::max();
  unsigned int m_nLayers = std::numeric_limits<unsigned int>::max();
  unsigned int m_nIPhi = std::numeric_limits<unsigned int>::max();
  unsigned int m_nIT = std::numeric_limits<unsigned int>::max();
  float m_SDLayer = std::numeric_limits<float>::quiet_NaN();
  float m_SDIPhi = std::numeric_limits<float>::quiet_NaN();
  float m_SDIT = std::numeric_limits<float>::quiet_NaN();
  float m_SDWeightedLayer = std::numeric_limits<float>::quiet_NaN();
  float m_SDWeightedIPhi = std::numeric_limits<float>::quiet_NaN();
  float m_SDWeightedIT = std::numeric_limits<float>::quiet_NaN();

  ClassDefOverride(LaserClusterv3, 1)
};

#endif //TRACKBASE_LASERCLUSTERV3_H
