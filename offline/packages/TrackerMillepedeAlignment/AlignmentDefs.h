#ifndef ALIGNMENTDEFS_H
#define ALIGNMENTDEFS_H

#include <trackbase/ActsGeometry.h>

namespace AlignmentDefs
{
  enum siliconGrp {snsr, stv, brrl};
  enum tpcGrp {htst, sctr, tp};
  enum mmsGrp {tl, mm};
  
  //! Map relating Acts::VolumeID to sPHENIX layer 
  std::map<unsigned int, unsigned int> base_layer_map = { {10, 0}, {12,3}, {14,7}, {16,55} };  
  int nsensors_stave[7] = {9,9,9,4,4,4,4};

  int getTpcRegion(int layer);
  void getGlobalLabels(Surface surf, int glbl_label[]);
  int getLabelBase(Acts::GeometryIdentifier id);
  void printBuffers(int index, Acts::Vector2 residual, 
		    Acts::Vector2 clus_sigma, float lcl_derivative[], 
		    float glbl_derivative[], int glbl_label[]);

  std::map<unsigned int, float> m_layerMisalignment;

  void setLayerUncInflation(const unsigned int layer, 
			    const float unc)
  {
    m_layerMisalignment.insert(std::make_pair(layer, unc));
  }

  static const int NGL = 6;
  static const int NLC = 5;
  siliconGrp si_grp = siliconGrp::snsr;
  tpcGrp tpc_grp = tpcGrp::htst;
  mmsGrp mms_grp = mmsGrp::tl;


}
#endif
