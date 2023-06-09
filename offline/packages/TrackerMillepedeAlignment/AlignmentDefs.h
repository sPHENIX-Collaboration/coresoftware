#ifndef ALIGNMENTDEFS_H
#define ALIGNMENTDEFS_H

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

namespace AlignmentDefs
{
  enum siliconGrp
  {
    snsr,
    stv,
    brrl
  };
  enum tpcGrp
  {
    htst,
    sctr,
    tp
  };
  enum mmsGrp
  {
    tl,
    mm
  };

  static constexpr int NGL = 6;
  static constexpr int NLC = 5;

  //! Map relating Acts::VolumeID to sPHENIX layer
  static const std::map<unsigned int, unsigned int> base_layer_map = {{10, 0}, {12, 3}, {14, 7}, {16, 55}};
  static constexpr int nsensors_stave[7] = {9, 9, 9, 4, 4, 4, 4};
  static constexpr int nstaves_layer_intt[7] = {12, 12, 16, 16};

  int getTpcRegion(int layer);
  void getSiliconGlobalLabels(Surface surf, int glbl_label[], siliconGrp grp);
  void getTpcGlobalLabels(Surface surf, TrkrDefs::cluskey cluskey, int glbl_label[], tpcGrp grp);
  void getMMGlobalLabels(Surface surf, int glbl_label[], mmsGrp grp);

  int getLabelBase(Acts::GeometryIdentifier id, TrkrDefs::cluskey cluskey, int group);

  void printBuffers(int index, Acts::Vector2 residual,
                    Acts::Vector2 clus_sigma, float lcl_derivative[],
                    float glbl_derivative[], int glbl_label[]);

}  // namespace AlignmentDefs
#endif
