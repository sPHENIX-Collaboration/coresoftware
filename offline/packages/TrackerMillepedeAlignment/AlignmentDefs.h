#ifndef ALIGNMENTDEFS_H
#define ALIGNMENTDEFS_H

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

namespace AlignmentDefs
{
  enum mvtxGrp
  {
    snsr,
    stv,
    mvtxlyr,
    clamshl
  };

  enum inttGrp
  {
    chp,
    lad,
    inttlyr,
    inttbrl
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

  static constexpr int NGLVTX = 3;
  static constexpr int NGL = 6;
  static constexpr int NLC = 5;

  //! Map relating Acts::VolumeID to sPHENIX layer
  static const std::map<unsigned int, unsigned int> base_layer_map = {{10, 0}, {12, 3}, {14, 7}, {16, 55}};
  static constexpr int nsensors_stave[7] = {9, 9, 9, 4, 4, 4, 4};
  static constexpr int nstaves_layer_intt[4] = {12, 12, 16, 16};
  static constexpr int nstaves_layer_mvtx[3] = {12, 16, 20};

  static constexpr int clamshell_stave_list[3][2][10] = {  // [layer][clamshell][stave]
      {{3, 4, 5, 6, 7, 8, 0, 0, 0, 0},
       {0, 1, 2, 9, 10, 11, 0, 0, 0, 0}},
      {{4, 5, 6, 7, 8, 9, 10, 11, 0, 0},
       {0, 1, 2, 3, 12, 13, 14, 15, 0, 0}},
      {{5, 6, 7, 8, 9, 10, 11, 12, 13, 14},
       {0, 1, 2, 3, 4, 15, 16, 17, 18, 19}}};

  static constexpr int glbl_vtx_label[NGLVTX] = {60000000, 60000001, 60000002};

  int getTpcRegion(int layer);
  void getMvtxGlobalLabels(const Surface& surf, int glbl_label[], mvtxGrp grp);
  void getInttGlobalLabels(const Surface& surf, int glbl_label[], inttGrp grp);
  void getTpcGlobalLabels(const Surface& surf, TrkrDefs::cluskey cluskey, int glbl_label[], tpcGrp grp);
  void getMMGlobalLabels(const Surface& surf, int glbl_label[], mmsGrp grp);
  int getMvtxClamshell(int layer, int stave);

  std::vector<int> getAllMvtxGlobalLabels(int grp);
  std::vector<int> getAllInttGlobalLabels(int grp);
  std::vector<int> getAllTpcGlobalLabels(int grp);
  std::vector<int> makeLabelsFromBase(std::vector<int>& label_base);

  int getLabelBase(Acts::GeometryIdentifier id, TrkrDefs::cluskey cluskey, int group);

  void printBuffers(int index, Acts::Vector2 residual,
                    Acts::Vector2 clus_sigma, float lcl_derivative[],
                    float glbl_derivative[], int glbl_label[]);

}  // namespace AlignmentDefs
#endif
