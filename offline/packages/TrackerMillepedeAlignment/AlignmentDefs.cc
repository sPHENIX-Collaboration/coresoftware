#include "AlignmentDefs.h"
#include <trackbase/TpcDefs.h>

void AlignmentDefs::getMvtxGlobalLabels(const Surface& surf, int glbl_label[], AlignmentDefs::mvtxGrp grp)
{
  Acts::GeometryIdentifier id = surf->geometryId();
  int group = 0;
  switch (grp)
  {
  case AlignmentDefs::mvtxGrp::snsr:
    group = 0;
    break;
  case AlignmentDefs::mvtxGrp::stv:
    group = 1;
    break;
  case AlignmentDefs::mvtxGrp::mvtxlyr:
    group = 2;
    break;

  case AlignmentDefs::mvtxGrp::clamshl:
    group = 3;
    break;
  }

  int label_base = getLabelBase(id, 0, group);
  for (int i = 0; i < NGL; ++i)
  {
    glbl_label[i] = label_base + i;
  }
}

void AlignmentDefs::getInttGlobalLabels(const Surface& surf, int glbl_label[], AlignmentDefs::inttGrp grp)
{
  Acts::GeometryIdentifier id = surf->geometryId();
  int group = 0;
  switch (grp)
  {
  case AlignmentDefs::inttGrp::chp:
    group = 4;
    break;
  case AlignmentDefs::inttGrp::lad:
    group = 5;
    break;
  case AlignmentDefs::inttGrp::inttlyr:
    group = 6;
    break;
  case AlignmentDefs::inttGrp::inttbrl:
    group = 7;
    break;
  }

  int label_base = getLabelBase(id, 0, group);
  for (int i = 0; i < NGL; ++i)
  {
    glbl_label[i] = label_base + i;
  }
}

void AlignmentDefs::getTpcGlobalLabels(const Surface& surf, TrkrDefs::cluskey cluskey, int glbl_label[], AlignmentDefs::tpcGrp grp)
{
  Acts::GeometryIdentifier id = surf->geometryId();
  int group = 0;
  switch (grp)
  {
  case AlignmentDefs::tpcGrp::htst:
    group = 8;
    break;
  case AlignmentDefs::tpcGrp::sctr:
    group = 9;
    break;
  case AlignmentDefs::tpcGrp::tp:
    group = 10;
    break;
  }

  int label_base = getLabelBase(id, cluskey, group);
  for (int i = 0; i < NGL; ++i)
  {
    glbl_label[i] = label_base + i;
  }
}
void AlignmentDefs::getMMGlobalLabels(const Surface& surf, int glbl_label[], AlignmentDefs::mmsGrp grp)
{
  Acts::GeometryIdentifier id = surf->geometryId();
  int group = 0;
  switch (grp)
  {
  case AlignmentDefs::mmsGrp::tl:
    group = 11;
    break;
  case AlignmentDefs::mmsGrp::mm:
    group = 12;
    break;
  }

  int label_base = getLabelBase(id, 0, group);
  for (int i = 0; i < NGL; ++i)
  {
    glbl_label[i] = label_base + i;
  }
}

int AlignmentDefs::getTpcRegion(int layer)
{
  int region = 0;
  if (layer > 22 && layer < 39)
  {
    region = 1;
  }
  if (layer > 38 && layer < 55)
  {
    region = 2;
  }

  return region;
}

int AlignmentDefs::getMvtxClamshell(int layer, int stave)
{
  for (int istave = 0; istave < nstaves_layer_mvtx[layer]; ++istave)
  {
    for (int ishell = 0; ishell < 2; ++ishell)
    {
      int stave_ref = clamshell_stave_list[layer][ishell][istave];
      if (stave == stave_ref)
      {
        return ishell;
      }
    }
  }

  std::cout << " AlignemntDefs::getMvtxClamshell: did not find stave " << stave << std::endl;
  return 0;
}

int AlignmentDefs::getLabelBase(Acts::GeometryIdentifier id, TrkrDefs::cluskey cluskey, int group)
{
  unsigned int volume = id.volume();
  unsigned int acts_layer = id.layer();
  unsigned int layer = base_layer_map.find(volume)->second + acts_layer / 2 - 1;
  unsigned int sensor = id.sensitive() - 1;  // Acts starts at 1

  int label_base = 1;  // Mille wants to start at 1

  // decide what level of grouping we want
  if (layer < 3)
  {
    if (group == 0)
    {
      // every sensor has a different label
      int stave = sensor / nsensors_stave[layer];
      label_base += layer * 1000000 + stave * 10000 + sensor * 10;
      return label_base;
    }
    if (group == 1)
    {
      // layer and stave, assign all sensors to the stave number
      int stave = sensor / nsensors_stave[layer];
      label_base += layer * 1000000 + stave * 10000;
      /*
        std::cout << id << std::endl;
        std::cout << "    label_base " << label_base << " volume " << volume << " acts_layer " << acts_layer
        << " layer " << layer << " stave " << stave << " sensor " << sensor << std::endl;
      */
      return label_base;
    }
    if (group == 2)
    {
      // layer only, assign all sensors to sensor 0 in each clamshell
      int stave = sensor / nsensors_stave[layer];
      int clamshell = getMvtxClamshell(layer, stave);
      label_base += layer * 1000000 + clamshell * 10000;
      //	std::cout << " mvtx group 2 layer " << layer << " sensor " << sensor << " stave " << stave
      //	  << " clamshell " << clamshell << " label_base " << label_base << std::endl;

      return label_base;
    }
    if (group == 3)
    {
      // group by half-barrel, or clamshell
      // Assume for now low staves are in clamshell 0 - check!!!
      int stave = sensor / nsensors_stave[layer];
      int breakat = nstaves_layer_mvtx[layer] / 2;
      int clamshell = 1;
      if (stave < breakat)
      {
        clamshell = 0;
      }
      label_base += 0 * 1000000 + clamshell * 10000;
      //	std::cout << " mvtx group 3 layer " << layer << " sensor " << sensor << " clamshell " << clamshell << " label_base " << label_base << std::endl;

      return label_base;
    }
  }
  else if (layer > 2 && layer < 7)
  {
    // calculating the stave number from the sensor is different between the INTT and MVTX
    // There are 4 sensors/stave, but they are mapped to staves in a strange way
    // sensors 1-> (nstaves/layer)*2 are in staves 1->nstaves/layer in pairs
    // sensors (nstaves/layer * 2) +1->nstaves/layer)*2 are in staves 1->nstaves/layer in pairs

    int stave;
    unsigned int breakat = nstaves_layer_intt[layer - 3] * 2;
    if (sensor < breakat)
    {
      stave = sensor / 2;  // staves 0 -> (nstaves/layer) -1
    }
    else
    {
      stave = (sensor - breakat) / 2;  // staves 0 -> (nstaves/layer) -1
    }

    if (group == 4)
    {
      // every sensor has a different label
      label_base += layer * 1000000 + stave * 10000 + sensor * 10;
      return label_base;
    }
    if (group == 5)
    {
      // layer and stave, assign all sensors to the stave number
      label_base += layer * 1000000 + stave * 10000;
      /*
      std::cout << "    "  << id << std::endl;
      std::cout << "    label_base " << label_base << " volume " << volume << " acts_layer " << acts_layer
                << " layer " << layer << " breakat " << breakat << " stave " << stave << " sensor " << sensor << std::endl;
      */
      return label_base;
    }
    if (group == 6)
    {
      // layer only, assign all sensors to sensor 0 for this layer
      label_base += layer * 1000000 + 0;
      return label_base;
    }
    if (group == 7)
    {
      // entire INTT
      // assign all sensors to layer 3
      label_base += 3 * 1000000 + 0;
      return label_base;
    }
  }
  else if (layer > 6 && layer < 55)
  {
    if (group == 8)
    {
      // want every hitset (layer, sector, side) to have a separate label
      // each group of 12 subsurfaces (sensors) is in a single hitset
      int hitset = sensor / 12;  // 0-11 on side 0, 12-23 on side 1
      label_base += layer * 1000000 + hitset * 10000;
      return label_base;
    }
    if (group == 9)
    {
      // group all tpc layers in each region and sector, assign layer 7 and side and sector number to all layers and hitsets
      int side = TpcDefs::getSide(cluskey);
      int sector = TpcDefs::getSectorId(cluskey);
      ;
      int region = getTpcRegion(layer);  // inner, mid, outer

      // for a given layer there are only 12 sectors x 2 sides
      // The following gives the sectors in the inner, mid, outer regions unique group labels
      label_base += 7 * 1000000 + (region * 24 + side * 12 + sector) * 10000;
      /*
      std::cout << " sensor " << sensor << " sector " << sector << " region " << region
                << " side " << side << " label base " << label_base << std::endl;
      std::cout << " Volume " << volume << " acts_layer " << acts_layer << " base_layer " <<  base_layer_map.find(volume)->second << " layer " << layer << std::endl;
      */

      return label_base;
    }
    if (group == 10)
    {
      // all tpc layers and all sectors, assign layer 7 and sensor 0 to all layers and sensors
      label_base += 7 * 1000000 + 0;
      return label_base;
    }
  }
  else
  {
    if (group == 11)
    {
      // every tile has different label
      int tile = sensor;
      label_base += layer * 1000000 + tile * 10000 + sensor * 10;
      return label_base;
    }
    if (group == 12)
    {
      // assign layer 55 and tile 0 to all
      label_base += 55 * 1000000 + 0;
      return label_base;
    }
  }
  return -1;
}

std::vector<int> AlignmentDefs::getAllMvtxGlobalLabels(int grp)
{
  std::vector<int> label_base;

  if (grp == mvtxGrp::clamshl)
  {
    for (int ishl = 0; ishl < 2; ++ishl)
    {
      int label = 1 + 0 * 1000000 + ishl * 10000;
      label_base.push_back(label);
    }
  }
  else if (grp == mvtxGrp::mvtxlyr)
  {
    for (int ishl = 0; ishl < 2; ++ishl)
    {
      for (int ilyr = 0; ilyr < 3; ++ilyr)
      {
        int label = 1 + ilyr * 1000000 + ishl * 10000;
        label_base.push_back(label);
      }
    }
  }
  else if (grp == mvtxGrp::stv)
  {
    for (int ilyr = 0; ilyr < 3; ++ilyr)
    {
      for (int istv = 0; istv < nstaves_layer_mvtx[ilyr]; ++istv)
      {
        int label = 1 + ilyr * 1000000 + istv * 10000;
        label_base.push_back(label);
      }
    }
  }
  else if (grp == mvtxGrp::snsr)
  {
    for (int ilyr = 0; ilyr < 3; ++ilyr)
    {
      for (int istv = 0; istv < nstaves_layer_mvtx[ilyr]; ++istv)
      {
        for (int isnsr = 0; isnsr < nsensors_stave[ilyr]; ++isnsr)
        {
          int label = 1 + ilyr * 1000000 + istv * 10000 + isnsr * 10;
          label_base.push_back(label);
        }
      }
    }
  }

  auto labels = makeLabelsFromBase(label_base);

  return labels;
}

std::vector<int> AlignmentDefs::getAllInttGlobalLabels(int grp)
{
  std::vector<int> label_base;

  if (grp == inttGrp::inttbrl)
  {
    int label = 1 + 3 * 1000000 + 0;
    label_base.push_back(label);
  }
  else if (grp == inttGrp::inttlyr)
  {
    for (int ilyr = 3; ilyr < 7; ++ilyr)
    {
      int label = 1 + ilyr * 1000000 + 0;
      label_base.push_back(label);
    }
  }
  else if (grp == inttGrp::lad)
  {
    for (int ilyr = 3; ilyr < 7; ++ilyr)
    {
      for (int istv = 0; istv < nstaves_layer_intt[ilyr - 3]; ++istv)
      {
        int label = 1 + ilyr * 1000000 + istv * 10000;
        label_base.push_back(label);
      }
    }
  }
  else if (grp == inttGrp::chp)
  {
    for (int ilyr = 3; ilyr < 7; ++ilyr)
    {
      for (int istv = 0; istv < nstaves_layer_intt[ilyr - 3]; ++istv)
      {
        for (int isnsr = 0; isnsr < nsensors_stave[ilyr]; ++isnsr)
        {
          int label = 1 + ilyr * 1000000 + istv * 10000 + isnsr * 10;
          label_base.push_back(label);
        }
      }
    }
  }

  auto labels = makeLabelsFromBase(label_base);

  return labels;
}

std::vector<int> AlignmentDefs::getAllTpcGlobalLabels(int grp)
{
  std::vector<int> label_base;

  if (grp == tpcGrp::sctr)
  {
    for (int isec = 0; isec < 72; ++isec)
    {
      int label = 1 + 7 * 1000000 + isec * 10000;
      label_base.push_back(label);
    }
  }

  auto labels = makeLabelsFromBase(label_base);

  return labels;
}

std::vector<int> AlignmentDefs::makeLabelsFromBase(std::vector<int>& label_base)
{
  std::vector<int> labels;
  for (int ilbl : label_base)
  {
    for (int ipar = 0; ipar < 6; ++ipar)
    {
      int label_plus = ilbl + ipar;
      labels.push_back(label_plus);
    }
  }

  return labels;
}

void AlignmentDefs::printBuffers(int index, Acts::Vector2 residual, Acts::Vector2 clus_sigma, float lcl_derivative[], float glbl_derivative[], int glbl_label[])
{
  std::cout << " float buffer: "
            << " residual "
            << "  " << residual(index);
  for (int il = 0; il < NLC; ++il)
  {
    if (lcl_derivative[il] != 0)
    {
      std::cout << " lcl_deriv[" << il << "] " << lcl_derivative[il] << "  ";
    }
  }
  std::cout << " sigma "
            << "  " << clus_sigma(index) << "  ";
  for (int ig = 0; ig < NGL; ++ig)
  {
    if (glbl_derivative[ig] != 0)
    {
      std::cout << " glbl_deriv[" << ig << "] " << glbl_derivative[ig] << "  ";
    }
  }
  std::cout << " int buffer: "
            << " 0 "
            << " 0 "
            << " ";  // spacer, rmeas placeholder
  for (int il = 0; il < NLC; ++il)
  {
    if (lcl_derivative[il] != 0)
    {
      std::cout << " lcl_label[" << il << "] " << il + 1 << "  ";
    }
  }
  std::cout << " 0 "
            << "  ";
  for (int ig = 0; ig < NGL; ++ig)
  {
    if (glbl_derivative[ig] != 0)
    {
      std::cout << " glbl_label[" << ig << "] " << glbl_label[ig] << "  ";
    }
  }
  std::cout << " end of meas " << std::endl;
}
