#include "AlignmentDefs.h"
#include <trackbase/TpcDefs.h>

void AlignmentDefs::getSiliconGlobalLabels(Surface surf, int glbl_label[], AlignmentDefs::siliconGrp grp)
{
  Acts::GeometryIdentifier id = surf->geometryId();
  int group = 0;
  switch (grp)
  {
  case AlignmentDefs::siliconGrp::snsr:
    group = 0;
    break;
  case AlignmentDefs::siliconGrp::stv:
    group = 1;
    break;
  case AlignmentDefs::siliconGrp::brrl:
    group = 2;
    break;
  }

  int label_base = getLabelBase(id, 0, group);
  for (int i = 0; i < NGL; ++i)
  {
    glbl_label[i] = label_base + i;
  }
}

void AlignmentDefs::getTpcGlobalLabels(Surface surf, TrkrDefs::cluskey cluskey, int glbl_label[], AlignmentDefs::tpcGrp grp)
{
  Acts::GeometryIdentifier id = surf->geometryId();
  int group = 0;
  switch (grp)
  {
  case AlignmentDefs::tpcGrp::htst:
    group = 3;
    break;
  case AlignmentDefs::tpcGrp::sctr:
    group = 4;
    break;
  case AlignmentDefs::tpcGrp::tp:
    group = 5;
    break;
  }

  int label_base = getLabelBase(id, cluskey, group);
  for (int i = 0; i < NGL; ++i)
  {
    glbl_label[i] = label_base + i;
  }
}
void AlignmentDefs::getMMGlobalLabels(Surface surf, int glbl_label[], AlignmentDefs::mmsGrp grp)
{
  Acts::GeometryIdentifier id = surf->geometryId();
  int group = 0;
  switch (grp)
  {
  case AlignmentDefs::mmsGrp::tl:
    group = 6;
    break;
  case AlignmentDefs::mmsGrp::mm:
    group = 7;
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
    region = 1;
  if (layer > 38 && layer < 55)
    region = 2;

  return region;
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
      // layer only, assign all sensors to sensor 0
      label_base += layer * 1000000 + 0;
    return label_base;
  }
  else if(layer > 2 && layer < 7)
    {
      // calculating the stave number from the sensor is different between the INTT and MVTX
      // There are 4 sensors/stave, but they are mapped to staves in a strange way
      // sensors 1-> (nstaves/layer)*2 are in staves 1->nstaves/layer in pairs
      // sensors (nstaves/layer * 2) +1->nstaves/layer)*2 are in staves 1->nstaves/layer in pairs

      int stave;
      unsigned int breakat = nstaves_layer_intt[layer-3] * 2;
      if(sensor < breakat)
	{
	  stave = sensor/2;  // staves 0 -> (nstaves/layer) -1
	}
      else
	{
	  stave = (sensor - breakat)/2;  //staves 0 -> (nstaves/layer) -1
	}

      if (group == 0)
	{
	  // every sensor has a different label
	  label_base += layer * 1000000 + stave * 10000 + sensor * 10;
	  return label_base;
	}
      if (group == 1)
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
      if (group == 2)
	// layer only, assign all sensors to sensor 0
	label_base += layer * 1000000 + 0;
      return label_base;
    }
  else if (layer > 6 && layer < 55)
  {
    if (group == 3)
      {
	// want every hitset (layer, sector, side) to have a separate label
      // each group of 12 subsurfaces (sensors) is in a single hitset
      int hitset = sensor / 12;  // 0-11 on side 0, 12-23 on side 1
      label_base += layer * 1000000 + hitset * 10000;
      return label_base;
    }
    if (group == 4)
    {
      // group all tpc layers in each region and sector, assign layer 7 and side and sector number to all layers and hitsets
      int side = TpcDefs::getSide(cluskey);
      int sector = TpcDefs::getSectorId(cluskey);;
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
    if (group == 5)
    {
      // all tpc layers and all sectors, assign layer 7 and sensor 0 to all layers and sensors
      label_base += 7 * 1000000 + 0;
      return label_base;
    }
  }
  else
  {
    if (group == 6)
    {
      // every tile has different label
      int tile = sensor;
      label_base += layer * 1000000 + tile * 10000 + sensor * 10;
      return label_base;
    }
    if (group == 7)
    {
      // assign layer 55 and tile 0 to all
      label_base += 55 * 1000000 + 0;
      return label_base;
    }
  }
  return -1;
}

void AlignmentDefs::printBuffers(int index, Acts::Vector2 residual, Acts::Vector2 clus_sigma, float lcl_derivative[], float glbl_derivative[], int glbl_label[])
{
  std::cout << " float buffer: "
            << " residual "
            << "  " << residual(index);
  for (int il = 0; il < NLC; ++il)
  {
    if (lcl_derivative[il] != 0) std::cout << " lcl_deriv[" << il << "] " << lcl_derivative[il] << "  ";
  }
  std::cout << " sigma "
            << "  " << clus_sigma(index) << "  ";
  for (int ig = 0; ig < NGL; ++ig)
  {
    if (glbl_derivative[ig] != 0) std::cout << " glbl_deriv[" << ig << "] " << glbl_derivative[ig] << "  ";
  }
  std::cout << " int buffer: "
            << " 0 "
            << " 0 "
            << " ";  // spacer, rmeas placeholder
  for (int il = 0; il < NLC; ++il)
  {
    if (lcl_derivative[il] != 0) std::cout << " lcl_label[" << il << "] " << il + 1 << "  ";
  }
  std::cout << " 0 "
            << "  ";
  for (int ig = 0; ig < NGL; ++ig)
  {
    if (glbl_derivative[ig] != 0) std::cout << " glbl_label[" << ig << "] " << glbl_label[ig] << "  ";
  }
  std::cout << " end of meas " << std::endl;
}
