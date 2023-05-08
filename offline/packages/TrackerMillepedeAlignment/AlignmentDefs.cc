#include "AlignmentDefs.h"


void AlignmentDefs::getGlobalLabels(Surface surf, int glbl_label[])
{
  Acts::GeometryIdentifier id = surf->geometryId();
  int label_base = getLabelBase(id);   // This value depends on how the surfaces are grouped
  for(int i=0; i<NGL; i++)
    {
      glbl_label[i] = label_base + i;
    }

}
int AlignmentDefs::getTpcRegion(int layer)
{
  int region = 0;
  if(layer > 23 && layer < 39)
    region = 1;
  if(layer > 38 && layer < 55)
    region = 2;

  return region;  
}
int AlignmentDefs::getLabelBase(Acts::GeometryIdentifier id)
{
  
  unsigned int volume = id.volume(); 
  unsigned int acts_layer = id.layer();
  unsigned int layer = base_layer_map.find(volume)->second + acts_layer / 2 -1;
  unsigned int sensor = id.sensitive() - 1;  // Acts starts at 1

  int label_base = 1;  // Mille wants to start at 1

  // decide what level of grouping we want
  if(layer < 7)
    {
      if(si_grp == siliconGrp::snsr)
	{
	  // every sensor has a different label
	  int stave = sensor / nsensors_stave[layer];
	  label_base += layer*1000000  + stave*10000 + sensor*10;
	  return label_base;
	}
      if(si_grp == siliconGrp::stv)
	{
	  // layer and stave, assign all sensors to the stave number
	  int stave = sensor / nsensors_stave[layer];
	  label_base += layer*1000000 + stave*10000;
	  return label_base;
	}
      if(si_grp == siliconGrp::brrl)
	// layer only, assign all sensors to sensor 0 
	label_base += layer*1000000 + 0;
      return label_base;
    }
  else if(layer > 6 && layer < 55)
    {
      if(tpc_grp == tpcGrp::htst)
	{
	  // want every hitset (layer, sector, side) to have a separate label
	  // each group of 12 subsurfaces (sensors) is in a single hitset
	  int hitset = sensor/12; // 0-11 on side 0, 12-23 on side 1
	  label_base += layer*1000000 + hitset*10000;
	  return label_base;
	}
      if(tpc_grp == tpcGrp::sctr)
	{
	  // group all tpc layers in each region and sector, assign layer 7 and side and sector number to all layers and hitsets
	  int side = sensor / 144; // 0-143 on side 0, 144-287 on side 1
	  int sector = (sensor - side *144) / 12; 
	  // for a given layer there are only 12 sectors x 2 sides
	  // The following gives the sectors in the inner, mid, outer regions unique group labels
	  int region = getTpcRegion(layer);  // inner, mid, outer
	  label_base += 7*1000000 + (region * 24 + side*12 + sector) *10000; 
	  return label_base;
	}
      if(tpc_grp == tpcGrp::tp)
	{
	  // all tpc layers and all sectors, assign layer 7 and sensor 0 to all layers and sensors
	  label_base += 7*1000000 + 0;
	  return label_base;
	}
    }
  else
    {
      if(mms_grp == mmsGrp::tl)
	{
	  // every tile has different label
	  int tile = sensor;
	  label_base += layer*1000000 + tile*10000+sensor*10;
	  return label_base;
	}
      if(mms_grp == mmsGrp::mm)
	{
	  // assign layer 55 and tile 0 to all
	  label_base += 55*1000000 + 0;	  
	  return label_base;
	}
    }

  return -1;
}

void AlignmentDefs::printBuffers(int index, Acts::Vector2 residual, Acts::Vector2 clus_sigma, float lcl_derivative[], float glbl_derivative[], int glbl_label[])
{
    std::cout << " float buffer: " << " residual " << "  " << residual(index);
  for (int il=0;il<NLC;++il) { if(lcl_derivative[il] != 0) std::cout << " lcl_deriv["<< il << "] " << lcl_derivative[il] << "  ";  }
  std::cout  << " sigma " << "  " << clus_sigma(index) << "  ";
  for (int ig=0;ig<NGL;++ig) { if(glbl_derivative[ig] != 0)  std::cout << " glbl_deriv["<< ig << "] " << glbl_derivative[ig] << "  ";  }
  std::cout << " int buffer: " << " 0 " << " 0 " << " ";  // spacer, rmeas placeholder
  for (int il=0;il<NLC;++il) { if(lcl_derivative[il] != 0) std::cout << " lcl_label["<< il << "] " << il+1 << "  ";  }
  std::cout << " 0 " << "  ";
  for (int ig=0;ig<NGL;++ig) { if(glbl_derivative[ig] != 0) std::cout << " glbl_label["<< ig << "] " << glbl_label[ig] << "  ";  }
  std::cout << " end of meas " << std::endl;
}

