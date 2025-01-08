#include "DumpPHG4CylinderGeomContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom_Spacalv1.h>  // for PHG4CylinderGeom_...
#include <g4detectors/PHG4CylinderGeom_Spacalv3.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<PHG4CylinderGeomContainer>;

DumpPHG4CylinderGeomContainer::DumpPHG4CylinderGeomContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHG4CylinderGeomContainer::process_Node(PHNode *myNode)
{
  PHG4CylinderGeomContainer *phg4geomcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    phg4geomcontainer = thisNode->getData();
  }
  if (phg4geomcontainer)
  {
    PHG4CylinderGeomContainer::ConstIterator hiter;
    PHG4CylinderGeomContainer::ConstRange geom_begin_end = phg4geomcontainer->get_begin_end();
    *fout << "num layers: " << phg4geomcontainer->get_NLayers() << std::endl;
    for (hiter = geom_begin_end.first; hiter != geom_begin_end.second; hiter++)
    {
      *fout << "layer: " << hiter->second->get_layer() << std::endl;
      *fout << "radius: " << hiter->second->get_radius() << std::endl;
      *fout << "thickness: " << hiter->second->get_thickness() << std::endl;
      *fout << "zmin: " << hiter->second->get_zmin() << std::endl;
      *fout << "zmax: " << hiter->second->get_zmax() << std::endl;
      *fout << "nscint: " << hiter->second->get_nscint() << std::endl;
      *fout << "tiltangle: " << hiter->second->get_tiltangle() << std::endl;
      *fout << "strip_y_spacing: " << hiter->second->get_strip_y_spacing() << std::endl;
      *fout << "strip_z_spacing: " << hiter->second->get_strip_z_spacing(0) << std::endl;
      *fout << "strip_tilt: " << hiter->second->get_strip_tilt() << std::endl;
      *fout << "N_strip_columns: " << hiter->second->get_N_strip_columns() << std::endl;
      *fout << "N_strips_per_column: " << hiter->second->get_N_strips_per_column() << std::endl;
      *fout << "N_sensors_in_layer: " << hiter->second->get_N_sensors_in_layer() << std::endl;
      *fout << "pixel_z: " << hiter->second->get_pixel_z() << std::endl;
      *fout << "pixel_x: " << hiter->second->get_pixel_x() << std::endl;
      *fout << "pixel_thickness: " << hiter->second->get_pixel_thickness() << std::endl;
      PHG4CylinderGeom_Spacalv1 *layergeomv1 = dynamic_cast<PHG4CylinderGeom_Spacalv1 *>(hiter->second);
      if (layergeomv1)
      {
        const PHG4CylinderGeom_Spacalv3::sector_map_t &sector_map = layergeomv1->get_sector_map();
        *fout << "xpos: " << layergeomv1->get_xpos() << std::endl;
        *fout << "ypos: " << layergeomv1->get_ypos() << std::endl;
        *fout << "zpos: " << layergeomv1->get_zpos() << std::endl;
        *fout << "fiber_clading_thickness: " << layergeomv1->get_fiber_clading_thickness() << std::endl;
        *fout << "fiber_core_diameter: " << layergeomv1->get_fiber_core_diameter() << std::endl;
        *fout << "fiber_distance: " << layergeomv1->get_fiber_distance() << std::endl;
        *fout << "absorber_mat: " << layergeomv1->get_absorber_mat() << std::endl;
        *fout << "fiber_clading_mat: " << layergeomv1->get_fiber_clading_mat() << std::endl;
        *fout << "fiber_core_mat: " << layergeomv1->get_fiber_core_mat() << std::endl;
        *fout << "virualize_fiber: " << layergeomv1->is_virualize_fiber() << std::endl;
        for (auto sectormapiter : sector_map)
        {
          *fout << "sector " << sectormapiter.first << ", rotation: " << sectormapiter.second << std::endl;
        }
      }
      PHG4CylinderGeom_Spacalv3 *layergeomv3 = dynamic_cast<PHG4CylinderGeom_Spacalv3 *>(hiter->second);
      if (layergeomv3)
      {
        *fout << "sidewall_outer_torr: " << layergeomv3->get_sidewall_outer_torr() << std::endl;
        *fout << "sidewall_thickness: " << layergeomv3->get_sidewall_thickness() << std::endl;
        *fout << "sidewall_mat: " << layergeomv3->get_sidewall_mat() << std::endl;
        *fout << "max_phi_bin_in_sec: " << layergeomv3->get_max_phi_bin_in_sec() << std::endl;
        *fout << "divider_mat: " << layergeomv3->get_divider_mat() << std::endl;
        *fout << "divider_width: " << layergeomv3->get_divider_width() << std::endl;

        const PHG4CylinderGeom_Spacalv3::tower_map_t &tower_map = layergeomv3->get_sector_tower_map();
        for (const auto &towermapiter : tower_map)
        {
          *fout << "tower " << towermapiter.first << ", id: " << towermapiter.second.id << std::endl;
          *fout << "tower " << towermapiter.first << ", pDz: " << towermapiter.second.pDz << std::endl;
          *fout << "tower " << towermapiter.first << ", pDy1: " << towermapiter.second.pDy1 << std::endl;
          *fout << "tower " << towermapiter.first << ", pDx1: " << towermapiter.second.pDx1 << std::endl;
          *fout << "tower " << towermapiter.first << ", pDx2: " << towermapiter.second.pDx2 << std::endl;
          *fout << "tower " << towermapiter.first << ", pDy2: " << towermapiter.second.pDy2 << std::endl;
          *fout << "tower " << towermapiter.first << ", pDx3: " << towermapiter.second.pDx3 << std::endl;
          *fout << "tower " << towermapiter.first << ", pDx4: " << towermapiter.second.pDx4 << std::endl;
          *fout << "tower " << towermapiter.first << ", pTheta: " << towermapiter.second.pTheta << std::endl;
          *fout << "tower " << towermapiter.first << ", pPhi: " << towermapiter.second.pPhi << std::endl;
          *fout << "tower " << towermapiter.first << ", pAlp1: " << towermapiter.second.pAlp1 << std::endl;
          *fout << "tower " << towermapiter.first << ", pAlp2: " << towermapiter.second.pAlp2 << std::endl;
          *fout << "tower " << towermapiter.first << ", pRotationAngleX: " << towermapiter.second.pRotationAngleX << std::endl;
          *fout << "tower " << towermapiter.first << ", centralX: " << towermapiter.second.centralX << std::endl;
          *fout << "tower " << towermapiter.first << ", centralY: " << towermapiter.second.centralY << std::endl;
          *fout << "tower " << towermapiter.first << ", centralZ: " << towermapiter.second.centralZ << std::endl;
          *fout << "tower " << towermapiter.first << ", ModuleSkinThickness: " << towermapiter.second.ModuleSkinThickness << std::endl;
          *fout << "tower " << towermapiter.first << ", NFiberX: " << towermapiter.second.NFiberX << std::endl;
          *fout << "tower " << towermapiter.first << ", NFiberY: " << towermapiter.second.NFiberY << std::endl;
          *fout << "tower " << towermapiter.first << ", NSubtowerX: " << towermapiter.second.NSubtowerX << std::endl;
          *fout << "tower " << towermapiter.first << ", NSubtowerY: " << towermapiter.second.NSubtowerY << std::endl;
          *fout << "tower " << towermapiter.first << ", LightguideHeight: " << towermapiter.second.LightguideHeight << std::endl;
          *fout << "tower " << towermapiter.first << ", LightguideTaperRatio: " << towermapiter.second.LightguideTaperRatio << std::endl;
          *fout << "tower " << towermapiter.first << ", LightguideMaterial: " << towermapiter.second.LightguideMaterial << std::endl;
        }
      }
    }
  }
  return 0;
}
