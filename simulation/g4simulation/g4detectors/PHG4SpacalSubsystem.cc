/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2014/08/12 03:49:12 $$
 */
#include "PHG4SpacalSubsystem.h"

#include "PHG4CylinderGeom_Spacalv1.h"  // for PHG4CylinderGeom_Spacalv1
#include "PHG4FullProjSpacalDetector.h"
#include "PHG4FullProjTiltedSpacalDetector.h"
#include "PHG4SpacalDetector.h"
#include "PHG4SpacalDisplayAction.h"
#include "PHG4SpacalSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <TSystem.h>

#include <iostream>  // for operator<<, basic_ostream
#include <sstream>

class PHG4Detector;

//_______________________________________________________________________
PHG4SpacalSubsystem::PHG4SpacalSubsystem(const std::string& na, const int lyr)
  : PHG4DetectorSubsystem(na, lyr)
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4SpacalSubsystem::~PHG4SpacalSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4SpacalSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  // create hit list only for active layers
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  // create detector
  //  _geom.set_layer( layer );
  //  if (lengthViaRapidityCoverage)
  //    {
  //       const double half_length =  PHG4Utils::GetLengthForRapidityCoverage(_geom.get_max_radius());
  //       _geom.set_zmin(-half_length);
  //       _geom.set_zmax(+half_length);
  //    }
  // create display settings before detector (detector adds its volumes to it)
  m_DisplayAction = new PHG4SpacalDisplayAction(Name());
  switch (GetParams()->get_int_param("config"))
  {
  case PHG4CylinderGeom_Spacalv1::kNonProjective:
    if (Verbosity() > 0) std::cout << "PHG4SpacalSubsystem::InitRun - use PHG4SpacalDetector" << std::endl;
    detector_ = new PHG4SpacalDetector(this, topNode, Name(), GetParams(), GetLayer());
    break;

  case PHG4CylinderGeom_Spacalv1::kFullProjective_2DTaper:
  case PHG4CylinderGeom_Spacalv1::kFullProjective_2DTaper_SameLengthFiberPerTower:
    if (Verbosity() > 0) std::cout << "PHG4SpacalSubsystem::InitRun - use PHG4FullProjSpacalDetector" << std::endl;
    detector_ = new PHG4FullProjSpacalDetector(this, topNode, Name(), GetParams(), GetLayer());
    break;

  case PHG4CylinderGeom_Spacalv1::kFullProjective_2DTaper_Tilted:
  case PHG4CylinderGeom_Spacalv1::kFullProjective_2DTaper_Tilted_SameLengthFiberPerTower:
    if (Verbosity() > 0) std::cout << "PHG4SpacalSubsystem::InitRun - use PHG4FullProjTiltedSpacalDetector" << std::endl;
    detector_ = new PHG4FullProjTiltedSpacalDetector(this, topNode, Name(), GetParams(), GetLayer());
    break;

  default:
    std::cout << "PHG4SpacalSubsystem::InitRun - unknown option exiting" << std::endl;
    gSystem->Exit(1);
    break;
  }

  detector_->SetActive(GetParams()->get_int_param("active"));
  detector_->SetAbsorberActive(GetParams()->get_int_param("absorberactive"));
  detector_->SuperDetector(SuperDetector());
  detector_->OverlapCheck(CheckOverlap());
  detector_->CosmicSetup(CosmicSetup());
  // the geometry object is set during detector construction, we need it for the
  // display to extract the visibility setting for logical volumes
  PHG4SpacalDisplayAction* DispAct = dynamic_cast<PHG4SpacalDisplayAction*>(m_DisplayAction);
  DispAct->SetGeom(detector_->get_geom());

  if (GetParams()->get_int_param("active"))
  {
    std::set<std::string> nodes;
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode* DetNode = dstNode;
    if (SuperDetector() != "NONE" && !SuperDetector().empty())
    {
      PHNodeIterator iter_dst(dstNode);
      DetNode = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode", SuperDetector()));
      if (!DetNode)
      {
        DetNode = new PHCompositeNode(SuperDetector());
        dstNode->addNode(DetNode);
      }
    }
    std::string detector_suffix = SuperDetector();
    if (detector_suffix == "NONE" || detector_suffix.empty())
    {
      detector_suffix = Name();
    }

    m_HitNodeName = "G4HIT_" + detector_suffix;
    nodes.insert(m_HitNodeName);
    m_AbsorberNodeName = "G4HIT_ABSORBER_" + detector_suffix;
    if (GetParams()->get_int_param("absorberactive"))
    {
      nodes.insert(m_AbsorberNodeName);
    }
    for (const auto& nodename : nodes)
    {
      PHG4HitContainer* g4_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(nodename);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, nodename, "PHObject"));
      }
      g4_hits->AddLayer(GetLayer());
    }

    steppingAction_ = new PHG4SpacalSteppingAction(detector_);
    steppingAction_->SetHitNodeName("G4HIT", m_HitNodeName);
    steppingAction_->SetHitNodeName("G4HIT_ABSORBER", m_AbsorberNodeName);
  }
  return 0;
}

//_______________________________________________________________________
int PHG4SpacalSubsystem::process_event(PHCompositeNode* topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
  {
    steppingAction_->SetInterfacePointers(topNode);
  }
  return 0;
}

//_______________________________________________________________________
PHG4Detector* PHG4SpacalSubsystem::GetDetector() const
{
  return detector_;
}

//_______________________________________________________________________
void PHG4SpacalSubsystem::Print(const std::string& what) const
{
  detector_->Print(what);
  return;
}

void PHG4SpacalSubsystem::SetDefaultParameters()
{
  set_default_double_param("xpos", 0.);  // translation in 3D
  set_default_double_param("ypos", 0.);  // translation in 3D
  set_default_double_param("zpos", 0.);  // translation in 3D

  set_default_double_param("thickness", 21.00000);
  set_default_double_param("radius", 90.);
  set_default_double_param("zmin", -149.470000);
  set_default_double_param("zmax", 149.470000);
  set_default_int_param("azimuthal_n_sec", 256);

  set_default_int_param("construction_verbose", 0.);
  set_default_int_param("azimuthal_seg_visible", 0.);
  set_default_int_param("virualize_fiber", 0.);
  set_default_int_param("config", static_cast<int>(PHG4CylinderGeom_Spacalv1::kNonProjective));

  set_default_double_param("divider_width", 0);       // radial size of the divider between blocks. <=0 means no dividers
  set_default_string_param("divider_mat", "G4_AIR");  // materials of the divider. G4_AIR is equivalent to not installing one in the term of material distribution

  return;
}
