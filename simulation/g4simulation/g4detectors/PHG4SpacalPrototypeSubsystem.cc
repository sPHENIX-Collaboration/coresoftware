// $$Id: PHG4SpacalPrototypeSubsystem.cc,v 1.2 2014/08/12 03:49:12 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2014/08/12 03:49:12 $$
 */
#include "PHG4SpacalPrototypeSubsystem.h"

#include "PHG4SpacalPrototypeDetector.h"
#include "PHG4ProjSpacalDetector.h"
#include "PHG4FullProjSpacalDetector.h"
#include "PHG4CylinderGeom.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4SpacalPrototypeSteppingAction.h"
#include "PHG4EventActionClearZeroEdep.h"

#include <g4main/PHG4Utils.h>
#include <g4main/PHG4PhenixDetector.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <pdbcalbase/PdbParameterMap.h>
#include <Geant4/globals.hh>

#include <sstream>
#include <cassert>

using namespace std;

//_______________________________________________________________________
PHG4SpacalPrototypeSubsystem::PHG4SpacalPrototypeSubsystem(
    const std::string &na) :
    detector_(NULL), steppingAction_(NULL), eventAction_(NULL), //
    active(0), absorberactive(0), //
    detector_type(na), superdetector("NONE"),//
    useDB(false),Params(na)
{
  Name(na);
}

//_______________________________________________________________________
int
PHG4SpacalPrototypeSubsystem::InitRun(PHCompositeNode* topNode)
{
  // create hit list only for active layers
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "DST"));

  // update the parameters on the node tree
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  string g4geonodename = "G4GEO_" + superdetector;
  parNode->addNode(new PHDataNode<PHG4Parameters>(new PHG4Parameters(superdetector), g4geonodename));

  PHG4Parameters *construction_params = findNode::getClass<PHG4Parameters>(
      parNode, g4geonodename);
  assert(construction_params);

  // start with default then fill it up.
  SetDefaultParameters(construction_params);

  // save updated persistant copy on node tree
  const string paramnodename = "G4GEOPARAM_" + superdetector;

  if (useDB)
    {
      // use DB

      int iret = construction_params->ReadFromDB();
      if (iret)
        {
          cout
              << "PHG4SpacalPrototypeSubsystem::InitRun - problem reading from DB"
              << endl;

          return Fun4AllReturnCodes::ABORTRUN;
        }
    }
  else
    {
      PdbParameterMap *nodeparams = findNode::getClass<PdbParameterMap>(topNode,
          paramnodename);
      if (nodeparams)
        {
          construction_params->FillFrom(nodeparams);
        }
    }

  // additional user set parameters
  construction_params->FillFrom(&Params);

//   this step is moved to after detector construction
//   save updated persistant copy on node tree
  construction_params->SaveToNodeTree(parNode, paramnodename);

  if (verbosity > 0)
    cout
        << "PHG4SpacalPrototypeSubsystem::InitRun - use PHG4SpacalPrototypeDetector"
        << endl;
  detector_ = new PHG4SpacalPrototypeDetector(topNode, Name());

  detector_->SetActive(active);
  detector_->SetAbsorberActive(absorberactive);
  detector_->SuperDetector(superdetector);
  detector_->OverlapCheck(overlapcheck);

  if (active)
    {
      ostringstream nodename;
      if (superdetector != "NONE")
        {
          nodename << "G4HIT_" << superdetector;
        }
      else
        {
          nodename << "G4HIT_" << detector_type;
        }
      PHG4HitContainer* cylinder_hits = findNode::getClass<PHG4HitContainer>(
          topNode, nodename.str().c_str());
      if (!cylinder_hits)
        {
          dstNode->addNode(
              new PHIODataNode<PHObject>(
                  cylinder_hits = new PHG4HitContainer(nodename.str()),
                  nodename.str().c_str(), "PHObject"));
        }
      cylinder_hits->AddLayer(0);
      PHG4EventActionClearZeroEdep *evtac = new PHG4EventActionClearZeroEdep(
          topNode, nodename.str());
      if (absorberactive)
        {
          nodename.str("");
          if (superdetector != "NONE")
            {
              nodename << "G4HIT_ABSORBER_" << superdetector;
            }
          else
            {
              nodename << "G4HIT_ABSORBER_" << detector_type ;
            }
          PHG4HitContainer* cylinder_hits =
              findNode::getClass<PHG4HitContainer>(topNode,
                  nodename.str().c_str());
          if (!cylinder_hits)
            {
              dstNode->addNode(
                  new PHIODataNode<PHObject>(cylinder_hits =
                      new PHG4HitContainer(nodename.str()),
                      nodename.str().c_str(), "PHObject"));
            }
          cylinder_hits->AddLayer(0);
          evtac->AddNode(nodename.str());
        }
      eventAction_ = evtac;
      steppingAction_ = new PHG4SpacalPrototypeSteppingAction(detector_);
    }

  return 0;

}

//_______________________________________________________________________
int
PHG4SpacalPrototypeSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector*
PHG4SpacalPrototypeSubsystem::GetDetector(void) const
{
  return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction*
PHG4SpacalPrototypeSubsystem::GetSteppingAction(void) const
{
  return steppingAction_;
}

void
PHG4SpacalPrototypeSubsystem::Print(const std::string &what) const
{
  detector_->Print(what);
  return;
}

void
PHG4SpacalPrototypeSubsystem::SetDefaultParameters(PHG4Parameters * param)
{
  assert(param);

//  param->set_double_param("radius", 95);
//  param->set_double_param("zmin", -40);
//  param->set_double_param("zmax", 40);
//  param->set_double_param("thickness", 16.6);
//
//
//  param->set_string_param("absorber_mat","Spacal_W_Epoxy");
//  param->set_string_param("fiber_clading_mat","PMMA");
//  param->set_string_param("fiber_core_mat","G4_POLYSTYRENE");
//
//  param->set_double_param("xpos", 0);
//  param->set_double_param("ypos", 0);
//  param->set_double_param("zpos", 0);
//
//
//  param->set_double_param("fiber_clading_thickness", 0.003 / 2);
//  param->set_double_param("fiber_core_diameter", 0.047 - (0.003 / 2) * 2);
//  param->set_double_param("fiber_distance",  0.1);
//
//  param->set_int_param("virualize_fiber",  0);
//  param->set_int_param("azimuthal_seg_visible",  0);
//  param->set_int_param("construction_verbose",  verbosity);
//
//  param->set_int_param("azimuthal_n_sec", 8);
//
//  param->set_double_param("assembly_spacing",  0.0001);
//
//  param->set_double_param("sidewall_thickness",  0.075000);
//  param->set_double_param("sidewall_outer_torr",  0.030000);
//  param->set_string_param("sidewall_mat",  "SS310");

//  param->set_int_param("max_phi_bin_in_sec", 1);
//
//  param->set_int_param("init_default_sector_map", 1);


}
