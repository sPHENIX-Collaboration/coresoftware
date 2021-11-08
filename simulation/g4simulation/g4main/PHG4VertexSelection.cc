// This is the new trackbase container version

/*!
 * \file PHG4VertexSelection.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHG4VertexSelection.h"

#include "PHG4TruthInfoContainer.h"
#include "PHG4VtxPoint.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <iostream>  // for operator<<, basic_ostream, endl

//____________________________________________________________________________
PHG4VertexSelection::PHG4VertexSelection(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
{
  InitializeParameters();
}

//____________________________________________________________________________
int PHG4VertexSelection::InitRun(PHCompositeNode */*topNode*/)
{
  UpdateParametersWithMacro();

  // load parameters
  m_vertex_zcut = get_double_param("vertex_zcut");

  // printout
  std::cout
      << "PHG4VertexSelection::InitRun\n"
      << " m_vertex_zcut: " << m_vertex_zcut << " cm\n"
      << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________
int PHG4VertexSelection::process_event(PHCompositeNode *topNode)
{
  // g4 truth info
  auto g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  // main vertex
  const auto main_vertex_id = g4truthinfo->GetPrimaryVertexIndex();
  const auto vertex = g4truthinfo->GetPrimaryVtx(main_vertex_id);
  if (!vertex) return false;

  // check vertex position along the beam
  return (m_vertex_zcut > 0 && std::abs(vertex->get_z()) > m_vertex_zcut) ? Fun4AllReturnCodes::DISCARDEVENT : Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
void PHG4VertexSelection::SetDefaultParameters()
{
  set_default_double_param("vertex_zcut", 10);
}
