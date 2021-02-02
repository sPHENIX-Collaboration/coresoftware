// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4MAIN_PHG4VERTEXSELECTION_H
#define G4MAIN_PHG4VERTEXSELECTION_H

/*!
 * \file PHG4VertexSelection.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/SubsysReco.h>
#include <phparameter/PHParameterInterface.h>

#include <string>  // for string

class PHCompositeNode;

class PHG4VertexSelection : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4VertexSelection(const std::string &name = "PHG4VertexSelection");

  //! run initialization
  int InitRun(PHCompositeNode *) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  //! parameters
  void SetDefaultParameters() override;

 private:
  // z vertex cut (cm)
  double m_vertex_zcut = 10;
};

#endif
