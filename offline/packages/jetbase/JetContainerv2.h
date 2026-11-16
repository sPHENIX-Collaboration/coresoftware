#ifndef JETBASE_JETCONTAINERV2_H
#define JETBASE_JETCONTAINERV2_H

#include "JetContainerv1.h"

#include <iostream>

// JetContainerv1 holding Jetv3 jets: each jet carries its corrected pT
// (Jet::get_pt_calib(), get_pt_insitu[_up|_down](), get_pt_smear[_up|_down]()), so JetCalib and the
// later correction modules fill the jets of this node instead of creating a new node.
// Everything else (jet access, properties, sources, rho, vertex) is JetContainerv1; code that asks
// for a JetContainerv1 also gets this container.
class JetContainerv2 : public JetContainerv1
{
 public:
  JetContainerv2();
  ~JetContainerv2() override = default;
  void identify(std::ostream& os = std::cout) const override;

 private:
  ClassDefOverride(JetContainerv2, 1);
};

#endif  // JETBASE_JETCONTAINERV2_H
