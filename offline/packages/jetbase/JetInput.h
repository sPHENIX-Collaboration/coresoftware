#ifndef JETBASE_JETINPUT_H
#define JETBASE_JETINPUT_H

#include "Jet.h"

#include <globalvertex/GlobalVertex.h>

#include <iostream>
#include <limits>
#include <string>
#include <vector>

class PHCompositeNode;

class JetInput
{
 public:
  virtual ~JetInput() {}

  virtual void identify(std::ostream& os = std::cout)
  {
    os << "JetInput base class" << std::endl;
  }

  virtual Jet::SRC get_src() { return Jet::VOID; }

  virtual std::vector<Jet*> get_input(PHCompositeNode* /*topNode*/)
  {
    return std::vector<Jet*>();
  }
  virtual int Verbosity() const { return m_Verbosity; }
  virtual void Verbosity(int i) { m_Verbosity = i; }

  // vertex actually used by this input in the last get_input() call.
  // has_zvertex() is true for inputs which use a z-vertex in their kinematics
  // (e.g. TowerJetInput, ClusterJetInput). For those, get_vertex_z() returns
  // the z value actually used (0 when the vertex was NaN or missing), and
  // get_vertex_type() returns the GlobalVertex::VTXTYPE name ("MBD", "SVTX", ...)
  // or "UNDEFINED" when no vertex type selection was applied.
  virtual bool has_zvertex() const { return false; }
  virtual std::string get_vertex_type() const { return ""; }
  virtual float get_vertex_z() const { return std::numeric_limits<float>::quiet_NaN(); }

  // name of a GlobalVertex::VTXTYPE, used by get_vertex_type().
  // A type without a case here is reported as "VTXTYPE_<number>", so this
  // does NOT need to be edited when GlobalVertex gains a new type --
  // add a case only if you want a nicer name for it.
  static std::string get_vtxtype_name(GlobalVertex::VTXTYPE type)
  {
    switch (type)
    {
    case GlobalVertex::UNDEFINED:
      return "UNDEFINED";
    case GlobalVertex::TRUTH:
      return "TRUTH";
    case GlobalVertex::SMEARED:
      return "SMEARED";
    case GlobalVertex::CALO:
      return "CALO";
    case GlobalVertex::MBD:
      return "MBD";
    case GlobalVertex::MBD_CALO:
      return "MBD_CALO";
    case GlobalVertex::SVTX:
      return "SVTX";
    case GlobalVertex::SVTX_MBD:
      return "SVTX_MBD";
    default:
      return "VTXTYPE_" + std::to_string(static_cast<int>(type));
    }
  }

 protected:
  JetInput()
    : m_Verbosity(0)
  {
  }

 private:
  int m_Verbosity;
};

#endif
