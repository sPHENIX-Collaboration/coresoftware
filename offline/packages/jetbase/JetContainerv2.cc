#include "JetContainerv2.h"
#include "Jetv3.h"

#include <phool/phool.h>  // for PHWHERE

#include <format>
#include <limits>
#include <string>

JetContainerv2::JetContainerv2()
  : m_clones(new TClonesArray("Jetv3", 50))
{
}

void JetContainerv2::identify(std::ostream& os) const
{
  os << "JetContainerv2 (Jetv3 jets with corrected pT): size = " << m_clones->GetEntriesFast() << std::endl;
  if (!m_vertex_type.empty())
  {
    os << "  reconstructed with vertex type " << m_vertex_type << " at z = " << m_vertex_z << std::endl;
  }
  os << "  Contains jets with the following properties:" << std::endl;
  print_property_types(os);
  return;
}

JetContainerv2::~JetContainerv2()
{
  JetContainerv2::Reset();
  delete m_clones;
}

void JetContainerv2::Reset()
{
  m_clones->Clear("C");
  m_njets = 0;
  m_RhoMedian = std::numeric_limits<float>::quiet_NaN();
  m_has_zvertex = false;
  m_vertex_type.clear();
  m_vertex_z = std::numeric_limits<float>::quiet_NaN();
}

Jet* JetContainerv2::add_jet()
{
  auto* jet = (Jet*) m_clones->ConstructedAt(m_njets++, "C");
  jet->resize_properties(m_psize);
  return jet;
}

Jet* JetContainerv2::get_jet(unsigned int ijet)
{
  if (ijet < m_njets)
  {
    return (Jet*) m_clones->At(ijet);
  }
  return nullptr;
}

Jet* JetContainerv2::get_UncheckedAt(unsigned int index)
{
  return (Jet*) m_clones->UncheckedAt(index);
}

// ----------------------------------------------------------------------------------------
//  Interface for adding, setting, and getting jet properties
// ----------------------------------------------------------------------------------------

void JetContainerv2::print_jets(std::ostream& os)
{
  os << " No. of jets: " << m_njets;
  if (!std::isnan(m_RhoMedian))
  {
    os << " rho median " << m_RhoMedian;
  }
  os << std::endl;

  int ijet = 0;
  for (auto* jet : *this)
  {
    os << std::format("  jet({:2}) : pT({:6.2f})  eta({:6.2f})  phi({:6.2f})  pT_calib({:6.2f})  pT_insitu({:6.2f})  pT_smear({:6.2f})",
                      ijet, jet->get_pt(), jet->get_eta(), jet->get_phi(), jet->get_pt_calib(), jet->get_pt_insitu(), jet->get_pt_smear());
    ++ijet;
    for (auto prop : m_pindex)
    {
      os << std::format("  {:>8}({:6.2f})", str_Jet_PROPERTY(prop.first), jet->get_property(prop.second));
    }
    os << std::endl;
  }
  os << std::endl;
}

void JetContainerv2::print_property_types(std::ostream& os) const
{
  os << " Jet properties in Jet vectors: " << std::endl;
  int i = 0;
  for (auto p : m_pindex)
  {
    os << " (" << i++ << ") -> " << str_Jet_PROPERTY(p.first) << std::endl;
  }
  return;
}

// Add properties to the jets.
size_t JetContainerv2::add_property(Jet::PROPERTY prop)
{
  auto [iter, is_new] = m_pindex.try_emplace(prop, static_cast<Jet::PROPERTY>(m_psize));
  if (is_new)
  {
    ++m_psize;
    resize_jet_pvecs();
  }
  return m_psize;
}

size_t JetContainerv2::add_property(std::set<Jet::PROPERTY> props)
{
  bool added = false;
  for (auto prop : props)
  {
    auto [iter, is_new] = m_pindex.try_emplace(prop, static_cast<Jet::PROPERTY>(m_psize));
    if (is_new)
    {
      ++m_psize;
      added = true;
    }
  }
  if (added)
  {
    resize_jet_pvecs();
  }
  return m_psize;
}

// get the index for a given property
Jet::PROPERTY JetContainerv2::property_index(Jet::PROPERTY prop)
{
  if (!has_property(prop))
  {
    add_property(prop);
  }
  return m_pindex[prop];
}

Jet::IterJetTCA JetContainerv2::begin()
{
  return Jet::IterJetTCA(m_clones);
}

Jet::IterJetTCA JetContainerv2::end()  // dummy implementation -- don't anticipate that it will ever be checked
{
  auto rval = Jet::IterJetTCA(m_clones);
  rval.index = rval.size;
  return rval;
}

void JetContainerv2::resize_jet_pvecs()
{
  for (auto* jet : *this)
  {
    jet->resize_properties(m_psize);
  }
  return;
}

std::string JetContainerv2::str_Jet_PROPERTY(Jet::PROPERTY prop)
{
  switch (prop)
  {
  case Jet::PROPERTY::prop_JetCharge:
    return "JetCharge";
  case Jet::PROPERTY::prop_BFrac:
    return "BFrac";
  case Jet::PROPERTY::prop_SeedD:
    return "SeedD";
  case Jet::PROPERTY::prop_SeedItr:
    return "SeedItr";
  case Jet::PROPERTY::prop_zg:
    return "zg";
  case Jet::PROPERTY::prop_Rg:
    return "Rg";
  case Jet::PROPERTY::prop_mu:
    return "mu";
  case Jet::PROPERTY::prop_gamma:
    return "gamma";
  case Jet::PROPERTY::prop_JetHadronFlavor:
    return "JetHadronFlavor";
  case Jet::PROPERTY::prop_JetHadronZT:
    return "JetHadronZT";
  case Jet::PROPERTY::prop_area:
    return "area";
  case Jet::PROPERTY::prop_t:
    return "t";
  default:
    return "no_property";
  }
  return "";
}
