#include "JetContainerv1.h"
#include "Jetv2.h"

#include <phool/phool.h>  // for PHWHERE
                         
#include <string>

JetContainerv1::JetContainerv1()
{
  m_clones = new TClonesArray("Jetv2", 50);
}

void JetContainerv1::identify(std::ostream& os) const
{
  os << "JetContainerv1: size = " << m_clones->GetEntriesFast() << std::endl
     << "  Contains jets with the following properties:" << std::endl;
  print_property_types(os);
  return;
}

JetContainerv1::JetContainerv1(const JetContainer& rhs)
  : m_njets{rhs.size()}
  , m_pindex{rhs.property_indices()}
  , m_psize{rhs.size_properties()}
  , m_RhoMedian{rhs.get_rho_median()}
{
  m_clones = rhs.clone_data();
  m_current_jet = nullptr;
  if (m_clones->GetEntriesFast() > 0) m_current_jet = (Jet*) m_clones->UncheckedAt(0);

  for (auto src = rhs.begin_src(); src != rhs.end_src(); ++src)
  {
    m_src.insert(*src);
  }

  m_jetpar_R = rhs.get_jetpar_R();
}

JetContainerv1::~JetContainerv1()
{
  JetContainerv1::Reset();
  delete m_clones;
}

void JetContainerv1::Reset()
{
  m_clones->Clear("C");
  m_njets = 0;
  m_current_jet = nullptr;
  m_RhoMedian = NAN;
}

Jet* JetContainerv1::add_jet()
{
  m_current_jet = (Jet*) m_clones->ConstructedAt(m_njets++, "C");
  m_current_jet->resize_properties(m_psize);
  return m_current_jet;
}

Jet* JetContainerv1::get_jet(unsigned int ijet)
{
  if (ijet < m_njets)
  {
    m_current_jet = (Jet*) m_clones->At(ijet);
    return m_current_jet;
  }
  else
  {
    return nullptr;
  }
}

Jet* JetContainerv1::get_UncheckedAt(unsigned int index)
{
  m_current_jet = (Jet*) m_clones->UncheckedAt(index);
  return m_current_jet;
}

// ----------------------------------------------------------------------------------------
//  Interface for adding, setting, and getting jet properties
// ----------------------------------------------------------------------------------------

void JetContainerv1::print_jets(std::ostream& os)
{
  os << " No. of jets: " << m_njets;
  if (!isnan(m_RhoMedian)) os << " rho median " << m_RhoMedian;
  os << std::endl;

  int ijet = 0;
  for (auto jet : *this)
  {
    os << Form("  jet(%2i) : pT(%6.2f)  eta(%6.2f)  phi(%6.2f)",
               ijet, jet->get_pt(), jet->get_eta(), jet->get_phi());
    ++ijet;
    unsigned int i = 0;
    for (auto prop : m_pindex)
    {
      os << Form("  %8s(%6.2f)", str_Jet_PROPERTY(prop.first).c_str(), jet->get_property(prop.second));
      i++;
    }
    os << std::endl;
  }
  os << std::endl;
}

void JetContainerv1::print_property_types(std::ostream& os) const
{
  os << " Jet properties in Jet vectors: " << std::endl;
  int i = 0;
  for (auto p : m_pindex)
  {
    os << " (" << i++ << ") -> " << str_Jet_PROPERTY(p.first) << std::endl;
  }
  return;
}

/* bool JetContainerv1::has_property(Jet::PROPERTY prop) const */
/* { */
/*   if (m_pindex.find(prop) == m_pindex.end()) */
/*   { */
/*     return false; */
/*   } */
/*   return true; */
/* } */

/* void JetContainerv1::print_missing_prop(Jet::PROPERTY prop) const */
/* { */
/*   std::cout << "JetContainerv1::find_prop_index - ERROR - property " */
/*             << prop << "(" << str_Jet_PROPERTY(prop) << ") not found" << std::endl; */
/* } */

// Add properties to the jets.
size_t JetContainerv1::add_property(Jet::PROPERTY prop)
{
  auto [iter, is_new] = m_pindex.try_emplace(prop, static_cast<Jet::PROPERTY>(m_psize));
  if (is_new)
  {
    ++m_psize;
    resize_jet_pvecs();
  }
  return m_psize;
}

size_t JetContainerv1::add_property(std::set<Jet::PROPERTY> props)
{
  bool added = false;
  for (auto prop : props) {
    auto [iter, is_new] = m_pindex.try_emplace(prop, static_cast<Jet::PROPERTY>(m_psize));
    if (is_new) {
      ++m_psize;
      added = true;
    }
  }
  if (added) {
    resize_jet_pvecs();
  }
  return m_psize;
}

// get the index for a given property
Jet::PROPERTY JetContainerv1::property_index(Jet::PROPERTY prop)
{
  if (has_property(prop))
  {
    return m_pindex[prop];
  }
  else
  {
     std::cout << "JetContainerv1::poperty_index - ERROR - property "
            << prop << "(" << str_Jet_PROPERTY(prop) << ") not found." << std::endl;
    std::cout << " Returning 1000" << std::endl;
    return static_cast<Jet::PROPERTY>(1000);
  }
}

Jet::IterJetTCA JetContainerv1::begin()
{
  if (m_clones->GetEntriesFast() > 1) m_current_jet = get_UncheckedAt(0);
  return Jet::IterJetTCA(m_clones, m_current_jet);
}

Jet::IterJetTCA JetContainerv1::end()
{
  return begin();
}

void JetContainerv1::resize_jet_pvecs()
{
  for (auto jet : *this)
  {
    jet->resize_properties(m_psize);
  }
  return;
}

std::string JetContainerv1::str_Jet_PROPERTY(Jet::PROPERTY prop) const
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
  default:
    return "no_property";
  }
  return "";
}
