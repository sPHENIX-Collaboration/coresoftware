#include "JetContainerv1.h"
#include "Jetv2.h"
#include <phool/phool.h>  // for PHWHERE
#include <string>

JetContainerv1::JetContainerv1() {
    m_clones = new TClonesArray("Jetv2", 50);
}

void JetContainerv1::identify(std::ostream& os) const {
    os << "JetContainerv1: size = " << m_clones->GetEntriesFast() << std::endl
       << "  Contains jets with the following properties:" << std::endl;
    print_property_types(os);
    return;
}

JetContainerv1::JetContainerv1(const JetContainer &rhs)
  : m_njets         { rhs.size()                      }
  , m_pvec          { rhs.vec_jet_properties()        }
  , m_pindex        { rhs.prop_indices_map()          }
  , m_psize         { rhs.vec_jet_properties().size() }
  , m_sel_index     { rhs.get_index_single()          }
  , m_sel_index_vec { rhs.get_index_vec()             }
  , m_is_sorted     { rhs.is_sorted()                 }
  , m_RhoMedian     { rhs.get_rho_median()            }
{
  m_clones = rhs.clone_data();
  m_current_jet = nullptr;
  if (m_clones->GetEntriesFast()>0) m_current_jet = (Jetv2*) m_clones->UncheckedAt(0);

  for (auto src = rhs.begin_src(); src != rhs.end_src(); ++src) {
    m_src.insert(*src);
  }

  m_jetpar_R = rhs.get_jetpar_R();
}

JetContainerv1::~JetContainerv1()
{
  JetContainerv1::Reset();
}

void JetContainerv1::Reset() {
    m_clones->Clear();
    m_njets = 0;
    m_is_sorted = false;
    m_current_jet = nullptr;
    m_RhoMedian = NAN;
}

Jetv2* JetContainerv1::add_jet() {
    m_current_jet = (Jetv2*) m_clones->ConstructedAt(m_njets++);
    m_current_jet->resize_properties(m_psize);
    return m_current_jet;
}

Jetv2* JetContainerv1::get_jet(unsigned int ijet) {
    if (ijet < m_njets) {
        return (Jetv2*) m_clones->At(ijet);
    } else {
        return nullptr;
    }
}

Jetv2* JetContainerv1::get_UncheckedAt(unsigned int index) {
    return (Jetv2*) m_clones->UncheckedAt(index);
} 

// ----------------------------------------------------------------------------------------
//  Interface for adding, setting, and getting jet properties
// ----------------------------------------------------------------------------------------

void JetContainerv1::print_jets(std::ostream& os) {
  os << " No. of jets: " << m_njets;
  if (!isnan(m_RhoMedian)) os << " rho median " << m_RhoMedian;
  os << std::endl;

  int ijet = 0;
  for (auto jet : iter_jets()) {
    os << Form("  jet(%2i) : pT(%6.2f)  eta(%6.2f)  phi(%6.2f)", 
        ijet, jet->get_pt(), jet->get_eta(), jet->get_phi());
    ++ijet;
    unsigned int i {0};
    for (auto prop : m_pvec) {
      os << Form("  %8s(%6.2f)", str_Jet_PROPERTY(prop).c_str(), jet->get_prop_by_index(i));
      i++;
    }
    os << std::endl;
  }
  os << std::endl;
}

void JetContainerv1::print_property_types(std::ostream& os) const {
    os << " Jet properties in Jetv2 vectors: " << std::endl;
    int i {0};
    for (auto p : m_pvec) {
        os << " ("<<i++<<") -> " << str_Jet_PROPERTY(p) << std::endl;
    }
    return;
}

bool JetContainerv1::has_property(Jet::PROPERTY prop) const {
    if (m_pindex.find(prop) == m_pindex.end()) {
        return false;
    }
    return true;
}

void JetContainerv1::print_missing_prop(Jet::PROPERTY prop) const {
    std::cout << "JetContainerv1::find_prop_index - ERROR - property " 
      << prop << "(" << str_Jet_PROPERTY(prop) << ") not found" << std::endl;
}

// Add properties to the jets. 
size_t JetContainerv1::add_property(Jet::PROPERTY prop) {
    auto emplace = m_pindex.try_emplace(prop, m_psize);
    if (!emplace.second) emplace.first->second += m_psize;
    resize_jet_pvecs();
    return m_pvec.size();
}

size_t JetContainerv1::add_property(std::set<Jet::PROPERTY> prop) {
    for (auto p : prop) { add_property(p); }
    resize_jet_pvecs();
    return m_pvec.size();
}

// Get and set values of properties by index (always on current_jet)
unsigned int JetContainerv1::find_prop_index(Jet::PROPERTY prop) {
    if (has_property(prop)) {
        return m_pindex[prop];
    } else {
        print_missing_prop(prop);
        std::cout << " Returning UINT_MAX" << std::endl;
        return UINT_MAX;
    }
}

void JetContainerv1::set_sorted_by(Jet::SORT sort, bool is_inverse, Jet::PROPERTY sorting_prop) { 
  m_is_sorted = true; 
  m_sorted_by = sort; 
  m_sorted_inverse= is_inverse;
  m_sorting_prop = sorting_prop;
}

void JetContainerv1::print_sorted_by(std::ostream& os) {
  std::string inversely = (m_sorted_inverse ? "large to small" : "small to large");
  switch (m_sorted_by) {
    case Jet::SORT::PT:
    case Jet::SORT::E:
    case Jet::SORT::P:
    case Jet::SORT::MASS:
    case Jet::SORT::MASS2:
    case Jet::SORT::ETA:
      os << " Jets in container sorted by " << str_Jet_SORT(m_sorted_by) << " " << inversely << std::endl;
      break;
    case Jet::SORT::PROPERTY: // for AREA and others
      os << " Jets in container sorted by PROPERTY " << str_Jet_PROPERTY(m_sorting_prop) << " " << inversely << std::endl;
      break;
    default:
      os << " Jet sorting information not set for jet container. " << std::endl;
      break;
  }
}

float JetContainerv1::get_prop_by_index(unsigned int index) const {
    return m_current_jet->get_prop_by_index(index);
}

void JetContainerv1::set_prop_by_index(unsigned int index, float value) {
    m_current_jet->set_prop_by_index(index, value);
}

// Get and set values of properties by selections (always on current_jet)
void JetContainerv1::select_property(Jet::PROPERTY prop) {
    m_sel_index_vec.resize(0); // can't use the vector
    if (!has_property(prop)) {
        print_missing_prop(prop);
        std::cout << " Cannot select this property " << std::endl;
        return;
    }
    m_sel_index = m_pindex[prop];
}

void JetContainerv1::select_property(std::vector<Jet::PROPERTY> vprop) {
    int i {0};
    m_sel_index_vec.resize(vprop.size());
    for (auto& v : vprop) {
        if (!has_property(v)) {
            print_missing_prop(v);
            std::cout << " Fatal error: cannot select this property " << std::endl;
            assert(false);
        }
        m_sel_index_vec[i++] = m_pindex[v];
    }
    m_sel_index = m_sel_index_vec[0];
}

float JetContainerv1::get_selected_property(unsigned int index) {
    if (index==0) return m_current_jet->get_prop_by_index(m_sel_index);
    else return m_current_jet->get_prop_by_index(m_sel_index_vec[index]);
}

void JetContainerv1::set_selected_property(float value, unsigned int index) {
    if (index==0) m_current_jet->set_prop_by_index(m_sel_index, value);
    else m_current_jet->set_prop_by_index(m_sel_index_vec[index], value);
}

void JetContainerv1::sort_jets(Jet::SORT isort, bool descending) {
  // check that it's an ok sorting criteria:
  m_is_sorted = true;
  m_sorted_by = isort;
  m_sorted_inverse = descending;
  switch (isort) {
    case Jet::SORT::PT:
    case Jet::SORT::E:
    case Jet::SORT::P:
    case Jet::SORT::MASS:
    case Jet::SORT::MASS2:
    case Jet::SORT::ETA:
      for (auto jet : iter_jets()) {
          jet->set_sort_criteria(isort, descending);
      }
      m_clones->UnSort();
      m_clones->Sort();
      m_is_sorted = true;
      break;

    case Jet::SORT::AREA:
      sort_jets(Jet::PROPERTY::prop_area, descending);
      break;

    default:
      std::cout << PHWHERE << std::endl
        << "Error in options passed to sort_jets. See options in code." << std::endl;
  }
  return;
}

void JetContainerv1::sort_jets(Jet::PROPERTY prop, bool descending) {
  if (has_property(prop)) {
    m_is_sorted = true;
    m_sorted_by = Jet::SORT::PROPERTY;
    m_sorted_inverse = descending;
    m_sorting_prop = prop;
      auto index = m_pindex[prop];
      for (auto jet : iter_jets()) {
          jet->set_sort_criteria(Jet::SORT::PROPERTY, descending, index);
      }
      m_clones->UnSort();
      m_clones->Sort();
      m_is_sorted = true;
  } else {
    m_is_sorted = false;
      std::cout << PHWHERE << std::endl
        << "Error: property selected for sorting not set for these jets. " << std::endl;
  }
}

// ----------------------------------------------------------------------------------------
//  Sorting the jets
// ----------------------------------------------------------------------------------------
/* void JetContainerv1::sort_jets(Jet::SORT iSort, bool desc) */ 
/* { */
/*     int neg = (desc) ? -1 : 1; */
/*     if (iSort == Jet::SORT::PT) { */
/*       sort_jets( [neg](const Jetv2* lhs, const Jetv2* rhs)->bool{ return lhs->get_pt() == rhs->get_pt();}, */
/*                  [neg](const Jetv2* lhs, const Jetv2* rhs)->int { return neg*(lhs->get_pt()-rhs->get_pt());} ); */
/*     } else if (iSort == Jet::SORT::E) { */
/*       sort_jets( [neg](const Jetv2* lhs, const Jetv2* rhs)->bool{ return lhs->get_e() == rhs->get_e();}, */
/*                  [neg](const Jetv2* lhs, const Jetv2* rhs)->int { return neg*(lhs->get_e()-rhs->get_e());} ); */
/*     } else if (iSort == Jet::SORT::P) { */
/*       sort_jets( [neg](const Jetv2* lhs, const Jetv2* rhs)->bool{ return lhs->get_p() == rhs->get_p();}, */
/*                  [neg](const Jetv2* lhs, const Jetv2* rhs)->int { return neg*(lhs->get_p()-rhs->get_p());} ); */
/*     } else if (iSort == Jet::SORT::MASS) { */
/*       sort_jets( [neg](const Jetv2* lhs, const Jetv2* rhs)->bool{ return lhs->get_mass() == rhs->get_mass();}, */
/*                  [neg](const Jetv2* lhs, const Jetv2* rhs)->int { return neg*(lhs->get_mass()-rhs->get_mass());} ); */
/*     } else if (iSort == Jet::SORT::MASS2) { */
/*       sort_jets( [neg](const Jetv2* lhs, const Jetv2* rhs)->bool{ return lhs->get_mass2() == rhs->get_mass2();}, */
/*                  [neg](const Jetv2* lhs, const Jetv2* rhs)->int { return neg*(lhs->get_mass2()-rhs->get_mass2());} ); */
/*     } else if (iSort == Jet::SORT::AREA) { */
/*       sort_jets(Jet::PROPERTY::prop_area, desc); */
/*       return; */
/*     } */ 

/*     std::cout << "JetContainerv1::sort_jets - unrecognized sort type " << ((int)iSort) << std::endl; */
/*     return; */
/* } */

/* void JetContainerv1::sort_jets(Jet::PROPERTY iProp, bool desc) { */
/*     if (!has_property(iProp)) { */
/*         std::cout << " Because can't find property, aborting jet sorting." << std::endl; */
/*         return; */
/*     } */

/*     auto pin = m_pindex[iProp]; */
/*     int neg = (desc) ? -1 : 1; */

/*     // Note: any jets with value FLT_MAX will be sorted to the end of the list automatically */
/*     sort_jets( */ 
/*       //lambda function for equal */
/*       [pin](const Jetv2* lhs, const Jetv2* rhs)->bool */
/*         { return lhs->get_prop_by_index(pin) == rhs->get_prop_by_index(pin); }, */
/*       //lambda function for compare */
/*       [pin,neg] (const Jetv2* lhs, const Jetv2* rhs)->int */
/*         { if      (lhs->get_prop_by_index(pin) == NAN) return -neg; // FIXME -- think this logic is correct, but should check */
/*           else if (rhs->get_prop_by_index(pin) == NAN) return neg; */
/*           return neg*(lhs->get_prop_by_index(pin)-rhs->get_prop_by_index(pin)); */ 
/*         } */
/*     ); */
/*     return; */
/* } */

/* void JetContainerv1::sort_jets( */ 
/*   std::function<bool(const Jetv2*, const Jetv2*)> fnIsEqual, */
/*   std::function<int (const Jetv2*, const Jetv2*)> fnCompare */
/* ) { */
/* } */

void JetContainerv1::resize_jet_pvecs() {
    for (auto j : iter_jets()) {
        j->resize_properties(m_psize);
    }
    return;
}

std::string JetContainerv1::str_Jet_PROPERTY(Jet::PROPERTY prop) const {
  switch(prop) {
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

std::string JetContainerv1::str_Jet_SORT(Jet::SORT sort) const {
  switch(sort) {
    case Jet::SORT::NO_SORT:
      return "no sorting";
    case Jet::SORT::PT:
      return "PT";
    case Jet::SORT::E:
      return "E";
    case Jet::SORT::P:
      return "P";
    case Jet::SORT::MASS:
      return "MASS";
    case Jet::SORT::AREA:
      return "AREA";
    case Jet::SORT::PROPERTY:
      return "PROPERTY";
    case Jet::SORT::MASS2:
      return "MASS2";
    case Jet::SORT::ETA:
      return "ETA";
  }
  return "";
}
