#include "JetContainerv1.h"
#include <phool/phool.h>  // for PHWHERE
#include <string>
#include "Jetv2.h"

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
  if (m_clones->GetEntriesFast()>0) m_current_jet = (Jet*) m_clones->UncheckedAt(0);

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
    m_clones->Clear("C");
    m_njets = 0;
    m_is_sorted = false;
    m_current_jet = nullptr;
    m_RhoMedian = NAN;
}

Jet* JetContainerv1::add_jet() {
    m_current_jet = (Jet*) m_clones->ConstructedAt(m_njets++,"C");
    m_current_jet->resize_properties(m_psize);
    return m_current_jet;
}

void JetContainerv1::add_component(Jet::SRC src, unsigned int id) {
  m_current_jet->insert_comp(src, id);
}

Jet* JetContainerv1::get_jet(unsigned int ijet) {
    if (ijet < m_njets) {
      m_current_jet = (Jet*) m_clones->At(ijet);
      return m_current_jet;
    } else {
        return nullptr;
    }
}

Jet* JetContainerv1::get_UncheckedAt(unsigned int index) {
    m_current_jet = (Jet*) m_clones->UncheckedAt(index);
    return m_current_jet;
} 

// ----------------------------------------------------------------------------------------
//  Interface for adding, setting, and getting jet properties
// ----------------------------------------------------------------------------------------

void JetContainerv1::print_jets(std::ostream& os) {
  os << " No. of jets: " << m_njets;
  if (!isnan(m_RhoMedian)) os << " rho median " << m_RhoMedian;
  os << std::endl;

  int ijet = 0;
  for (auto jet : *this) {
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
    os << " Jet properties in Jet vectors: " << std::endl;
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
  auto [iter, is_new] = m_pindex.try_emplace(prop, m_psize);
  if (is_new) {
    m_pvec.push_back(prop);
    ++m_psize;
    resize_jet_pvecs();
  }
  return m_psize;
}

size_t JetContainerv1::add_property(std::set<Jet::PROPERTY> prop) {
  auto old_size = m_psize;
  for (auto p : prop) { add_property(p); }
  if (old_size != m_psize) { resize_jet_pvecs(); }
  return m_psize;
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

void JetContainerv1::set_sorted_by(Jet::SORT sort, Jet::SORT_ORDER _order, 
    Jet::PROPERTY sorting_prop) { 
  m_is_sorted = true; 
  m_sort = sort;
  m_sort_order = _order;
  m_sort_prop = sorting_prop;
}

void JetContainerv1::print_sorted_by(std::ostream& os) {
  std::string sort_order = (
      m_sort_order == Jet::SORT_ORDER::DESCENDING ? 
      "large to small" : "small to large");
  switch (m_sort) {
    case Jet::SORT::PT:
    case Jet::SORT::E:
    case Jet::SORT::P:
    case Jet::SORT::MASS:
    case Jet::SORT::MASS2:
    case Jet::SORT::ETA:
      os << " Jets in container sorted by " << 
        str_Jet_SORT(m_sort) << " " << sort_order << std::endl;
      break;
    case Jet::SORT::PROPERTY: // for AREA and others
      os << " Jets in container sorted by PROPERTY " << 
        str_Jet_PROPERTY(m_sort_prop) << " " << sort_order << std::endl;
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

Jet::IterJetTCA JetContainerv1::begin() {
  if (m_clones->GetEntriesFast() > 1) m_current_jet = get_UncheckedAt(0);
  return Jet::IterJetTCA(m_clones, m_current_jet); 
}

Jet::IterJetTCA JetContainerv1::end() {
  return begin();
}

void JetContainerv1::sort_jets(
    Jet::PROPERTY _prop,
    Jet::SORT_ORDER _order) 
{
  sort_jets(Jet::SORT::PROPERTY, _order, _prop);
}

void JetContainerv1::sort_jets(
    Jet::SORT _sort
    , Jet::SORT_ORDER _order
    , Jet::PROPERTY _prop) {

  SortFnJetv2 sort_fn;
  sort_fn.order = _order;
  if (_sort == Jet::SORT::PROPERTY) {
    if (!has_property(_prop)) {
      std::cout << PHWHERE << std::endl;
      std::cout << "Jets don't have property " << str_Jet_PROPERTY(_prop)
        << " not sorting jets." << std::endl;
      return;
    }
    sort_fn.sort = Jet::SORT::PROPERTY;
    sort_fn.index = find_prop_index(_prop);
  } else if (
      _sort == Jet::SORT::PT
   || _sort == Jet::SORT::E
   || _sort == Jet::SORT::P
   || _sort == Jet::SORT::MASS
   || _sort == Jet::SORT::MASS2
   || _sort == Jet::SORT::ETA
  ) {
    sort_fn.sort = _sort;
  } else {
      std::cout << PHWHERE << std::endl;
      std::cout << "Jets can't sort on " << str_Jet_SORT(_sort)
        << ", just PT, E, P, MASS, MASS2, and ETA (or a calculated jet property)." << std::endl
        << " -> Not sorting jets" << std::endl;
      return;
  }
  // do the actual sorting with pointers to the sorting functor and TClonesArray->Sort()
  for (auto jet : (*this)) {
    jet->set_sort_ptr(&sort_fn); // 
  }

  m_clones->UnSort();
  m_clones->Sort();
  m_is_sorted = true;
  m_sort = _sort;
  m_sort_order = _order;
  m_sort_prop = _prop;
  return;
}

/* void JetContainerv1::sort_jets(Jet::PROPERTY prop, Jet::SORT_ORDER _order) { */
/*   if (has_property(prop)) { */
/*     m_sortopt = { Jet::SORT::PROPERTY, prop, _order, m_pindex[prop] }; */
/*     /1* m_sortopt.criteria = Jet::SORT::PROPERTY; *1/ */
/*     /1* m_sortopt.order = _order; *1/ */
/*     /1* m_sortopt.property = prop; *1/ */
/*     /1* m_sortopt.prop_index = m_pindex[prop]; *1/ */

/*     for (auto jet : *this) { */
/*         jet->set_sort_selptr(&m_sortopt); */
/*     } */
/*     m_clones->UnSort(); */
/*     m_clones->Sort(); */
/*     m_is_sorted = true; */
/*   } else { */
/*     m_is_sorted = false; */
/*       std::cout << PHWHERE << std::endl */
/*         << "Error: property selected for sorting not set for these jets. " << std::endl; */
/*   } */
/* } */

void JetContainerv1::resize_jet_pvecs() {
  for (auto jet : *this) {
      jet->resize_properties(m_psize);
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
