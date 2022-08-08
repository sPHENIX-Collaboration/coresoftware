#include <HepMC3/AttributeFeature.h>
#include <HepMC3/Relatives.h>
#include <HepMC3/Selector.h>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <sstream> // __str__
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <functional>
#include <string>
#include <HepMC3/Version.h>
#include <HepMC3/Relatives.h>
#include <HepMC3/Selector.h>
#include <src/search_binders.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

void bind_pyHepMC3search_0(std::function< pybind11::module &(std::string const &namespace_) > &M)
{

	 binder::search_binder(M("HepMC3"));
	{ // HepMC3::AttributeFeature file:HepMC3/AttributeFeature.h line:22
		pybind11::class_<HepMC3::AttributeFeature, std::shared_ptr<HepMC3::AttributeFeature>> cl(M("HepMC3"), "AttributeFeature", "");
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("name") );

		cl.def( pybind11::init( [](HepMC3::AttributeFeature const &o){ return new HepMC3::AttributeFeature(o); } ) );
		cl.def("exists", (class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)> (HepMC3::AttributeFeature::*)() const) &HepMC3::AttributeFeature::exists, "existence\n\nC++: HepMC3::AttributeFeature::exists() const --> class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)>");
		cl.def("__eq__", (class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)> (HepMC3::AttributeFeature::*)(std::string) const) &HepMC3::AttributeFeature::operator==, "equality operator\n\nC++: HepMC3::AttributeFeature::operator==(std::string) const --> class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)>", pybind11::arg("rhs"));
		cl.def("assign", (class HepMC3::AttributeFeature & (HepMC3::AttributeFeature::*)(const class HepMC3::AttributeFeature &)) &HepMC3::AttributeFeature::operator=, "C++: HepMC3::AttributeFeature::operator=(const class HepMC3::AttributeFeature &) --> class HepMC3::AttributeFeature &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::_parents file:HepMC3/Relatives.h line:76
		pybind11::class_<HepMC3::_parents, std::shared_ptr<HepMC3::_parents>> cl(M("HepMC3"), "_parents", "Provides operator to find the parent particles of a Vertex or Particle\n\n Note you would usually not instantiate this directly, but wrap it in a RelativesInterface");
		cl.def( pybind11::init( [](){ return new HepMC3::_parents(); } ) );
	}
	{ // HepMC3::_children file:HepMC3/Relatives.h line:99
		pybind11::class_<HepMC3::_children, std::shared_ptr<HepMC3::_children>> cl(M("HepMC3"), "_children", "Provides operator to find the child particles of a Vertex or Particle\n\n Note you would usually not instantiate this directly, but wrap it in a RelativesInterface");
		cl.def( pybind11::init( [](){ return new HepMC3::_children(); } ) );
	}
	{ // HepMC3::Relatives file:HepMC3/Relatives.h line:182
		pybind11::class_<HepMC3::Relatives, std::shared_ptr<HepMC3::Relatives>> cl(M("HepMC3"), "Relatives", "Define a common interface that all Relatives objects will satisfy\n         Relatives provides an operator to get the relatives of a range of different\n         GenObject types.  The following are examples\n\n         Relatives::ANCESTORS(GenParticlePtr);// returns ancestors of the particle\n         Descendants descendants;\n         descendants(GenVertexPtr);// descendants of the vertex\n         vector<Relatives*> relations = {&Relatives::CHILDREN, &Relatives::DESCENDANTS, &Relatives::PARENTS, new Ancestors()}; // make a vector of Relatives\n\n         You can also define your own relation and wrap it in the Relatives interface using\n         Relatives * relo = new RelativesInterface<MyRelationClass>();");
		cl.def("assign", (class HepMC3::Relatives & (HepMC3::Relatives::*)(const class HepMC3::Relatives &)) &HepMC3::Relatives::operator=, "C++: HepMC3::Relatives::operator=(const class HepMC3::Relatives &) --> class HepMC3::Relatives &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::Selector file:HepMC3/Selector.h line:57
		pybind11::class_<HepMC3::Selector, std::shared_ptr<HepMC3::Selector>> cl(M("HepMC3"), "Selector", "Selector is an interface to \"standard\" Features that are valid\n  for both integral and floating point comparisons\n\n  You would use this in preference to the more general\n  Feature<> templated type.  A Selector is constructed from a\n  function to extract features from particles, e.g.\n\n  ConstSelectorPtr status = std::make_shared<SelectorWrapper<int> >([](ConstParticlePtr p)->int{return p->status();});\n  ConstSelectorPtr pt = std::make_shared<SelectorWrapper<double> >([](ConstParticlePtr p)->double{return p->momentum().pt();});\n\n  You can then use the Selector to construct Filter functions that\n  evaluate on particles, e.g.\n  Filter is_stable = (*status) == 1;\n  bool stable = is_stable(p);\n  bool beam = (*status == 4)(p);\n\n  StandardSelector contains a few standard Selectors already defined, e.g.\n\n  ConstGenParticlePtr p;\n  (StandardSelector::STATUS == 1)(p);\n  (StandardSelector::PT > 15.)(p);\n  (abs(StandardSelector::RAPIDITY) < 2.5)(p);\n\n  you can also combined them e.g.\n\n  Filter myCuts = (StandardSelector::PT > 15.) && (*abs(StandardSelector::RAPIDITY) < 2.5) || (StandardSelector::PT > 100.);\n  bool passCuts = myCuts(p);");
		cl.def("__eq__", (class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)> (HepMC3::Selector::*)(int) const) &HepMC3::Selector::operator==, "C++: HepMC3::Selector::operator==(int) const --> class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)>", pybind11::arg("value"));
		cl.def("__eq__", (class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)> (HepMC3::Selector::*)(double) const) &HepMC3::Selector::operator==, "C++: HepMC3::Selector::operator==(double) const --> class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)>", pybind11::arg("value"));
		cl.def("__ne__", (class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)> (HepMC3::Selector::*)(int) const) &HepMC3::Selector::operator!=, "C++: HepMC3::Selector::operator!=(int) const --> class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)>", pybind11::arg("value"));
		cl.def("__ne__", (class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)> (HepMC3::Selector::*)(double) const) &HepMC3::Selector::operator!=, "C++: HepMC3::Selector::operator!=(double) const --> class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)>", pybind11::arg("value"));
		cl.def("abs", (class std::shared_ptr<const class HepMC3::Selector> (HepMC3::Selector::*)() const) &HepMC3::Selector::abs, "C++: HepMC3::Selector::abs() const --> class std::shared_ptr<const class HepMC3::Selector>");
		cl.def_static("ATTRIBUTE", (class HepMC3::AttributeFeature (*)(const std::string &)) &HepMC3::Selector::ATTRIBUTE, "C++: HepMC3::Selector::ATTRIBUTE(const std::string &) --> class HepMC3::AttributeFeature", pybind11::arg("name"));
		cl.def("assign", (class HepMC3::Selector & (HepMC3::Selector::*)(const class HepMC3::Selector &)) &HepMC3::Selector::operator=, "C++: HepMC3::Selector::operator=(const class HepMC3::Selector &) --> class HepMC3::Selector &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		binder::custom_Selector_binder(cl);
	}
}
