#include <HepMC3/AttributeFeature.h>
#include <HepMC3/Selector.h>
#include <functional>
#include <iterator>
#include <memory>
#include <sstream> // __str__
#include <string>

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

void bind_pyHepMC3search_1(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// HepMC3::abs(const class HepMC3::Selector &) file:HepMC3/Selector.h line:161
	M("HepMC3").def("abs", (class std::shared_ptr<const class HepMC3::Selector> (*)(const class HepMC3::Selector &)) &HepMC3::abs, "ConstSelectorPtr abs\n\nC++: HepMC3::abs(const class HepMC3::Selector &) --> class std::shared_ptr<const class HepMC3::Selector>", pybind11::arg("input"));

	{ // HepMC3::StandardSelector file:HepMC3/Selector.h line:174
		pybind11::class_<HepMC3::StandardSelector, std::shared_ptr<HepMC3::StandardSelector>, HepMC3::Selector> cl(M("HepMC3"), "StandardSelector", "StandardSelector ");
		cl.def("assign", (class HepMC3::StandardSelector & (HepMC3::StandardSelector::*)(const class HepMC3::StandardSelector &)) &HepMC3::StandardSelector::operator=, "C++: HepMC3::StandardSelector::operator=(const class HepMC3::StandardSelector &) --> class HepMC3::StandardSelector &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
