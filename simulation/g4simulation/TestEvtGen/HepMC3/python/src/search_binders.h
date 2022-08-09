#ifndef SEARCH_BINDERS_H
#define SEARCH_BINDERS_H
#include <HepMC3/Relatives.h>
#include <HepMC3/Selector.h>
#include <pybind11/pybind11.h>
namespace binder {
	void custom_Selector_binder( pybind11::class_<HepMC3::StandardSelector, std::shared_ptr<HepMC3::StandardSelector>> cl);
	void	search_binder(pybind11::module &M);
} // namespace binder
#endif
