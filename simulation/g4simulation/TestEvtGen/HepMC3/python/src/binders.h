#ifndef BINDERS_H
#define BINDERS_H

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenHeavyIon.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/LHEF.h>
#include <HepMC3/HEPEVT_Wrapper_Runtime.h>
#include <pybind11/pybind11.h>
namespace binder {
void custom_HEPEVT_Wrapper_Runtime_binder(pybind11::class_<HepMC3::HEPEVT_Wrapper_Runtime, std::shared_ptr<HepMC3::HEPEVT_Wrapper_Runtime>> cl);
void custom_GenEvent_binder(pybind11::class_<HepMC3::GenEvent, std::shared_ptr<HepMC3::GenEvent>> cl);
void custom_GenParticle_binder(pybind11::class_<HepMC3::GenParticle, std::shared_ptr<HepMC3::GenParticle>> cl);
void custom_GenVertex_binder(pybind11::class_<HepMC3::GenVertex, std::shared_ptr<HepMC3::GenVertex>> cl);

void custom_GenRunInfo_binder(pybind11::class_<HepMC3::GenRunInfo, std::shared_ptr<HepMC3::GenRunInfo>> cl);
void custom_Units_binder(pybind11::class_<HepMC3::Units, std::shared_ptr<HepMC3::Units>> cl);

void custom_FourVector_binder(pybind11::class_<HepMC3::FourVector, std::shared_ptr<HepMC3::FourVector>> cl);
template <typename T>  void custom_T_binder (pybind11::class_<T, std::shared_ptr<T>> cl)
{
//cl.def("print", (void (T::*)(std::ostream &) const) &T::print, "Print the object", pybind11::arg("file"));
cl.def("print", [](T const &o, pybind11::object  & a1) -> void { std::stringstream b;  o.print(b); a1.attr("write")(pybind11::str(b.str().c_str())); }, "Print the object", pybind11::arg("file"));
}
void custom_LHEFTagBase_binder (pybind11::class_<LHEF::TagBase, std::shared_ptr<LHEF::TagBase>> cl);
void	print_binder(pybind11::module &M);

} // namespace binder

#endif
