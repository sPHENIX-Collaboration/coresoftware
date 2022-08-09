#include <HepMC3/Attribute.h>
#include <HepMC3/Data/GenRunInfoData.h>
#include <HepMC3/GenCrossSection.h>
#include <HepMC3/GenHeavyIon.h>
#include <HepMC3/GenRunInfo.h>
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
#include <HepMC3/Reader.h>
#include <HepMC3/Writer.h>
#include <HepMC3/Print.h>
#include <src/stl_binders.hpp>
#include <src/binders.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

// HepMC3::GenCrossSection file:HepMC3/GenCrossSection.h line:42
struct PyCallBack_HepMC3_GenCrossSection : public HepMC3::GenCrossSection {
	using HepMC3::GenCrossSection::GenCrossSection;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::GenCrossSection *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return GenCrossSection::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::GenCrossSection *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return GenCrossSection::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::GenCrossSection *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return Attribute::init();
	}
	bool init(const class HepMC3::GenRunInfo & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::GenCrossSection *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return Attribute::init(a0);
	}
};

void bind_pyHepMC3_7(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // HepMC3::GenCrossSection file:HepMC3/GenCrossSection.h line:42
		pybind11::class_<HepMC3::GenCrossSection, std::shared_ptr<HepMC3::GenCrossSection>, PyCallBack_HepMC3_GenCrossSection, HepMC3::Attribute> cl(M("HepMC3"), "GenCrossSection", "");
		cl.def( pybind11::init( [](PyCallBack_HepMC3_GenCrossSection const &o){ return new PyCallBack_HepMC3_GenCrossSection(o); } ) );
		cl.def( pybind11::init( [](HepMC3::GenCrossSection const &o){ return new HepMC3::GenCrossSection(o); } ) );
		cl.def( pybind11::init( [](){ return new HepMC3::GenCrossSection(); }, [](){ return new PyCallBack_HepMC3_GenCrossSection(); } ) );
		cl.def("from_string", (bool (HepMC3::GenCrossSection::*)(const std::string &)) &HepMC3::GenCrossSection::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::GenCrossSection::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::GenCrossSection::*)(std::string &) const) &HepMC3::GenCrossSection::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::GenCrossSection::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("set_cross_section", [](HepMC3::GenCrossSection &o, const double & a0, const double & a1) -> void { return o.set_cross_section(a0, a1); }, "", pybind11::arg("xs"), pybind11::arg("xs_err"));
		cl.def("set_cross_section", [](HepMC3::GenCrossSection &o, const double & a0, const double & a1, const long & a2) -> void { return o.set_cross_section(a0, a1, a2); }, "", pybind11::arg("xs"), pybind11::arg("xs_err"), pybind11::arg("n_acc"));
		cl.def("set_cross_section", (void (HepMC3::GenCrossSection::*)(const double &, const double &, const long &, const long &)) &HepMC3::GenCrossSection::set_cross_section, "Set all fields \n\nC++: HepMC3::GenCrossSection::set_cross_section(const double &, const double &, const long &, const long &) --> void", pybind11::arg("xs"), pybind11::arg("xs_err"), pybind11::arg("n_acc"), pybind11::arg("n_att"));
		cl.def("set_accepted_events", (void (HepMC3::GenCrossSection::*)(const long &)) &HepMC3::GenCrossSection::set_accepted_events, "Set the number of accepted events\n\nC++: HepMC3::GenCrossSection::set_accepted_events(const long &) --> void", pybind11::arg("n_acc"));
		cl.def("set_attempted_events", (void (HepMC3::GenCrossSection::*)(const long &)) &HepMC3::GenCrossSection::set_attempted_events, "Set the number of attempted events\n\nC++: HepMC3::GenCrossSection::set_attempted_events(const long &) --> void", pybind11::arg("n_att"));
		cl.def("get_accepted_events", (long (HepMC3::GenCrossSection::*)() const) &HepMC3::GenCrossSection::get_accepted_events, "Get the number of accepted events\n\nC++: HepMC3::GenCrossSection::get_accepted_events() const --> long");
		cl.def("get_attempted_events", (long (HepMC3::GenCrossSection::*)() const) &HepMC3::GenCrossSection::get_attempted_events, "Get the number of attempted events\n\nC++: HepMC3::GenCrossSection::get_attempted_events() const --> long");
		cl.def("set_xsec", (void (HepMC3::GenCrossSection::*)(const std::string &, const double &)) &HepMC3::GenCrossSection::set_xsec, "Set the cross section  corresponding to the weight\n        named \n     \n\nC++: HepMC3::GenCrossSection::set_xsec(const std::string &, const double &) --> void", pybind11::arg("wName"), pybind11::arg("xs"));
		cl.def("set_xsec", (void (HepMC3::GenCrossSection::*)(const int &, const double &)) &HepMC3::GenCrossSection::set_xsec, "Set the cross section corresponding to the weight with\n        index \n     \n\nC++: HepMC3::GenCrossSection::set_xsec(const int &, const double &) --> void", pybind11::arg("indx"), pybind11::arg("xs"));
		cl.def("set_xsec_err", (void (HepMC3::GenCrossSection::*)(const std::string &, const double &)) &HepMC3::GenCrossSection::set_xsec_err, "Set the cross section error corresponding to the weight\n        named \n     \n\nC++: HepMC3::GenCrossSection::set_xsec_err(const std::string &, const double &) --> void", pybind11::arg("wName"), pybind11::arg("xs_err"));
		cl.def("set_xsec_err", (void (HepMC3::GenCrossSection::*)(const int &, const double &)) &HepMC3::GenCrossSection::set_xsec_err, "Set the cross section error corresponding to the weight\n        with index \n     \n\nC++: HepMC3::GenCrossSection::set_xsec_err(const int &, const double &) --> void", pybind11::arg("indx"), pybind11::arg("xs_err"));
		cl.def("xsec", (double (HepMC3::GenCrossSection::*)(const std::string &) const) &HepMC3::GenCrossSection::xsec, "Get the cross section corresponding to the weight named\n        \n     \n\nC++: HepMC3::GenCrossSection::xsec(const std::string &) const --> double", pybind11::arg("wName"));
		cl.def("xsec", [](HepMC3::GenCrossSection const &o) -> double { return o.xsec(); }, "");
		cl.def("xsec", (double (HepMC3::GenCrossSection::*)(const int &) const) &HepMC3::GenCrossSection::xsec, "Get the cross section corresponding to the weight with index\n        \n     \n\nC++: HepMC3::GenCrossSection::xsec(const int &) const --> double", pybind11::arg("indx"));
		cl.def("xsec_err", (double (HepMC3::GenCrossSection::*)(const std::string &) const) &HepMC3::GenCrossSection::xsec_err, "Get the cross section error corresponding to the weight\n        named \n     \n\nC++: HepMC3::GenCrossSection::xsec_err(const std::string &) const --> double", pybind11::arg("wName"));
		cl.def("xsec_err", [](HepMC3::GenCrossSection const &o) -> double { return o.xsec_err(); }, "");
		cl.def("xsec_err", (double (HepMC3::GenCrossSection::*)(const int &) const) &HepMC3::GenCrossSection::xsec_err, "Get the cross section error corresponding to the weight\n        with index \n     \n\nC++: HepMC3::GenCrossSection::xsec_err(const int &) const --> double", pybind11::arg("indx"));
		cl.def("__eq__", (bool (HepMC3::GenCrossSection::*)(const class HepMC3::GenCrossSection &) const) &HepMC3::GenCrossSection::operator==, "C++: HepMC3::GenCrossSection::operator==(const class HepMC3::GenCrossSection &) const --> bool", pybind11::arg(""));
		cl.def("__ne__", (bool (HepMC3::GenCrossSection::*)(const class HepMC3::GenCrossSection &) const) &HepMC3::GenCrossSection::operator!=, "C++: HepMC3::GenCrossSection::operator!=(const class HepMC3::GenCrossSection &) const --> bool", pybind11::arg(""));
		cl.def("is_valid", (bool (HepMC3::GenCrossSection::*)() const) &HepMC3::GenCrossSection::is_valid, "C++: HepMC3::GenCrossSection::is_valid() const --> bool");
		cl.def("assign", (class HepMC3::GenCrossSection & (HepMC3::GenCrossSection::*)(const class HepMC3::GenCrossSection &)) &HepMC3::GenCrossSection::operator=, "C++: HepMC3::GenCrossSection::operator=(const class HepMC3::GenCrossSection &) --> class HepMC3::GenCrossSection &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::Units file: line:26
		pybind11::class_<HepMC3::Units, std::shared_ptr<HepMC3::Units>> cl(M("HepMC3"), "Units", "");
		cl.def( pybind11::init( [](){ return new HepMC3::Units(); } ) );

		pybind11::enum_<HepMC3::Units::MomentumUnit>(cl, "MomentumUnit", pybind11::arithmetic(), "Momentum units ")
			.value("MEV", HepMC3::Units::MEV)
			.value("GEV", HepMC3::Units::GEV)
			.export_values();


		pybind11::enum_<HepMC3::Units::LengthUnit>(cl, "LengthUnit", pybind11::arithmetic(), "Position units ")
			.value("MM", HepMC3::Units::MM)
			.value("CM", HepMC3::Units::CM)
			.export_values();

		cl.def_static("momentum_unit", (enum HepMC3::Units::MomentumUnit (*)(const std::string &)) &HepMC3::Units::momentum_unit, "Get momentum unit based on its name\n\nC++: HepMC3::Units::momentum_unit(const std::string &) --> enum HepMC3::Units::MomentumUnit", pybind11::arg("name"));
		cl.def_static("length_unit", (enum HepMC3::Units::LengthUnit (*)(const std::string &)) &HepMC3::Units::length_unit, "Get length unit based on its name\n\nC++: HepMC3::Units::length_unit(const std::string &) --> enum HepMC3::Units::LengthUnit", pybind11::arg("name"));
		cl.def_static("name", (std::string (*)(enum HepMC3::Units::MomentumUnit)) &HepMC3::Units::name, "Get name of momentum unit \n\nC++: HepMC3::Units::name(enum HepMC3::Units::MomentumUnit) --> std::string", pybind11::arg("u"));
		cl.def_static("name", (std::string (*)(enum HepMC3::Units::LengthUnit)) &HepMC3::Units::name, "Get name of length unit \n\nC++: HepMC3::Units::name(enum HepMC3::Units::LengthUnit) --> std::string", pybind11::arg("u"));

		 binder::custom_Units_binder(cl);
	}
}
