#include <HepMC3/Attribute.h>
#include <HepMC3/Data/GenRunInfoData.h>
#include <HepMC3/FourVector.h>
#include <HepMC3/GenCrossSection.h>
#include <HepMC3/GenHeavyIon.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/LHEFAttributes.h>
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

// HepMC3::HEPRUPAttribute file:HepMC3/LHEFAttributes.h line:26
struct PyCallBack_HepMC3_HEPRUPAttribute : public HepMC3::HEPRUPAttribute {
	using HepMC3::HEPRUPAttribute::HEPRUPAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::HEPRUPAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HEPRUPAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::HEPRUPAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HEPRUPAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::HEPRUPAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::HEPRUPAttribute *>(this), "init");
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

// HepMC3::HEPEUPAttribute file:HepMC3/LHEFAttributes.h line:68
struct PyCallBack_HepMC3_HEPEUPAttribute : public HepMC3::HEPEUPAttribute {
	using HepMC3::HEPEUPAttribute::HEPEUPAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::HEPEUPAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HEPEUPAttribute::from_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::HEPEUPAttribute *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HEPEUPAttribute::init();
	}
	bool init(const class HepMC3::GenRunInfo & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::HEPEUPAttribute *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HEPEUPAttribute::init(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::HEPEUPAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HEPEUPAttribute::to_string(a0);
	}
};

void bind_pyHepMC3_17(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // HepMC3::HEPRUPAttribute file:HepMC3/LHEFAttributes.h line:26
		pybind11::class_<HepMC3::HEPRUPAttribute, std::shared_ptr<HepMC3::HEPRUPAttribute>, PyCallBack_HepMC3_HEPRUPAttribute, HepMC3::Attribute> cl(M("HepMC3"), "HEPRUPAttribute", "Class for storing data for LHEF run information");
		cl.def( pybind11::init( [](){ return new HepMC3::HEPRUPAttribute(); }, [](){ return new PyCallBack_HepMC3_HEPRUPAttribute(); } ) );
		cl.def( pybind11::init<std::string>(), pybind11::arg("s") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_HEPRUPAttribute const &o){ return new PyCallBack_HepMC3_HEPRUPAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::HEPRUPAttribute const &o){ return new HepMC3::HEPRUPAttribute(o); } ) );
		cl.def_readwrite("heprup", &HepMC3::HEPRUPAttribute::heprup);
		cl.def_readwrite("tags", &HepMC3::HEPRUPAttribute::tags);
		cl.def("from_string", (bool (HepMC3::HEPRUPAttribute::*)(const std::string &)) &HepMC3::HEPRUPAttribute::from_string, "Fill class content from string \n\nC++: HepMC3::HEPRUPAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::HEPRUPAttribute::*)(std::string &) const) &HepMC3::HEPRUPAttribute::to_string, "Fill string from class content \n\nC++: HepMC3::HEPRUPAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("clear", (void (HepMC3::HEPRUPAttribute::*)()) &HepMC3::HEPRUPAttribute::clear, "Clear this object. \n\nC++: HepMC3::HEPRUPAttribute::clear() --> void");
		cl.def("assign", (class HepMC3::HEPRUPAttribute & (HepMC3::HEPRUPAttribute::*)(const class HepMC3::HEPRUPAttribute &)) &HepMC3::HEPRUPAttribute::operator=, "C++: HepMC3::HEPRUPAttribute::operator=(const class HepMC3::HEPRUPAttribute &) --> class HepMC3::HEPRUPAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::HEPEUPAttribute file:HepMC3/LHEFAttributes.h line:68
		pybind11::class_<HepMC3::HEPEUPAttribute, std::shared_ptr<HepMC3::HEPEUPAttribute>, PyCallBack_HepMC3_HEPEUPAttribute, HepMC3::Attribute> cl(M("HepMC3"), "HEPEUPAttribute", "Class for storing data for LHEF run information");
		cl.def( pybind11::init( [](){ return new HepMC3::HEPEUPAttribute(); }, [](){ return new PyCallBack_HepMC3_HEPEUPAttribute(); } ) );
		cl.def( pybind11::init<std::string>(), pybind11::arg("s") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_HEPEUPAttribute const &o){ return new PyCallBack_HepMC3_HEPEUPAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::HEPEUPAttribute const &o){ return new HepMC3::HEPEUPAttribute(o); } ) );
		cl.def_readwrite("hepeup", &HepMC3::HEPEUPAttribute::hepeup);
		cl.def_readwrite("tags", &HepMC3::HEPEUPAttribute::tags);
		cl.def("from_string", (bool (HepMC3::HEPEUPAttribute::*)(const std::string &)) &HepMC3::HEPEUPAttribute::from_string, "Fill class content from string \n\nC++: HepMC3::HEPEUPAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("init", (bool (HepMC3::HEPEUPAttribute::*)()) &HepMC3::HEPEUPAttribute::init, "Parse the XML-tags. \n\nC++: HepMC3::HEPEUPAttribute::init() --> bool");
		cl.def("init", (bool (HepMC3::HEPEUPAttribute::*)(const class HepMC3::GenRunInfo &)) &HepMC3::HEPEUPAttribute::init, "Dummy function. \n\nC++: HepMC3::HEPEUPAttribute::init(const class HepMC3::GenRunInfo &) --> bool", pybind11::arg(""));
		cl.def("to_string", (bool (HepMC3::HEPEUPAttribute::*)(std::string &) const) &HepMC3::HEPEUPAttribute::to_string, "Fill string from class content \n\nC++: HepMC3::HEPEUPAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("momentum", (class HepMC3::FourVector (HepMC3::HEPEUPAttribute::*)(int) const) &HepMC3::HEPEUPAttribute::momentum, "Get momentum \n\nC++: HepMC3::HEPEUPAttribute::momentum(int) const --> class HepMC3::FourVector", pybind11::arg("i"));
		cl.def("clear", (void (HepMC3::HEPEUPAttribute::*)()) &HepMC3::HEPEUPAttribute::clear, "Clear this object. \n\nC++: HepMC3::HEPEUPAttribute::clear() --> void");
		cl.def("assign", (class HepMC3::HEPEUPAttribute & (HepMC3::HEPEUPAttribute::*)(const class HepMC3::HEPEUPAttribute &)) &HepMC3::HEPEUPAttribute::operator=, "C++: HepMC3::HEPEUPAttribute::operator=(const class HepMC3::HEPEUPAttribute &) --> class HepMC3::HEPEUPAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
