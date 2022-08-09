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

// HepMC3::VectorULongLongAttribute file:HepMC3/Attribute.h line:956
struct PyCallBack_HepMC3_VectorULongLongAttribute : public HepMC3::VectorULongLongAttribute {
	using HepMC3::VectorULongLongAttribute::VectorULongLongAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorULongLongAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorULongLongAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorULongLongAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorULongLongAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorULongLongAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorULongLongAttribute *>(this), "init");
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

// HepMC3::VectorIntAttribute file:HepMC3/Attribute.h line:1001
struct PyCallBack_HepMC3_VectorIntAttribute : public HepMC3::VectorIntAttribute {
	using HepMC3::VectorIntAttribute::VectorIntAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorIntAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorIntAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorIntAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorIntAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorIntAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorIntAttribute *>(this), "init");
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

// HepMC3::VectorLongIntAttribute file:HepMC3/Attribute.h line:1046
struct PyCallBack_HepMC3_VectorLongIntAttribute : public HepMC3::VectorLongIntAttribute {
	using HepMC3::VectorLongIntAttribute::VectorLongIntAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorLongIntAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorLongIntAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorLongIntAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorLongIntAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorLongIntAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorLongIntAttribute *>(this), "init");
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

// HepMC3::VectorDoubleAttribute file:HepMC3/Attribute.h line:1091
struct PyCallBack_HepMC3_VectorDoubleAttribute : public HepMC3::VectorDoubleAttribute {
	using HepMC3::VectorDoubleAttribute::VectorDoubleAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorDoubleAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorDoubleAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorDoubleAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorDoubleAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorDoubleAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorDoubleAttribute *>(this), "init");
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

// HepMC3::VectorStringAttribute file:HepMC3/Attribute.h line:1137
struct PyCallBack_HepMC3_VectorStringAttribute : public HepMC3::VectorStringAttribute {
	using HepMC3::VectorStringAttribute::VectorStringAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorStringAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorStringAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorStringAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorStringAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorStringAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorStringAttribute *>(this), "init");
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

// HepMC3::GenHeavyIon file:HepMC3/GenHeavyIon.h line:28
struct PyCallBack_HepMC3_GenHeavyIon : public HepMC3::GenHeavyIon {
	using HepMC3::GenHeavyIon::GenHeavyIon;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::GenHeavyIon *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return GenHeavyIon::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::GenHeavyIon *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return GenHeavyIon::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::GenHeavyIon *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::GenHeavyIon *>(this), "init");
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

void bind_pyHepMC3_6(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // HepMC3::VectorULongLongAttribute file:HepMC3/Attribute.h line:956
		pybind11::class_<HepMC3::VectorULongLongAttribute, std::shared_ptr<HepMC3::VectorULongLongAttribute>, PyCallBack_HepMC3_VectorULongLongAttribute, HepMC3::Attribute> cl(M("HepMC3"), "VectorULongLongAttribute", "Attribute that holds a vector of unsigned long longegers of type  unsigned long long\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::VectorULongLongAttribute(); }, [](){ return new PyCallBack_HepMC3_VectorULongLongAttribute(); } ) );
		cl.def( pybind11::init<class std::vector<unsigned long long>>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_VectorULongLongAttribute const &o){ return new PyCallBack_HepMC3_VectorULongLongAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::VectorULongLongAttribute const &o){ return new HepMC3::VectorULongLongAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::VectorULongLongAttribute::*)(const std::string &)) &HepMC3::VectorULongLongAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::VectorULongLongAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::VectorULongLongAttribute::*)(std::string &) const) &HepMC3::VectorULongLongAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::VectorULongLongAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (class std::vector<unsigned long long> (HepMC3::VectorULongLongAttribute::*)() const) &HepMC3::VectorULongLongAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::VectorULongLongAttribute::value() const --> class std::vector<unsigned long long>");
		cl.def("set_value", (void (HepMC3::VectorULongLongAttribute::*)(const class std::vector<unsigned long long> &)) &HepMC3::VectorULongLongAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::VectorULongLongAttribute::set_value(const class std::vector<unsigned long long> &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::VectorULongLongAttribute & (HepMC3::VectorULongLongAttribute::*)(const class HepMC3::VectorULongLongAttribute &)) &HepMC3::VectorULongLongAttribute::operator=, "C++: HepMC3::VectorULongLongAttribute::operator=(const class HepMC3::VectorULongLongAttribute &) --> class HepMC3::VectorULongLongAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::VectorIntAttribute file:HepMC3/Attribute.h line:1001
		pybind11::class_<HepMC3::VectorIntAttribute, std::shared_ptr<HepMC3::VectorIntAttribute>, PyCallBack_HepMC3_VectorIntAttribute, HepMC3::Attribute> cl(M("HepMC3"), "VectorIntAttribute", "Attribute that holds a vector of integers of type  int\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::VectorIntAttribute(); }, [](){ return new PyCallBack_HepMC3_VectorIntAttribute(); } ) );
		cl.def( pybind11::init<class std::vector<int>>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_VectorIntAttribute const &o){ return new PyCallBack_HepMC3_VectorIntAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::VectorIntAttribute const &o){ return new HepMC3::VectorIntAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::VectorIntAttribute::*)(const std::string &)) &HepMC3::VectorIntAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::VectorIntAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::VectorIntAttribute::*)(std::string &) const) &HepMC3::VectorIntAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::VectorIntAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (class std::vector<int> (HepMC3::VectorIntAttribute::*)() const) &HepMC3::VectorIntAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::VectorIntAttribute::value() const --> class std::vector<int>");
		cl.def("set_value", (void (HepMC3::VectorIntAttribute::*)(const class std::vector<int> &)) &HepMC3::VectorIntAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::VectorIntAttribute::set_value(const class std::vector<int> &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::VectorIntAttribute & (HepMC3::VectorIntAttribute::*)(const class HepMC3::VectorIntAttribute &)) &HepMC3::VectorIntAttribute::operator=, "C++: HepMC3::VectorIntAttribute::operator=(const class HepMC3::VectorIntAttribute &) --> class HepMC3::VectorIntAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::VectorLongIntAttribute file:HepMC3/Attribute.h line:1046
		pybind11::class_<HepMC3::VectorLongIntAttribute, std::shared_ptr<HepMC3::VectorLongIntAttribute>, PyCallBack_HepMC3_VectorLongIntAttribute, HepMC3::Attribute> cl(M("HepMC3"), "VectorLongIntAttribute", "Attribute that holds a vector of integers of type  int\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::VectorLongIntAttribute(); }, [](){ return new PyCallBack_HepMC3_VectorLongIntAttribute(); } ) );
		cl.def( pybind11::init<class std::vector<long>>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_VectorLongIntAttribute const &o){ return new PyCallBack_HepMC3_VectorLongIntAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::VectorLongIntAttribute const &o){ return new HepMC3::VectorLongIntAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::VectorLongIntAttribute::*)(const std::string &)) &HepMC3::VectorLongIntAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::VectorLongIntAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::VectorLongIntAttribute::*)(std::string &) const) &HepMC3::VectorLongIntAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::VectorLongIntAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (class std::vector<long> (HepMC3::VectorLongIntAttribute::*)() const) &HepMC3::VectorLongIntAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::VectorLongIntAttribute::value() const --> class std::vector<long>");
		cl.def("set_value", (void (HepMC3::VectorLongIntAttribute::*)(const class std::vector<long> &)) &HepMC3::VectorLongIntAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::VectorLongIntAttribute::set_value(const class std::vector<long> &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::VectorLongIntAttribute & (HepMC3::VectorLongIntAttribute::*)(const class HepMC3::VectorLongIntAttribute &)) &HepMC3::VectorLongIntAttribute::operator=, "C++: HepMC3::VectorLongIntAttribute::operator=(const class HepMC3::VectorLongIntAttribute &) --> class HepMC3::VectorLongIntAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::VectorDoubleAttribute file:HepMC3/Attribute.h line:1091
		pybind11::class_<HepMC3::VectorDoubleAttribute, std::shared_ptr<HepMC3::VectorDoubleAttribute>, PyCallBack_HepMC3_VectorDoubleAttribute, HepMC3::Attribute> cl(M("HepMC3"), "VectorDoubleAttribute", "Attribute that holds a vector of FPs of type  double\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::VectorDoubleAttribute(); }, [](){ return new PyCallBack_HepMC3_VectorDoubleAttribute(); } ) );
		cl.def( pybind11::init<class std::vector<double>>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_VectorDoubleAttribute const &o){ return new PyCallBack_HepMC3_VectorDoubleAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::VectorDoubleAttribute const &o){ return new HepMC3::VectorDoubleAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::VectorDoubleAttribute::*)(const std::string &)) &HepMC3::VectorDoubleAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::VectorDoubleAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::VectorDoubleAttribute::*)(std::string &) const) &HepMC3::VectorDoubleAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::VectorDoubleAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (class std::vector<double> (HepMC3::VectorDoubleAttribute::*)() const) &HepMC3::VectorDoubleAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::VectorDoubleAttribute::value() const --> class std::vector<double>");
		cl.def("set_value", (void (HepMC3::VectorDoubleAttribute::*)(const class std::vector<double> &)) &HepMC3::VectorDoubleAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::VectorDoubleAttribute::set_value(const class std::vector<double> &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::VectorDoubleAttribute & (HepMC3::VectorDoubleAttribute::*)(const class HepMC3::VectorDoubleAttribute &)) &HepMC3::VectorDoubleAttribute::operator=, "C++: HepMC3::VectorDoubleAttribute::operator=(const class HepMC3::VectorDoubleAttribute &) --> class HepMC3::VectorDoubleAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::VectorStringAttribute file:HepMC3/Attribute.h line:1137
		pybind11::class_<HepMC3::VectorStringAttribute, std::shared_ptr<HepMC3::VectorStringAttribute>, PyCallBack_HepMC3_VectorStringAttribute, HepMC3::Attribute> cl(M("HepMC3"), "VectorStringAttribute", "Attribute that holds a vector of FPs of type  string\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::VectorStringAttribute(); }, [](){ return new PyCallBack_HepMC3_VectorStringAttribute(); } ) );
		cl.def( pybind11::init<class std::vector<std::string >>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_VectorStringAttribute const &o){ return new PyCallBack_HepMC3_VectorStringAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::VectorStringAttribute const &o){ return new HepMC3::VectorStringAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::VectorStringAttribute::*)(const std::string &)) &HepMC3::VectorStringAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::VectorStringAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::VectorStringAttribute::*)(std::string &) const) &HepMC3::VectorStringAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::VectorStringAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (class std::vector<std::string > (HepMC3::VectorStringAttribute::*)() const) &HepMC3::VectorStringAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::VectorStringAttribute::value() const --> class std::vector<std::string >");
		cl.def("set_value", (void (HepMC3::VectorStringAttribute::*)(const class std::vector<std::string > &)) &HepMC3::VectorStringAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::VectorStringAttribute::set_value(const class std::vector<std::string > &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::VectorStringAttribute & (HepMC3::VectorStringAttribute::*)(const class HepMC3::VectorStringAttribute &)) &HepMC3::VectorStringAttribute::operator=, "C++: HepMC3::VectorStringAttribute::operator=(const class HepMC3::VectorStringAttribute &) --> class HepMC3::VectorStringAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::GenHeavyIon file:HepMC3/GenHeavyIon.h line:28
		pybind11::class_<HepMC3::GenHeavyIon, std::shared_ptr<HepMC3::GenHeavyIon>, PyCallBack_HepMC3_GenHeavyIon, HepMC3::Attribute> cl(M("HepMC3"), "GenHeavyIon", "");
		cl.def( pybind11::init( [](){ return new HepMC3::GenHeavyIon(); }, [](){ return new PyCallBack_HepMC3_GenHeavyIon(); } ) );
		cl.def( pybind11::init( [](PyCallBack_HepMC3_GenHeavyIon const &o){ return new PyCallBack_HepMC3_GenHeavyIon(o); } ) );
		cl.def( pybind11::init( [](HepMC3::GenHeavyIon const &o){ return new HepMC3::GenHeavyIon(o); } ) );
		cl.def_readwrite("Ncoll_hard", &HepMC3::GenHeavyIon::Ncoll_hard);
		cl.def_readwrite("Npart_proj", &HepMC3::GenHeavyIon::Npart_proj);
		cl.def_readwrite("Npart_targ", &HepMC3::GenHeavyIon::Npart_targ);
		cl.def_readwrite("Ncoll", &HepMC3::GenHeavyIon::Ncoll);
		cl.def_readwrite("spectator_neutrons", &HepMC3::GenHeavyIon::spectator_neutrons);
		cl.def_readwrite("spectator_protons", &HepMC3::GenHeavyIon::spectator_protons);
		cl.def_readwrite("N_Nwounded_collisions", &HepMC3::GenHeavyIon::N_Nwounded_collisions);
		cl.def_readwrite("Nwounded_N_collisions", &HepMC3::GenHeavyIon::Nwounded_N_collisions);
		cl.def_readwrite("Nwounded_Nwounded_collisions", &HepMC3::GenHeavyIon::Nwounded_Nwounded_collisions);
		cl.def_readwrite("impact_parameter", &HepMC3::GenHeavyIon::impact_parameter);
		cl.def_readwrite("event_plane_angle", &HepMC3::GenHeavyIon::event_plane_angle);
		cl.def_readwrite("eccentricity", &HepMC3::GenHeavyIon::eccentricity);
		cl.def_readwrite("sigma_inel_NN", &HepMC3::GenHeavyIon::sigma_inel_NN);
		cl.def_readwrite("centrality", &HepMC3::GenHeavyIon::centrality);
		cl.def_readwrite("user_cent_estimate", &HepMC3::GenHeavyIon::user_cent_estimate);
		cl.def_readwrite("Nspec_proj_n", &HepMC3::GenHeavyIon::Nspec_proj_n);
		cl.def_readwrite("Nspec_targ_n", &HepMC3::GenHeavyIon::Nspec_targ_n);
		cl.def_readwrite("Nspec_proj_p", &HepMC3::GenHeavyIon::Nspec_proj_p);
		cl.def_readwrite("Nspec_targ_p", &HepMC3::GenHeavyIon::Nspec_targ_p);
		cl.def_readwrite("participant_plane_angles", &HepMC3::GenHeavyIon::participant_plane_angles);
		cl.def_readwrite("eccentricities", &HepMC3::GenHeavyIon::eccentricities);
		cl.def_readwrite("forceoldformat", &HepMC3::GenHeavyIon::forceoldformat);
		cl.def("from_string", (bool (HepMC3::GenHeavyIon::*)(const std::string &)) &HepMC3::GenHeavyIon::from_string, "Implementation of Attribute::from_string.\n\nC++: HepMC3::GenHeavyIon::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::GenHeavyIon::*)(std::string &) const) &HepMC3::GenHeavyIon::to_string, "Implementation of Attribute::to_string.\n\nC++: HepMC3::GenHeavyIon::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("__eq__", (bool (HepMC3::GenHeavyIon::*)(const class HepMC3::GenHeavyIon &) const) &HepMC3::GenHeavyIon::operator==, "Operator ==\n\nC++: HepMC3::GenHeavyIon::operator==(const class HepMC3::GenHeavyIon &) const --> bool", pybind11::arg(""));
		cl.def("__ne__", (bool (HepMC3::GenHeavyIon::*)(const class HepMC3::GenHeavyIon &) const) &HepMC3::GenHeavyIon::operator!=, "Operator !=\n\nC++: HepMC3::GenHeavyIon::operator!=(const class HepMC3::GenHeavyIon &) const --> bool", pybind11::arg(""));
		cl.def("set", [](HepMC3::GenHeavyIon &o, const int & a0, const int & a1, const int & a2, const int & a3, const int & a4, const int & a5) -> void { return o.set(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("nh"), pybind11::arg("np"), pybind11::arg("nt"), pybind11::arg("nc"), pybind11::arg("ns"), pybind11::arg("nsp"));
		cl.def("set", [](HepMC3::GenHeavyIon &o, const int & a0, const int & a1, const int & a2, const int & a3, const int & a4, const int & a5, const int & a6) -> void { return o.set(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("nh"), pybind11::arg("np"), pybind11::arg("nt"), pybind11::arg("nc"), pybind11::arg("ns"), pybind11::arg("nsp"), pybind11::arg("nnw"));
		cl.def("set", [](HepMC3::GenHeavyIon &o, const int & a0, const int & a1, const int & a2, const int & a3, const int & a4, const int & a5, const int & a6, const int & a7) -> void { return o.set(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("nh"), pybind11::arg("np"), pybind11::arg("nt"), pybind11::arg("nc"), pybind11::arg("ns"), pybind11::arg("nsp"), pybind11::arg("nnw"), pybind11::arg("nwn"));
		cl.def("set", [](HepMC3::GenHeavyIon &o, const int & a0, const int & a1, const int & a2, const int & a3, const int & a4, const int & a5, const int & a6, const int & a7, const int & a8) -> void { return o.set(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("nh"), pybind11::arg("np"), pybind11::arg("nt"), pybind11::arg("nc"), pybind11::arg("ns"), pybind11::arg("nsp"), pybind11::arg("nnw"), pybind11::arg("nwn"), pybind11::arg("nwnw"));
		cl.def("set", [](HepMC3::GenHeavyIon &o, const int & a0, const int & a1, const int & a2, const int & a3, const int & a4, const int & a5, const int & a6, const int & a7, const int & a8, const double & a9) -> void { return o.set(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("nh"), pybind11::arg("np"), pybind11::arg("nt"), pybind11::arg("nc"), pybind11::arg("ns"), pybind11::arg("nsp"), pybind11::arg("nnw"), pybind11::arg("nwn"), pybind11::arg("nwnw"), pybind11::arg("im"));
		cl.def("set", [](HepMC3::GenHeavyIon &o, const int & a0, const int & a1, const int & a2, const int & a3, const int & a4, const int & a5, const int & a6, const int & a7, const int & a8, const double & a9, const double & a10) -> void { return o.set(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("nh"), pybind11::arg("np"), pybind11::arg("nt"), pybind11::arg("nc"), pybind11::arg("ns"), pybind11::arg("nsp"), pybind11::arg("nnw"), pybind11::arg("nwn"), pybind11::arg("nwnw"), pybind11::arg("im"), pybind11::arg("pl"));
		cl.def("set", [](HepMC3::GenHeavyIon &o, const int & a0, const int & a1, const int & a2, const int & a3, const int & a4, const int & a5, const int & a6, const int & a7, const int & a8, const double & a9, const double & a10, const double & a11) -> void { return o.set(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("nh"), pybind11::arg("np"), pybind11::arg("nt"), pybind11::arg("nc"), pybind11::arg("ns"), pybind11::arg("nsp"), pybind11::arg("nnw"), pybind11::arg("nwn"), pybind11::arg("nwnw"), pybind11::arg("im"), pybind11::arg("pl"), pybind11::arg("ec"));
		cl.def("set", [](HepMC3::GenHeavyIon &o, const int & a0, const int & a1, const int & a2, const int & a3, const int & a4, const int & a5, const int & a6, const int & a7, const int & a8, const double & a9, const double & a10, const double & a11, const double & a12) -> void { return o.set(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); }, "", pybind11::arg("nh"), pybind11::arg("np"), pybind11::arg("nt"), pybind11::arg("nc"), pybind11::arg("ns"), pybind11::arg("nsp"), pybind11::arg("nnw"), pybind11::arg("nwn"), pybind11::arg("nwnw"), pybind11::arg("im"), pybind11::arg("pl"), pybind11::arg("ec"), pybind11::arg("s"));
		cl.def("set", [](HepMC3::GenHeavyIon &o, const int & a0, const int & a1, const int & a2, const int & a3, const int & a4, const int & a5, const int & a6, const int & a7, const int & a8, const double & a9, const double & a10, const double & a11, const double & a12, const double & a13) -> void { return o.set(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("nh"), pybind11::arg("np"), pybind11::arg("nt"), pybind11::arg("nc"), pybind11::arg("ns"), pybind11::arg("nsp"), pybind11::arg("nnw"), pybind11::arg("nwn"), pybind11::arg("nwnw"), pybind11::arg("im"), pybind11::arg("pl"), pybind11::arg("ec"), pybind11::arg("s"), pybind11::arg("cent"));
		cl.def("set", (void (HepMC3::GenHeavyIon::*)(const int &, const int &, const int &, const int &, const int &, const int &, const int &, const int &, const int &, const double &, const double &, const double &, const double &, const double &, const double &)) &HepMC3::GenHeavyIon::set, "Set all fields.\n\n HEPMC3_DEPRECATED(\"Set individual fields directly instead.\")\n \n\n Set all fields \n\nC++: HepMC3::GenHeavyIon::set(const int &, const int &, const int &, const int &, const int &, const int &, const int &, const int &, const int &, const double &, const double &, const double &, const double &, const double &, const double &) --> void", pybind11::arg("nh"), pybind11::arg("np"), pybind11::arg("nt"), pybind11::arg("nc"), pybind11::arg("ns"), pybind11::arg("nsp"), pybind11::arg("nnw"), pybind11::arg("nwn"), pybind11::arg("nwnw"), pybind11::arg("im"), pybind11::arg("pl"), pybind11::arg("ec"), pybind11::arg("s"), pybind11::arg("cent"), pybind11::arg("ucent"));
		cl.def("is_valid", (bool (HepMC3::GenHeavyIon::*)() const) &HepMC3::GenHeavyIon::is_valid, "Verify that the instance contains non-zero information.\n\n HEPMC3_DEPRECATED(\"Each filed now have default values meaning\n that they have not been set\")\n\nC++: HepMC3::GenHeavyIon::is_valid() const --> bool");
		cl.def("assign", (class HepMC3::GenHeavyIon & (HepMC3::GenHeavyIon::*)(const class HepMC3::GenHeavyIon &)) &HepMC3::GenHeavyIon::operator=, "C++: HepMC3::GenHeavyIon::operator=(const class HepMC3::GenHeavyIon &) --> class HepMC3::GenHeavyIon &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
