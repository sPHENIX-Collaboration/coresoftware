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

// HepMC3::ULongAttribute file:HepMC3/Attribute.h line:557
struct PyCallBack_HepMC3_ULongAttribute : public HepMC3::ULongAttribute {
	using HepMC3::ULongAttribute::ULongAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ULongAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ULongAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ULongAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ULongAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ULongAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ULongAttribute *>(this), "init");
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

// HepMC3::ULongLongAttribute file:HepMC3/Attribute.h line:599
struct PyCallBack_HepMC3_ULongLongAttribute : public HepMC3::ULongLongAttribute {
	using HepMC3::ULongLongAttribute::ULongLongAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ULongLongAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ULongLongAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ULongLongAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ULongLongAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ULongLongAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ULongLongAttribute *>(this), "init");
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

// HepMC3::BoolAttribute file:HepMC3/Attribute.h line:639
struct PyCallBack_HepMC3_BoolAttribute : public HepMC3::BoolAttribute {
	using HepMC3::BoolAttribute::BoolAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::BoolAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return BoolAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::BoolAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return BoolAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::BoolAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::BoolAttribute *>(this), "init");
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

// HepMC3::VectorCharAttribute file:HepMC3/Attribute.h line:682
struct PyCallBack_HepMC3_VectorCharAttribute : public HepMC3::VectorCharAttribute {
	using HepMC3::VectorCharAttribute::VectorCharAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorCharAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorCharAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorCharAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorCharAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorCharAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorCharAttribute *>(this), "init");
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

// HepMC3::VectorFloatAttribute file:HepMC3/Attribute.h line:727
struct PyCallBack_HepMC3_VectorFloatAttribute : public HepMC3::VectorFloatAttribute {
	using HepMC3::VectorFloatAttribute::VectorFloatAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorFloatAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorFloatAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorFloatAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorFloatAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorFloatAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorFloatAttribute *>(this), "init");
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

// HepMC3::VectorLongDoubleAttribute file:HepMC3/Attribute.h line:773
struct PyCallBack_HepMC3_VectorLongDoubleAttribute : public HepMC3::VectorLongDoubleAttribute {
	using HepMC3::VectorLongDoubleAttribute::VectorLongDoubleAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorLongDoubleAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorLongDoubleAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorLongDoubleAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorLongDoubleAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorLongDoubleAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorLongDoubleAttribute *>(this), "init");
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

// HepMC3::VectorLongLongAttribute file:HepMC3/Attribute.h line:820
struct PyCallBack_HepMC3_VectorLongLongAttribute : public HepMC3::VectorLongLongAttribute {
	using HepMC3::VectorLongLongAttribute::VectorLongLongAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorLongLongAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorLongLongAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorLongLongAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorLongLongAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorLongLongAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorLongLongAttribute *>(this), "init");
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

// HepMC3::VectorUIntAttribute file:HepMC3/Attribute.h line:865
struct PyCallBack_HepMC3_VectorUIntAttribute : public HepMC3::VectorUIntAttribute {
	using HepMC3::VectorUIntAttribute::VectorUIntAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorUIntAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorUIntAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorUIntAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorUIntAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorUIntAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorUIntAttribute *>(this), "init");
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

// HepMC3::VectorULongAttribute file:HepMC3/Attribute.h line:910
struct PyCallBack_HepMC3_VectorULongAttribute : public HepMC3::VectorULongAttribute {
	using HepMC3::VectorULongAttribute::VectorULongAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorULongAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorULongAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorULongAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return VectorULongAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorULongAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::VectorULongAttribute *>(this), "init");
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

void bind_pyHepMC3_5(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // HepMC3::ULongAttribute file:HepMC3/Attribute.h line:557
		pybind11::class_<HepMC3::ULongAttribute, std::shared_ptr<HepMC3::ULongAttribute>, PyCallBack_HepMC3_ULongAttribute, HepMC3::Attribute> cl(M("HepMC3"), "ULongAttribute", "Attribute that holds an unsigned long\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::ULongAttribute(); }, [](){ return new PyCallBack_HepMC3_ULongAttribute(); } ) );
		cl.def( pybind11::init<unsigned long>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_ULongAttribute const &o){ return new PyCallBack_HepMC3_ULongAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::ULongAttribute const &o){ return new HepMC3::ULongAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::ULongAttribute::*)(const std::string &)) &HepMC3::ULongAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::ULongAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::ULongAttribute::*)(std::string &) const) &HepMC3::ULongAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::ULongAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (unsigned long (HepMC3::ULongAttribute::*)() const) &HepMC3::ULongAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::ULongAttribute::value() const --> unsigned long");
		cl.def("set_value", (void (HepMC3::ULongAttribute::*)(const unsigned long &)) &HepMC3::ULongAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::ULongAttribute::set_value(const unsigned long &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::ULongAttribute & (HepMC3::ULongAttribute::*)(const class HepMC3::ULongAttribute &)) &HepMC3::ULongAttribute::operator=, "C++: HepMC3::ULongAttribute::operator=(const class HepMC3::ULongAttribute &) --> class HepMC3::ULongAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::ULongLongAttribute file:HepMC3/Attribute.h line:599
		pybind11::class_<HepMC3::ULongLongAttribute, std::shared_ptr<HepMC3::ULongLongAttribute>, PyCallBack_HepMC3_ULongLongAttribute, HepMC3::Attribute> cl(M("HepMC3"), "ULongLongAttribute", "Attribute that holds an unsigned long long\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::ULongLongAttribute(); }, [](){ return new PyCallBack_HepMC3_ULongLongAttribute(); } ) );
		cl.def( pybind11::init<unsigned long long>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_ULongLongAttribute const &o){ return new PyCallBack_HepMC3_ULongLongAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::ULongLongAttribute const &o){ return new HepMC3::ULongLongAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::ULongLongAttribute::*)(const std::string &)) &HepMC3::ULongLongAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::ULongLongAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::ULongLongAttribute::*)(std::string &) const) &HepMC3::ULongLongAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::ULongLongAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (unsigned long long (HepMC3::ULongLongAttribute::*)() const) &HepMC3::ULongLongAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::ULongLongAttribute::value() const --> unsigned long long");
		cl.def("set_value", (void (HepMC3::ULongLongAttribute::*)(const unsigned long long &)) &HepMC3::ULongLongAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::ULongLongAttribute::set_value(const unsigned long long &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::ULongLongAttribute & (HepMC3::ULongLongAttribute::*)(const class HepMC3::ULongLongAttribute &)) &HepMC3::ULongLongAttribute::operator=, "C++: HepMC3::ULongLongAttribute::operator=(const class HepMC3::ULongLongAttribute &) --> class HepMC3::ULongLongAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::BoolAttribute file:HepMC3/Attribute.h line:639
		pybind11::class_<HepMC3::BoolAttribute, std::shared_ptr<HepMC3::BoolAttribute>, PyCallBack_HepMC3_BoolAttribute, HepMC3::Attribute> cl(M("HepMC3"), "BoolAttribute", "Attribute that holds an Booleger implemented as an int\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::BoolAttribute(); }, [](){ return new PyCallBack_HepMC3_BoolAttribute(); } ) );
		cl.def( pybind11::init<bool>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_BoolAttribute const &o){ return new PyCallBack_HepMC3_BoolAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::BoolAttribute const &o){ return new HepMC3::BoolAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::BoolAttribute::*)(const std::string &)) &HepMC3::BoolAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::BoolAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::BoolAttribute::*)(std::string &) const) &HepMC3::BoolAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::BoolAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (bool (HepMC3::BoolAttribute::*)() const) &HepMC3::BoolAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::BoolAttribute::value() const --> bool");
		cl.def("set_value", (void (HepMC3::BoolAttribute::*)(const bool &)) &HepMC3::BoolAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::BoolAttribute::set_value(const bool &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::BoolAttribute & (HepMC3::BoolAttribute::*)(const class HepMC3::BoolAttribute &)) &HepMC3::BoolAttribute::operator=, "C++: HepMC3::BoolAttribute::operator=(const class HepMC3::BoolAttribute &) --> class HepMC3::BoolAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::VectorCharAttribute file:HepMC3/Attribute.h line:682
		pybind11::class_<HepMC3::VectorCharAttribute, std::shared_ptr<HepMC3::VectorCharAttribute>, PyCallBack_HepMC3_VectorCharAttribute, HepMC3::Attribute> cl(M("HepMC3"), "VectorCharAttribute", "Attribute that holds a vector of charegers of type  char\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::VectorCharAttribute(); }, [](){ return new PyCallBack_HepMC3_VectorCharAttribute(); } ) );
		cl.def( pybind11::init<class std::vector<char>>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_VectorCharAttribute const &o){ return new PyCallBack_HepMC3_VectorCharAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::VectorCharAttribute const &o){ return new HepMC3::VectorCharAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::VectorCharAttribute::*)(const std::string &)) &HepMC3::VectorCharAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::VectorCharAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::VectorCharAttribute::*)(std::string &) const) &HepMC3::VectorCharAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::VectorCharAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (class std::vector<char> (HepMC3::VectorCharAttribute::*)() const) &HepMC3::VectorCharAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::VectorCharAttribute::value() const --> class std::vector<char>");
		cl.def("set_value", (void (HepMC3::VectorCharAttribute::*)(const class std::vector<char> &)) &HepMC3::VectorCharAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::VectorCharAttribute::set_value(const class std::vector<char> &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::VectorCharAttribute & (HepMC3::VectorCharAttribute::*)(const class HepMC3::VectorCharAttribute &)) &HepMC3::VectorCharAttribute::operator=, "C++: HepMC3::VectorCharAttribute::operator=(const class HepMC3::VectorCharAttribute &) --> class HepMC3::VectorCharAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::VectorFloatAttribute file:HepMC3/Attribute.h line:727
		pybind11::class_<HepMC3::VectorFloatAttribute, std::shared_ptr<HepMC3::VectorFloatAttribute>, PyCallBack_HepMC3_VectorFloatAttribute, HepMC3::Attribute> cl(M("HepMC3"), "VectorFloatAttribute", "Attribute that holds a vector of floategers of type  float\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::VectorFloatAttribute(); }, [](){ return new PyCallBack_HepMC3_VectorFloatAttribute(); } ) );
		cl.def( pybind11::init<class std::vector<float>>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_VectorFloatAttribute const &o){ return new PyCallBack_HepMC3_VectorFloatAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::VectorFloatAttribute const &o){ return new HepMC3::VectorFloatAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::VectorFloatAttribute::*)(const std::string &)) &HepMC3::VectorFloatAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::VectorFloatAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::VectorFloatAttribute::*)(std::string &) const) &HepMC3::VectorFloatAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::VectorFloatAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (class std::vector<float> (HepMC3::VectorFloatAttribute::*)() const) &HepMC3::VectorFloatAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::VectorFloatAttribute::value() const --> class std::vector<float>");
		cl.def("set_value", (void (HepMC3::VectorFloatAttribute::*)(const class std::vector<float> &)) &HepMC3::VectorFloatAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::VectorFloatAttribute::set_value(const class std::vector<float> &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::VectorFloatAttribute & (HepMC3::VectorFloatAttribute::*)(const class HepMC3::VectorFloatAttribute &)) &HepMC3::VectorFloatAttribute::operator=, "C++: HepMC3::VectorFloatAttribute::operator=(const class HepMC3::VectorFloatAttribute &) --> class HepMC3::VectorFloatAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::VectorLongDoubleAttribute file:HepMC3/Attribute.h line:773
		pybind11::class_<HepMC3::VectorLongDoubleAttribute, std::shared_ptr<HepMC3::VectorLongDoubleAttribute>, PyCallBack_HepMC3_VectorLongDoubleAttribute, HepMC3::Attribute> cl(M("HepMC3"), "VectorLongDoubleAttribute", "Attribute that holds a vector of long doubleegers of type  long double\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::VectorLongDoubleAttribute(); }, [](){ return new PyCallBack_HepMC3_VectorLongDoubleAttribute(); } ) );
		cl.def( pybind11::init<class std::vector<long double>>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_VectorLongDoubleAttribute const &o){ return new PyCallBack_HepMC3_VectorLongDoubleAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::VectorLongDoubleAttribute const &o){ return new HepMC3::VectorLongDoubleAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::VectorLongDoubleAttribute::*)(const std::string &)) &HepMC3::VectorLongDoubleAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::VectorLongDoubleAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::VectorLongDoubleAttribute::*)(std::string &) const) &HepMC3::VectorLongDoubleAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::VectorLongDoubleAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (class std::vector<long double> (HepMC3::VectorLongDoubleAttribute::*)() const) &HepMC3::VectorLongDoubleAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::VectorLongDoubleAttribute::value() const --> class std::vector<long double>");
		cl.def("set_value", (void (HepMC3::VectorLongDoubleAttribute::*)(const class std::vector<long double> &)) &HepMC3::VectorLongDoubleAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::VectorLongDoubleAttribute::set_value(const class std::vector<long double> &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::VectorLongDoubleAttribute & (HepMC3::VectorLongDoubleAttribute::*)(const class HepMC3::VectorLongDoubleAttribute &)) &HepMC3::VectorLongDoubleAttribute::operator=, "C++: HepMC3::VectorLongDoubleAttribute::operator=(const class HepMC3::VectorLongDoubleAttribute &) --> class HepMC3::VectorLongDoubleAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::VectorLongLongAttribute file:HepMC3/Attribute.h line:820
		pybind11::class_<HepMC3::VectorLongLongAttribute, std::shared_ptr<HepMC3::VectorLongLongAttribute>, PyCallBack_HepMC3_VectorLongLongAttribute, HepMC3::Attribute> cl(M("HepMC3"), "VectorLongLongAttribute", "Attribute that holds a vector of long longegers of type  long long\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::VectorLongLongAttribute(); }, [](){ return new PyCallBack_HepMC3_VectorLongLongAttribute(); } ) );
		cl.def( pybind11::init<class std::vector<long long>>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_VectorLongLongAttribute const &o){ return new PyCallBack_HepMC3_VectorLongLongAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::VectorLongLongAttribute const &o){ return new HepMC3::VectorLongLongAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::VectorLongLongAttribute::*)(const std::string &)) &HepMC3::VectorLongLongAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::VectorLongLongAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::VectorLongLongAttribute::*)(std::string &) const) &HepMC3::VectorLongLongAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::VectorLongLongAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (class std::vector<long long> (HepMC3::VectorLongLongAttribute::*)() const) &HepMC3::VectorLongLongAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::VectorLongLongAttribute::value() const --> class std::vector<long long>");
		cl.def("set_value", (void (HepMC3::VectorLongLongAttribute::*)(const class std::vector<long long> &)) &HepMC3::VectorLongLongAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::VectorLongLongAttribute::set_value(const class std::vector<long long> &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::VectorLongLongAttribute & (HepMC3::VectorLongLongAttribute::*)(const class HepMC3::VectorLongLongAttribute &)) &HepMC3::VectorLongLongAttribute::operator=, "C++: HepMC3::VectorLongLongAttribute::operator=(const class HepMC3::VectorLongLongAttribute &) --> class HepMC3::VectorLongLongAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::VectorUIntAttribute file:HepMC3/Attribute.h line:865
		pybind11::class_<HepMC3::VectorUIntAttribute, std::shared_ptr<HepMC3::VectorUIntAttribute>, PyCallBack_HepMC3_VectorUIntAttribute, HepMC3::Attribute> cl(M("HepMC3"), "VectorUIntAttribute", "Attribute that holds a vector of unsigned integers of type  unsigned int\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::VectorUIntAttribute(); }, [](){ return new PyCallBack_HepMC3_VectorUIntAttribute(); } ) );
		cl.def( pybind11::init<class std::vector<unsigned int>>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_VectorUIntAttribute const &o){ return new PyCallBack_HepMC3_VectorUIntAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::VectorUIntAttribute const &o){ return new HepMC3::VectorUIntAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::VectorUIntAttribute::*)(const std::string &)) &HepMC3::VectorUIntAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::VectorUIntAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::VectorUIntAttribute::*)(std::string &) const) &HepMC3::VectorUIntAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::VectorUIntAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (class std::vector<unsigned int> (HepMC3::VectorUIntAttribute::*)() const) &HepMC3::VectorUIntAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::VectorUIntAttribute::value() const --> class std::vector<unsigned int>");
		cl.def("set_value", (void (HepMC3::VectorUIntAttribute::*)(const class std::vector<unsigned int> &)) &HepMC3::VectorUIntAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::VectorUIntAttribute::set_value(const class std::vector<unsigned int> &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::VectorUIntAttribute & (HepMC3::VectorUIntAttribute::*)(const class HepMC3::VectorUIntAttribute &)) &HepMC3::VectorUIntAttribute::operator=, "C++: HepMC3::VectorUIntAttribute::operator=(const class HepMC3::VectorUIntAttribute &) --> class HepMC3::VectorUIntAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::VectorULongAttribute file:HepMC3/Attribute.h line:910
		pybind11::class_<HepMC3::VectorULongAttribute, std::shared_ptr<HepMC3::VectorULongAttribute>, PyCallBack_HepMC3_VectorULongAttribute, HepMC3::Attribute> cl(M("HepMC3"), "VectorULongAttribute", "Attribute that holds a vector of unsigned longegers of type  unsigned long\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::VectorULongAttribute(); }, [](){ return new PyCallBack_HepMC3_VectorULongAttribute(); } ) );
		cl.def( pybind11::init<class std::vector<unsigned long>>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_VectorULongAttribute const &o){ return new PyCallBack_HepMC3_VectorULongAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::VectorULongAttribute const &o){ return new HepMC3::VectorULongAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::VectorULongAttribute::*)(const std::string &)) &HepMC3::VectorULongAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::VectorULongAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::VectorULongAttribute::*)(std::string &) const) &HepMC3::VectorULongAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::VectorULongAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (class std::vector<unsigned long> (HepMC3::VectorULongAttribute::*)() const) &HepMC3::VectorULongAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::VectorULongAttribute::value() const --> class std::vector<unsigned long>");
		cl.def("set_value", (void (HepMC3::VectorULongAttribute::*)(const class std::vector<unsigned long> &)) &HepMC3::VectorULongAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::VectorULongAttribute::set_value(const class std::vector<unsigned long> &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::VectorULongAttribute & (HepMC3::VectorULongAttribute::*)(const class HepMC3::VectorULongAttribute &)) &HepMC3::VectorULongAttribute::operator=, "C++: HepMC3::VectorULongAttribute::operator=(const class HepMC3::VectorULongAttribute &) --> class HepMC3::VectorULongAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
