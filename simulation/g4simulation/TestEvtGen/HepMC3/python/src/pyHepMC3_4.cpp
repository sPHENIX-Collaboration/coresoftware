#include <HepMC3/Attribute.h>
#include <HepMC3/Data/GenEventData.h>
#include <HepMC3/Data/GenRunInfoData.h>
#include <HepMC3/FourVector.h>
#include <HepMC3/GenCrossSection.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenHeavyIon.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/GenVertex.h>
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

// HepMC3::Attribute file:HepMC3/Attribute.h line:44
struct PyCallBack_HepMC3_Attribute : public HepMC3::Attribute {
	using HepMC3::Attribute::Attribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::Attribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Attribute::from_string\"");
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::Attribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::Attribute *>(this), "init");
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
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::Attribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Attribute::to_string\"");
	}
};

// HepMC3::IntAttribute file:HepMC3/Attribute.h line:157
struct PyCallBack_HepMC3_IntAttribute : public HepMC3::IntAttribute {
	using HepMC3::IntAttribute::IntAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::IntAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return IntAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::IntAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return IntAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::IntAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::IntAttribute *>(this), "init");
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

// HepMC3::LongAttribute file:HepMC3/Attribute.h line:198
struct PyCallBack_HepMC3_LongAttribute : public HepMC3::LongAttribute {
	using HepMC3::LongAttribute::LongAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::LongAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LongAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::LongAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LongAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::LongAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::LongAttribute *>(this), "init");
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

// HepMC3::DoubleAttribute file:HepMC3/Attribute.h line:241
struct PyCallBack_HepMC3_DoubleAttribute : public HepMC3::DoubleAttribute {
	using HepMC3::DoubleAttribute::DoubleAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::DoubleAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return DoubleAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::DoubleAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return DoubleAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::DoubleAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::DoubleAttribute *>(this), "init");
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

// HepMC3::FloatAttribute file:HepMC3/Attribute.h line:286
struct PyCallBack_HepMC3_FloatAttribute : public HepMC3::FloatAttribute {
	using HepMC3::FloatAttribute::FloatAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::FloatAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return FloatAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::FloatAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return FloatAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::FloatAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::FloatAttribute *>(this), "init");
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

// HepMC3::StringAttribute file:HepMC3/Attribute.h line:335
struct PyCallBack_HepMC3_StringAttribute : public HepMC3::StringAttribute {
	using HepMC3::StringAttribute::StringAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::StringAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return StringAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::StringAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return StringAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::StringAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::StringAttribute *>(this), "init");
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

// HepMC3::CharAttribute file:HepMC3/Attribute.h line:379
struct PyCallBack_HepMC3_CharAttribute : public HepMC3::CharAttribute {
	using HepMC3::CharAttribute::CharAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::CharAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return CharAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::CharAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return CharAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::CharAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::CharAttribute *>(this), "init");
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

// HepMC3::LongLongAttribute file:HepMC3/Attribute.h line:424
struct PyCallBack_HepMC3_LongLongAttribute : public HepMC3::LongLongAttribute {
	using HepMC3::LongLongAttribute::LongLongAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::LongLongAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LongLongAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::LongLongAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LongLongAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::LongLongAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::LongLongAttribute *>(this), "init");
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

// HepMC3::LongDoubleAttribute file:HepMC3/Attribute.h line:467
struct PyCallBack_HepMC3_LongDoubleAttribute : public HepMC3::LongDoubleAttribute {
	using HepMC3::LongDoubleAttribute::LongDoubleAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::LongDoubleAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LongDoubleAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::LongDoubleAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LongDoubleAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::LongDoubleAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::LongDoubleAttribute *>(this), "init");
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

// HepMC3::UIntAttribute file:HepMC3/Attribute.h line:514
struct PyCallBack_HepMC3_UIntAttribute : public HepMC3::UIntAttribute {
	using HepMC3::UIntAttribute::UIntAttribute;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::UIntAttribute *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UIntAttribute::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::UIntAttribute *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UIntAttribute::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::UIntAttribute *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::UIntAttribute *>(this), "init");
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

void bind_pyHepMC3_4(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // HepMC3::Attribute file:HepMC3/Attribute.h line:44
		pybind11::class_<HepMC3::Attribute, std::shared_ptr<HepMC3::Attribute>, PyCallBack_HepMC3_Attribute> cl(M("HepMC3"), "Attribute", "Base attribute class. ");
		cl.def( pybind11::init( [](){ return new PyCallBack_HepMC3_Attribute(); } ) );
		cl.def(pybind11::init<PyCallBack_HepMC3_Attribute const &>());
		cl.def("from_string", (bool (HepMC3::Attribute::*)(const std::string &)) &HepMC3::Attribute::from_string, "Fill class content from string.\n\nC++: HepMC3::Attribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("init", (bool (HepMC3::Attribute::*)()) &HepMC3::Attribute::init, "Optionally initialize the attribute after from_string.\n\nC++: HepMC3::Attribute::init() --> bool");
		cl.def("init", (bool (HepMC3::Attribute::*)(const class HepMC3::GenRunInfo &)) &HepMC3::Attribute::init, "Optionally initialize the attribute after from_string\n\n Is passed a reference to the GenRunInfo object to which the\n Attribute belongs.\n\nC++: HepMC3::Attribute::init(const class HepMC3::GenRunInfo &) --> bool", pybind11::arg(""));
		cl.def("to_string", (bool (HepMC3::Attribute::*)(std::string &) const) &HepMC3::Attribute::to_string, "Fill string from class content \n\nC++: HepMC3::Attribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("is_parsed", (bool (HepMC3::Attribute::*)() const) &HepMC3::Attribute::is_parsed, "Check if this attribute is parsed \n\nC++: HepMC3::Attribute::is_parsed() const --> bool");
		cl.def("unparsed_string", (const std::string & (HepMC3::Attribute::*)() const) &HepMC3::Attribute::unparsed_string, "Get unparsed string \n\nC++: HepMC3::Attribute::unparsed_string() const --> const std::string &", pybind11::return_value_policy::automatic);
		cl.def("event", (const class HepMC3::GenEvent * (HepMC3::Attribute::*)() const) &HepMC3::Attribute::event, "return the GenEvent to which this Attribute belongs, if at all. \n\nC++: HepMC3::Attribute::event() const --> const class HepMC3::GenEvent *", pybind11::return_value_policy::automatic);
		cl.def("particle", (class std::shared_ptr<class HepMC3::GenParticle> (HepMC3::Attribute::*)()) &HepMC3::Attribute::particle, "return the GenParticle to which this Attribute belongs, if at all. \n\nC++: HepMC3::Attribute::particle() --> class std::shared_ptr<class HepMC3::GenParticle>");
		cl.def("vertex", (class std::shared_ptr<class HepMC3::GenVertex> (HepMC3::Attribute::*)()) &HepMC3::Attribute::vertex, "return the GenVertex to which this Attribute belongs, if at all. \n\nC++: HepMC3::Attribute::vertex() --> class std::shared_ptr<class HepMC3::GenVertex>");
		cl.def("assign", (class HepMC3::Attribute & (HepMC3::Attribute::*)(const class HepMC3::Attribute &)) &HepMC3::Attribute::operator=, "C++: HepMC3::Attribute::operator=(const class HepMC3::Attribute &) --> class HepMC3::Attribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::IntAttribute file:HepMC3/Attribute.h line:157
		pybind11::class_<HepMC3::IntAttribute, std::shared_ptr<HepMC3::IntAttribute>, PyCallBack_HepMC3_IntAttribute, HepMC3::Attribute> cl(M("HepMC3"), "IntAttribute", "Attribute that holds an Integer implemented as an int\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::IntAttribute(); }, [](){ return new PyCallBack_HepMC3_IntAttribute(); } ) );
		cl.def( pybind11::init<int>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_IntAttribute const &o){ return new PyCallBack_HepMC3_IntAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::IntAttribute const &o){ return new HepMC3::IntAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::IntAttribute::*)(const std::string &)) &HepMC3::IntAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::IntAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::IntAttribute::*)(std::string &) const) &HepMC3::IntAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::IntAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (int (HepMC3::IntAttribute::*)() const) &HepMC3::IntAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::IntAttribute::value() const --> int");
		cl.def("set_value", (void (HepMC3::IntAttribute::*)(const int &)) &HepMC3::IntAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::IntAttribute::set_value(const int &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::IntAttribute & (HepMC3::IntAttribute::*)(const class HepMC3::IntAttribute &)) &HepMC3::IntAttribute::operator=, "C++: HepMC3::IntAttribute::operator=(const class HepMC3::IntAttribute &) --> class HepMC3::IntAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::LongAttribute file:HepMC3/Attribute.h line:198
		pybind11::class_<HepMC3::LongAttribute, std::shared_ptr<HepMC3::LongAttribute>, PyCallBack_HepMC3_LongAttribute, HepMC3::Attribute> cl(M("HepMC3"), "LongAttribute", "Attribute that holds an Integer implemented as an int\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::LongAttribute(); }, [](){ return new PyCallBack_HepMC3_LongAttribute(); } ) );
		cl.def( pybind11::init<long>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_LongAttribute const &o){ return new PyCallBack_HepMC3_LongAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::LongAttribute const &o){ return new HepMC3::LongAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::LongAttribute::*)(const std::string &)) &HepMC3::LongAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::LongAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::LongAttribute::*)(std::string &) const) &HepMC3::LongAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::LongAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (long (HepMC3::LongAttribute::*)() const) &HepMC3::LongAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::LongAttribute::value() const --> long");
		cl.def("set_value", (void (HepMC3::LongAttribute::*)(const long &)) &HepMC3::LongAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::LongAttribute::set_value(const long &) --> void", pybind11::arg("l"));
		cl.def("assign", (class HepMC3::LongAttribute & (HepMC3::LongAttribute::*)(const class HepMC3::LongAttribute &)) &HepMC3::LongAttribute::operator=, "C++: HepMC3::LongAttribute::operator=(const class HepMC3::LongAttribute &) --> class HepMC3::LongAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::DoubleAttribute file:HepMC3/Attribute.h line:241
		pybind11::class_<HepMC3::DoubleAttribute, std::shared_ptr<HepMC3::DoubleAttribute>, PyCallBack_HepMC3_DoubleAttribute, HepMC3::Attribute> cl(M("HepMC3"), "DoubleAttribute", "Attribute that holds a real number as a double.\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::DoubleAttribute(); }, [](){ return new PyCallBack_HepMC3_DoubleAttribute(); } ) );
		cl.def( pybind11::init<double>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_DoubleAttribute const &o){ return new PyCallBack_HepMC3_DoubleAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::DoubleAttribute const &o){ return new HepMC3::DoubleAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::DoubleAttribute::*)(const std::string &)) &HepMC3::DoubleAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::DoubleAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::DoubleAttribute::*)(std::string &) const) &HepMC3::DoubleAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::DoubleAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (double (HepMC3::DoubleAttribute::*)() const) &HepMC3::DoubleAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::DoubleAttribute::value() const --> double");
		cl.def("set_value", (void (HepMC3::DoubleAttribute::*)(const double &)) &HepMC3::DoubleAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::DoubleAttribute::set_value(const double &) --> void", pybind11::arg("d"));
		cl.def("assign", (class HepMC3::DoubleAttribute & (HepMC3::DoubleAttribute::*)(const class HepMC3::DoubleAttribute &)) &HepMC3::DoubleAttribute::operator=, "C++: HepMC3::DoubleAttribute::operator=(const class HepMC3::DoubleAttribute &) --> class HepMC3::DoubleAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::FloatAttribute file:HepMC3/Attribute.h line:286
		pybind11::class_<HepMC3::FloatAttribute, std::shared_ptr<HepMC3::FloatAttribute>, PyCallBack_HepMC3_FloatAttribute, HepMC3::Attribute> cl(M("HepMC3"), "FloatAttribute", "Attribute that holds a real number as a float.\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::FloatAttribute(); }, [](){ return new PyCallBack_HepMC3_FloatAttribute(); } ) );
		cl.def( pybind11::init<float>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_FloatAttribute const &o){ return new PyCallBack_HepMC3_FloatAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::FloatAttribute const &o){ return new HepMC3::FloatAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::FloatAttribute::*)(const std::string &)) &HepMC3::FloatAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::FloatAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::FloatAttribute::*)(std::string &) const) &HepMC3::FloatAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::FloatAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (float (HepMC3::FloatAttribute::*)() const) &HepMC3::FloatAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::FloatAttribute::value() const --> float");
		cl.def("set_value", (void (HepMC3::FloatAttribute::*)(const float &)) &HepMC3::FloatAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::FloatAttribute::set_value(const float &) --> void", pybind11::arg("f"));
		cl.def("assign", (class HepMC3::FloatAttribute & (HepMC3::FloatAttribute::*)(const class HepMC3::FloatAttribute &)) &HepMC3::FloatAttribute::operator=, "C++: HepMC3::FloatAttribute::operator=(const class HepMC3::FloatAttribute &) --> class HepMC3::FloatAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::StringAttribute file:HepMC3/Attribute.h line:335
		pybind11::class_<HepMC3::StringAttribute, std::shared_ptr<HepMC3::StringAttribute>, PyCallBack_HepMC3_StringAttribute, HepMC3::Attribute> cl(M("HepMC3"), "StringAttribute", "Attribute that holds a string\n\n  Default attribute constructed when reading input files.\n  It can be then parsed by other attributes or left as a string.\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::StringAttribute(); }, [](){ return new PyCallBack_HepMC3_StringAttribute(); } ) );
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("st") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_StringAttribute const &o){ return new PyCallBack_HepMC3_StringAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::StringAttribute const &o){ return new HepMC3::StringAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::StringAttribute::*)(const std::string &)) &HepMC3::StringAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::StringAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::StringAttribute::*)(std::string &) const) &HepMC3::StringAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::StringAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (std::string (HepMC3::StringAttribute::*)() const) &HepMC3::StringAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::StringAttribute::value() const --> std::string");
		cl.def("set_value", (void (HepMC3::StringAttribute::*)(const std::string &)) &HepMC3::StringAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::StringAttribute::set_value(const std::string &) --> void", pybind11::arg("s"));
		cl.def("assign", (class HepMC3::StringAttribute & (HepMC3::StringAttribute::*)(const class HepMC3::StringAttribute &)) &HepMC3::StringAttribute::operator=, "C++: HepMC3::StringAttribute::operator=(const class HepMC3::StringAttribute &) --> class HepMC3::StringAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::CharAttribute file:HepMC3/Attribute.h line:379
		pybind11::class_<HepMC3::CharAttribute, std::shared_ptr<HepMC3::CharAttribute>, PyCallBack_HepMC3_CharAttribute, HepMC3::Attribute> cl(M("HepMC3"), "CharAttribute", "Attribute that holds an Chareger implemented as an int\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::CharAttribute(); }, [](){ return new PyCallBack_HepMC3_CharAttribute(); } ) );
		cl.def( pybind11::init<char>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_CharAttribute const &o){ return new PyCallBack_HepMC3_CharAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::CharAttribute const &o){ return new HepMC3::CharAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::CharAttribute::*)(const std::string &)) &HepMC3::CharAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::CharAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::CharAttribute::*)(std::string &) const) &HepMC3::CharAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::CharAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (char (HepMC3::CharAttribute::*)() const) &HepMC3::CharAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::CharAttribute::value() const --> char");
		cl.def("set_value", (void (HepMC3::CharAttribute::*)(const char &)) &HepMC3::CharAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::CharAttribute::set_value(const char &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::CharAttribute & (HepMC3::CharAttribute::*)(const class HepMC3::CharAttribute &)) &HepMC3::CharAttribute::operator=, "C++: HepMC3::CharAttribute::operator=(const class HepMC3::CharAttribute &) --> class HepMC3::CharAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::LongLongAttribute file:HepMC3/Attribute.h line:424
		pybind11::class_<HepMC3::LongLongAttribute, std::shared_ptr<HepMC3::LongLongAttribute>, PyCallBack_HepMC3_LongLongAttribute, HepMC3::Attribute> cl(M("HepMC3"), "LongLongAttribute", "Attribute that holds an Integer implemented as an int\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::LongLongAttribute(); }, [](){ return new PyCallBack_HepMC3_LongLongAttribute(); } ) );
		cl.def( pybind11::init<long long>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_LongLongAttribute const &o){ return new PyCallBack_HepMC3_LongLongAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::LongLongAttribute const &o){ return new HepMC3::LongLongAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::LongLongAttribute::*)(const std::string &)) &HepMC3::LongLongAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::LongLongAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::LongLongAttribute::*)(std::string &) const) &HepMC3::LongLongAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::LongLongAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (long long (HepMC3::LongLongAttribute::*)() const) &HepMC3::LongLongAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::LongLongAttribute::value() const --> long long");
		cl.def("set_value", (void (HepMC3::LongLongAttribute::*)(const long long &)) &HepMC3::LongLongAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::LongLongAttribute::set_value(const long long &) --> void", pybind11::arg("l"));
		cl.def("assign", (class HepMC3::LongLongAttribute & (HepMC3::LongLongAttribute::*)(const class HepMC3::LongLongAttribute &)) &HepMC3::LongLongAttribute::operator=, "C++: HepMC3::LongLongAttribute::operator=(const class HepMC3::LongLongAttribute &) --> class HepMC3::LongLongAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::LongDoubleAttribute file:HepMC3/Attribute.h line:467
		pybind11::class_<HepMC3::LongDoubleAttribute, std::shared_ptr<HepMC3::LongDoubleAttribute>, PyCallBack_HepMC3_LongDoubleAttribute, HepMC3::Attribute> cl(M("HepMC3"), "LongDoubleAttribute", "Attribute that holds a real number as a double.\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::LongDoubleAttribute(); }, [](){ return new PyCallBack_HepMC3_LongDoubleAttribute(); } ) );
		cl.def( pybind11::init<long double>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_LongDoubleAttribute const &o){ return new PyCallBack_HepMC3_LongDoubleAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::LongDoubleAttribute const &o){ return new HepMC3::LongDoubleAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::LongDoubleAttribute::*)(const std::string &)) &HepMC3::LongDoubleAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::LongDoubleAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::LongDoubleAttribute::*)(std::string &) const) &HepMC3::LongDoubleAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::LongDoubleAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (long double (HepMC3::LongDoubleAttribute::*)() const) &HepMC3::LongDoubleAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::LongDoubleAttribute::value() const --> long double");
		cl.def("set_value", (void (HepMC3::LongDoubleAttribute::*)(const long double &)) &HepMC3::LongDoubleAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::LongDoubleAttribute::set_value(const long double &) --> void", pybind11::arg("d"));
		cl.def("assign", (class HepMC3::LongDoubleAttribute & (HepMC3::LongDoubleAttribute::*)(const class HepMC3::LongDoubleAttribute &)) &HepMC3::LongDoubleAttribute::operator=, "C++: HepMC3::LongDoubleAttribute::operator=(const class HepMC3::LongDoubleAttribute &) --> class HepMC3::LongDoubleAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::UIntAttribute file:HepMC3/Attribute.h line:514
		pybind11::class_<HepMC3::UIntAttribute, std::shared_ptr<HepMC3::UIntAttribute>, PyCallBack_HepMC3_UIntAttribute, HepMC3::Attribute> cl(M("HepMC3"), "UIntAttribute", "Attribute that holds an unsigned int\n\n  \n\n ");
		cl.def( pybind11::init( [](){ return new HepMC3::UIntAttribute(); }, [](){ return new PyCallBack_HepMC3_UIntAttribute(); } ) );
		cl.def( pybind11::init<unsigned int>(), pybind11::arg("val") );

		cl.def( pybind11::init( [](PyCallBack_HepMC3_UIntAttribute const &o){ return new PyCallBack_HepMC3_UIntAttribute(o); } ) );
		cl.def( pybind11::init( [](HepMC3::UIntAttribute const &o){ return new HepMC3::UIntAttribute(o); } ) );
		cl.def("from_string", (bool (HepMC3::UIntAttribute::*)(const std::string &)) &HepMC3::UIntAttribute::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::UIntAttribute::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::UIntAttribute::*)(std::string &) const) &HepMC3::UIntAttribute::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::UIntAttribute::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("value", (unsigned int (HepMC3::UIntAttribute::*)() const) &HepMC3::UIntAttribute::value, "get the value associated to this Attribute. \n\nC++: HepMC3::UIntAttribute::value() const --> unsigned int");
		cl.def("set_value", (void (HepMC3::UIntAttribute::*)(const unsigned int &)) &HepMC3::UIntAttribute::set_value, "set the value associated to this Attribute. \n\nC++: HepMC3::UIntAttribute::set_value(const unsigned int &) --> void", pybind11::arg("i"));
		cl.def("assign", (class HepMC3::UIntAttribute & (HepMC3::UIntAttribute::*)(const class HepMC3::UIntAttribute &)) &HepMC3::UIntAttribute::operator=, "C++: HepMC3::UIntAttribute::operator=(const class HepMC3::UIntAttribute &) --> class HepMC3::UIntAttribute &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
