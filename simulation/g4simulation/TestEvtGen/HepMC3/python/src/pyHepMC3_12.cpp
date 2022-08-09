#include <HepMC3/Attribute.h>
#include <HepMC3/Data/GenEventData.h>
#include <HepMC3/Data/GenParticleData.h>
#include <HepMC3/Data/GenRunInfoData.h>
#include <HepMC3/FourVector.h>
#include <HepMC3/GenCrossSection.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenHeavyIon.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/Reader.h>
#include <HepMC3/ReaderAscii.h>
#include <HepMC3/ReaderAsciiHepMC2.h>
#include <HepMC3/Writer.h>
#include <HepMC3/WriterAsciiHepMC2.h>
#include <functional>
#include <ios>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
#include <sstream> // __str__
#include <streambuf>
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

// HepMC3::ReaderAscii file:HepMC3/ReaderAscii.h line:29
struct PyCallBack_HepMC3_ReaderAscii : public HepMC3::ReaderAscii {
	using HepMC3::ReaderAscii::ReaderAscii;

	bool skip(const int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAscii *>(this), "skip");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderAscii::skip(a0);
	}
	bool read_event(class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAscii *>(this), "read_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderAscii::read_event(a0);
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAscii *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderAscii::failed();
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAscii *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReaderAscii::close();
	}
	class std::shared_ptr<class HepMC3::GenRunInfo> run_info() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAscii *>(this), "run_info");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class HepMC3::GenRunInfo>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class HepMC3::GenRunInfo>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class HepMC3::GenRunInfo>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class HepMC3::GenRunInfo>>(std::move(o));
		}
		return Reader::run_info();
	}
	void set_options(const class std::map<std::string, std::string > & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAscii *>(this), "set_options");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Reader::set_options(a0);
	}
	using _binder_ret_0 = std::map<std::string, std::string >;
	_binder_ret_0 get_options() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAscii *>(this), "get_options");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return Reader::get_options();
	}
	void set_run_info(class std::shared_ptr<class HepMC3::GenRunInfo> a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAscii *>(this), "set_run_info");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Reader::set_run_info(a0);
	}
};

// HepMC3::WriterAsciiHepMC2 file:HepMC3/WriterAsciiHepMC2.h line:26
struct PyCallBack_HepMC3_WriterAsciiHepMC2 : public HepMC3::WriterAsciiHepMC2 {
	using HepMC3::WriterAsciiHepMC2::WriterAsciiHepMC2;

	void write_event(const class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterAsciiHepMC2 *>(this), "write_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterAsciiHepMC2::write_event(a0);
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterAsciiHepMC2 *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return WriterAsciiHepMC2::failed();
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterAsciiHepMC2 *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterAsciiHepMC2::close();
	}
	void set_run_info(class std::shared_ptr<class HepMC3::GenRunInfo> a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterAsciiHepMC2 *>(this), "set_run_info");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Writer::set_run_info(a0);
	}
	class std::shared_ptr<class HepMC3::GenRunInfo> run_info() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterAsciiHepMC2 *>(this), "run_info");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class HepMC3::GenRunInfo>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class HepMC3::GenRunInfo>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class HepMC3::GenRunInfo>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class HepMC3::GenRunInfo>>(std::move(o));
		}
		return Writer::run_info();
	}
	void set_options(const class std::map<std::string, std::string > & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterAsciiHepMC2 *>(this), "set_options");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Writer::set_options(a0);
	}
	using _binder_ret_0 = std::map<std::string, std::string >;
	_binder_ret_0 get_options() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterAsciiHepMC2 *>(this), "get_options");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return Writer::get_options();
	}
};

// HepMC3::ReaderAsciiHepMC2 file:HepMC3/ReaderAsciiHepMC2.h line:30
struct PyCallBack_HepMC3_ReaderAsciiHepMC2 : public HepMC3::ReaderAsciiHepMC2 {
	using HepMC3::ReaderAsciiHepMC2::ReaderAsciiHepMC2;

	bool skip(const int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAsciiHepMC2 *>(this), "skip");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderAsciiHepMC2::skip(a0);
	}
	bool read_event(class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAsciiHepMC2 *>(this), "read_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderAsciiHepMC2::read_event(a0);
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAsciiHepMC2 *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderAsciiHepMC2::failed();
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAsciiHepMC2 *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReaderAsciiHepMC2::close();
	}
	class std::shared_ptr<class HepMC3::GenRunInfo> run_info() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAsciiHepMC2 *>(this), "run_info");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class HepMC3::GenRunInfo>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class HepMC3::GenRunInfo>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class HepMC3::GenRunInfo>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class HepMC3::GenRunInfo>>(std::move(o));
		}
		return Reader::run_info();
	}
	void set_options(const class std::map<std::string, std::string > & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAsciiHepMC2 *>(this), "set_options");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Reader::set_options(a0);
	}
	using _binder_ret_0 = std::map<std::string, std::string >;
	_binder_ret_0 get_options() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAsciiHepMC2 *>(this), "get_options");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return Reader::get_options();
	}
	void set_run_info(class std::shared_ptr<class HepMC3::GenRunInfo> a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAsciiHepMC2 *>(this), "set_run_info");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Reader::set_run_info(a0);
	}
};

void bind_pyHepMC3_12(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // HepMC3::ReaderAscii file:HepMC3/ReaderAscii.h line:29
		pybind11::class_<HepMC3::ReaderAscii, std::shared_ptr<HepMC3::ReaderAscii>, PyCallBack_HepMC3_ReaderAscii, HepMC3::Reader> cl(M("HepMC3"), "ReaderAscii", "");
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("filename") );

		cl.def("skip", (bool (HepMC3::ReaderAscii::*)(const int)) &HepMC3::ReaderAscii::skip, "skip events\n\nC++: HepMC3::ReaderAscii::skip(const int) --> bool", pybind11::arg(""));
		cl.def("read_event", (bool (HepMC3::ReaderAscii::*)(class HepMC3::GenEvent &)) &HepMC3::ReaderAscii::read_event, "Load event from file\n\n \n Event to be filled\n\nC++: HepMC3::ReaderAscii::read_event(class HepMC3::GenEvent &) --> bool", pybind11::arg("evt"));
		cl.def("failed", (bool (HepMC3::ReaderAscii::*)()) &HepMC3::ReaderAscii::failed, "Return status of the stream\n\nC++: HepMC3::ReaderAscii::failed() --> bool");
		cl.def("close", (void (HepMC3::ReaderAscii::*)()) &HepMC3::ReaderAscii::close, "Close file stream\n\nC++: HepMC3::ReaderAscii::close() --> void");
	}
	{ // HepMC3::WriterAsciiHepMC2 file:HepMC3/WriterAsciiHepMC2.h line:26
		pybind11::class_<HepMC3::WriterAsciiHepMC2, std::shared_ptr<HepMC3::WriterAsciiHepMC2>, PyCallBack_HepMC3_WriterAsciiHepMC2, HepMC3::Writer> cl(M("HepMC3"), "WriterAsciiHepMC2", "");
		cl.def( pybind11::init( [](const std::string & a0){ return new HepMC3::WriterAsciiHepMC2(a0); }, [](const std::string & a0){ return new PyCallBack_HepMC3_WriterAsciiHepMC2(a0); } ), "doc");
		cl.def( pybind11::init<const std::string &, class std::shared_ptr<class HepMC3::GenRunInfo>>(), pybind11::arg("filename"), pybind11::arg("run") );

		cl.def("write_event", (void (HepMC3::WriterAsciiHepMC2::*)(const class HepMC3::GenEvent &)) &HepMC3::WriterAsciiHepMC2::write_event, "Write event to file\n\n \n Event to be serialized\n\nC++: HepMC3::WriterAsciiHepMC2::write_event(const class HepMC3::GenEvent &) --> void", pybind11::arg("evt"));
		cl.def("write_run_info", (void (HepMC3::WriterAsciiHepMC2::*)()) &HepMC3::WriterAsciiHepMC2::write_run_info, "Write the GenRunInfo object to file.\n\nC++: HepMC3::WriterAsciiHepMC2::write_run_info() --> void");
		cl.def("failed", (bool (HepMC3::WriterAsciiHepMC2::*)()) &HepMC3::WriterAsciiHepMC2::failed, "Return status of the stream\n\nC++: HepMC3::WriterAsciiHepMC2::failed() --> bool");
		cl.def("close", (void (HepMC3::WriterAsciiHepMC2::*)()) &HepMC3::WriterAsciiHepMC2::close, "Close file stream\n\nC++: HepMC3::WriterAsciiHepMC2::close() --> void");
		cl.def("set_precision", (void (HepMC3::WriterAsciiHepMC2::*)(const int &)) &HepMC3::WriterAsciiHepMC2::set_precision, "Set output precision\n\n Available range is [2,24]. Default is 16.\n\nC++: HepMC3::WriterAsciiHepMC2::set_precision(const int &) --> void", pybind11::arg("prec"));
		cl.def("precision", (int (HepMC3::WriterAsciiHepMC2::*)() const) &HepMC3::WriterAsciiHepMC2::precision, "Return output precision\n\nC++: HepMC3::WriterAsciiHepMC2::precision() const --> int");
	}
	{ // HepMC3::ReaderAsciiHepMC2 file:HepMC3/ReaderAsciiHepMC2.h line:30
		pybind11::class_<HepMC3::ReaderAsciiHepMC2, std::shared_ptr<HepMC3::ReaderAsciiHepMC2>, PyCallBack_HepMC3_ReaderAsciiHepMC2, HepMC3::Reader> cl(M("HepMC3"), "ReaderAsciiHepMC2", "");
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("filename") );

		cl.def("skip", (bool (HepMC3::ReaderAsciiHepMC2::*)(const int)) &HepMC3::ReaderAsciiHepMC2::skip, "skip events\n\nC++: HepMC3::ReaderAsciiHepMC2::skip(const int) --> bool", pybind11::arg(""));
		cl.def("read_event", (bool (HepMC3::ReaderAsciiHepMC2::*)(class HepMC3::GenEvent &)) &HepMC3::ReaderAsciiHepMC2::read_event, "Implementation of Reader::read_event \n\nC++: HepMC3::ReaderAsciiHepMC2::read_event(class HepMC3::GenEvent &) --> bool", pybind11::arg("evt"));
		cl.def("failed", (bool (HepMC3::ReaderAsciiHepMC2::*)()) &HepMC3::ReaderAsciiHepMC2::failed, "Return status of the stream\n\nC++: HepMC3::ReaderAsciiHepMC2::failed() --> bool");
		cl.def("close", (void (HepMC3::ReaderAsciiHepMC2::*)()) &HepMC3::ReaderAsciiHepMC2::close, "Close file stream\n\nC++: HepMC3::ReaderAsciiHepMC2::close() --> void");
	}
}
