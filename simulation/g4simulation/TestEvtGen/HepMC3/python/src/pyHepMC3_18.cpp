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
#include <HepMC3/HEPEVT_Wrapper_Runtime.h>
#include <HepMC3/Reader.h>
#include <HepMC3/ReaderLHEF.h>
#include <HepMC3/ReaderPlugin.h>
#include <HepMC3/WriterPlugin.h>
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

// HepMC3::ReaderLHEF file:HepMC3/ReaderLHEF.h line:34
struct PyCallBack_HepMC3_ReaderLHEF : public HepMC3::ReaderLHEF {
	using HepMC3::ReaderLHEF::ReaderLHEF;

	bool skip(const int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderLHEF *>(this), "skip");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderLHEF::skip(a0);
	}
	bool read_event(class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderLHEF *>(this), "read_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderLHEF::read_event(a0);
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderLHEF *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReaderLHEF::close();
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderLHEF *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderLHEF::failed();
	}
	class std::shared_ptr<class HepMC3::GenRunInfo> run_info() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderLHEF *>(this), "run_info");
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
	void set_options(const std::map<std::string, std::string > & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderLHEF *>(this), "set_options");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderLHEF *>(this), "get_options");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderLHEF *>(this), "set_run_info");
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

// HepMC3::ReaderPlugin file:HepMC3/ReaderPlugin.h line:23
struct PyCallBack_HepMC3_ReaderPlugin : public HepMC3::ReaderPlugin {
	using HepMC3::ReaderPlugin::ReaderPlugin;

	bool read_event(class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderPlugin *>(this), "read_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderPlugin::read_event(a0);
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderPlugin *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReaderPlugin::close();
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderPlugin *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderPlugin::failed();
	}
	class std::shared_ptr<class HepMC3::GenRunInfo> run_info() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderPlugin *>(this), "run_info");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class HepMC3::GenRunInfo>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class HepMC3::GenRunInfo>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class HepMC3::GenRunInfo>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class HepMC3::GenRunInfo>>(std::move(o));
		}
		return ReaderPlugin::run_info();
	}
	void set_options(const class std::map<std::string, std::string > & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderPlugin *>(this), "set_options");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReaderPlugin::set_options(a0);
	}
	using _binder_ret_0 = std::map<std::string, std::string >;
	_binder_ret_0 get_options() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderPlugin *>(this), "get_options");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return ReaderPlugin::get_options();
	}
	void set_run_info(class std::shared_ptr<class HepMC3::GenRunInfo> a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderPlugin *>(this), "set_run_info");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReaderPlugin::set_run_info(a0);
	}
	bool skip(const int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderPlugin *>(this), "skip");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return Reader::skip(a0);
	}
};

// HepMC3::WriterPlugin file:HepMC3/WriterPlugin.h line:23
struct PyCallBack_HepMC3_WriterPlugin : public HepMC3::WriterPlugin {
	using HepMC3::WriterPlugin::WriterPlugin;

	void write_event(const class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterPlugin *>(this), "write_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterPlugin::write_event(a0);
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterPlugin *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterPlugin::close();
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterPlugin *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return WriterPlugin::failed();
	}
	class std::shared_ptr<class HepMC3::GenRunInfo> run_info() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterPlugin *>(this), "run_info");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class HepMC3::GenRunInfo>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class HepMC3::GenRunInfo>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class HepMC3::GenRunInfo>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class HepMC3::GenRunInfo>>(std::move(o));
		}
		return WriterPlugin::run_info();
	}
	void set_options(const class std::map<std::string, std::string > & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterPlugin *>(this), "set_options");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterPlugin::set_options(a0);
	}
	using _binder_ret_0 = std::map<std::string, std::string >;
	_binder_ret_0 get_options() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterPlugin *>(this), "get_options");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return WriterPlugin::get_options();
	}
	void set_run_info(class std::shared_ptr<class HepMC3::GenRunInfo> a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterPlugin *>(this), "set_run_info");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterPlugin::set_run_info(a0);
	}
};

void bind_pyHepMC3_18(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // HepMC3::ReaderLHEF file:HepMC3/ReaderLHEF.h line:34
		pybind11::class_<HepMC3::ReaderLHEF, std::shared_ptr<HepMC3::ReaderLHEF>, PyCallBack_HepMC3_ReaderLHEF, HepMC3::Reader> cl(M("HepMC3"), "ReaderLHEF", "");
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("filename") );

		cl.def("skip", (bool (HepMC3::ReaderLHEF::*)(const int)) &HepMC3::ReaderLHEF::skip, "skip events\n\nC++: HepMC3::ReaderLHEF::skip(const int) --> bool", pybind11::arg(""));
		cl.def("read_event", (bool (HepMC3::ReaderLHEF::*)(class HepMC3::GenEvent &)) &HepMC3::ReaderLHEF::read_event, "Reading event \n\nC++: HepMC3::ReaderLHEF::read_event(class HepMC3::GenEvent &) --> bool", pybind11::arg("ev"));
		cl.def("close", (void (HepMC3::ReaderLHEF::*)()) &HepMC3::ReaderLHEF::close, "Close \n\nC++: HepMC3::ReaderLHEF::close() --> void");
		cl.def("failed", (bool (HepMC3::ReaderLHEF::*)()) &HepMC3::ReaderLHEF::failed, "State \n\nC++: HepMC3::ReaderLHEF::failed() --> bool");
	}
	{ // HepMC3::ReaderPlugin file:HepMC3/ReaderPlugin.h line:23
		pybind11::class_<HepMC3::ReaderPlugin, std::shared_ptr<HepMC3::ReaderPlugin>, PyCallBack_HepMC3_ReaderPlugin, HepMC3::Reader> cl(M("HepMC3"), "ReaderPlugin", "");
		cl.def( pybind11::init<const std::string &, const std::string &, const std::string &>(), pybind11::arg("filename"), pybind11::arg("libname"), pybind11::arg("newreader") );

		cl.def("read_event", (bool (HepMC3::ReaderPlugin::*)(class HepMC3::GenEvent &)) &HepMC3::ReaderPlugin::read_event, "Reading event \n\nC++: HepMC3::ReaderPlugin::read_event(class HepMC3::GenEvent &) --> bool", pybind11::arg("ev"));
		cl.def("close", (void (HepMC3::ReaderPlugin::*)()) &HepMC3::ReaderPlugin::close, "Close \n\nC++: HepMC3::ReaderPlugin::close() --> void");
		cl.def("failed", (bool (HepMC3::ReaderPlugin::*)()) &HepMC3::ReaderPlugin::failed, "State \n\nC++: HepMC3::ReaderPlugin::failed() --> bool");
		cl.def("run_info", (class std::shared_ptr<class HepMC3::GenRunInfo> (HepMC3::ReaderPlugin::*)() const) &HepMC3::ReaderPlugin::run_info, "Get the global GenRunInfo object. \n\nC++: HepMC3::ReaderPlugin::run_info() const --> class std::shared_ptr<class HepMC3::GenRunInfo>");
		cl.def("set_options", (void (HepMC3::ReaderPlugin::*)(const class std::map<std::string, std::string > &)) &HepMC3::ReaderPlugin::set_options, "Set options \n\nC++: HepMC3::ReaderPlugin::set_options(const class std::map<std::string, std::string > &) --> void", pybind11::arg("options"));
		cl.def("get_options", (class std::map<std::string, std::string > (HepMC3::ReaderPlugin::*)() const) &HepMC3::ReaderPlugin::get_options, "Get options  \n\nC++: HepMC3::ReaderPlugin::get_options() const --> class std::map<std::string, std::string >");
		cl.def("set_run_info", (void (HepMC3::ReaderPlugin::*)(class std::shared_ptr<class HepMC3::GenRunInfo>)) &HepMC3::ReaderPlugin::set_run_info, "Set the global GenRunInfo object.\n\nC++: HepMC3::ReaderPlugin::set_run_info(class std::shared_ptr<class HepMC3::GenRunInfo>) --> void", pybind11::arg("run"));
	}
	{ // HepMC3::WriterPlugin file:HepMC3/WriterPlugin.h line:23
		pybind11::class_<HepMC3::WriterPlugin, std::shared_ptr<HepMC3::WriterPlugin>, PyCallBack_HepMC3_WriterPlugin, HepMC3::Writer> cl(M("HepMC3"), "WriterPlugin", "");
		cl.def( pybind11::init( [](const std::string & a0, const std::string & a1, const std::string & a2){ return new HepMC3::WriterPlugin(a0, a1, a2); }, [](const std::string & a0, const std::string & a1, const std::string & a2){ return new PyCallBack_HepMC3_WriterPlugin(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<const std::string &, const std::string &, const std::string &, class std::shared_ptr<class HepMC3::GenRunInfo>>(), pybind11::arg("filename"), pybind11::arg("libname"), pybind11::arg("newwriter"), pybind11::arg("run") );

		cl.def("write_event", (void (HepMC3::WriterPlugin::*)(const class HepMC3::GenEvent &)) &HepMC3::WriterPlugin::write_event, "Reading event \n\nC++: HepMC3::WriterPlugin::write_event(const class HepMC3::GenEvent &) --> void", pybind11::arg("ev"));
		cl.def("close", (void (HepMC3::WriterPlugin::*)()) &HepMC3::WriterPlugin::close, "Close \n\nC++: HepMC3::WriterPlugin::close() --> void");
		cl.def("failed", (bool (HepMC3::WriterPlugin::*)()) &HepMC3::WriterPlugin::failed, "State \n\nC++: HepMC3::WriterPlugin::failed() --> bool");
		cl.def("run_info", (class std::shared_ptr<class HepMC3::GenRunInfo> (HepMC3::WriterPlugin::*)() const) &HepMC3::WriterPlugin::run_info, "Get the global GenRunInfo object. \n\nC++: HepMC3::WriterPlugin::run_info() const --> class std::shared_ptr<class HepMC3::GenRunInfo>");
		cl.def("set_options", (void (HepMC3::WriterPlugin::*)(const class std::map<std::string, std::string > &)) &HepMC3::WriterPlugin::set_options, "Set options \n\nC++: HepMC3::WriterPlugin::set_options(const class std::map<std::string, std::string > &) --> void", pybind11::arg("options"));
		cl.def("get_options", (class std::map<std::string, std::string > (HepMC3::WriterPlugin::*)() const) &HepMC3::WriterPlugin::get_options, "Get options  \n\nC++: HepMC3::WriterPlugin::get_options() const --> class std::map<std::string, std::string >");
		cl.def("set_run_info", (void (HepMC3::WriterPlugin::*)(class std::shared_ptr<class HepMC3::GenRunInfo>)) &HepMC3::WriterPlugin::set_run_info, "Set the global GenRunInfo object.\n\nC++: HepMC3::WriterPlugin::set_run_info(class std::shared_ptr<class HepMC3::GenRunInfo>) --> void", pybind11::arg("run"));
	}
	{ // HepMC3::HEPEVT_Wrapper_Runtime file:HepMC3/HEPEVT_Wrapper_Runtime.h line:29
		pybind11::class_<HepMC3::HEPEVT_Wrapper_Runtime, std::shared_ptr<HepMC3::HEPEVT_Wrapper_Runtime>> cl(M("HepMC3"), "HEPEVT_Wrapper_Runtime", "");
		cl.def( pybind11::init( [](){ return new HepMC3::HEPEVT_Wrapper_Runtime(); } ) );
		cl.def( pybind11::init( [](HepMC3::HEPEVT_Wrapper_Runtime const &o){ return new HepMC3::HEPEVT_Wrapper_Runtime(o); } ) );
		cl.def("zero_everything", (void (HepMC3::HEPEVT_Wrapper_Runtime::*)()) &HepMC3::HEPEVT_Wrapper_Runtime::zero_everything, "Set all entries in HEPEVT to zero \n\nC++: HepMC3::HEPEVT_Wrapper_Runtime::zero_everything() --> void");
		cl.def("GenEvent_to_HEPEVT", (bool (HepMC3::HEPEVT_Wrapper_Runtime::*)(const class HepMC3::GenEvent *)) &HepMC3::HEPEVT_Wrapper_Runtime::GenEvent_to_HEPEVT, "Convert GenEvent to HEPEVT\n\nC++: HepMC3::HEPEVT_Wrapper_Runtime::GenEvent_to_HEPEVT(const class HepMC3::GenEvent *) --> bool", pybind11::arg("evt"));
		cl.def("HEPEVT_to_GenEvent", (bool (HepMC3::HEPEVT_Wrapper_Runtime::*)(class HepMC3::GenEvent *) const) &HepMC3::HEPEVT_Wrapper_Runtime::HEPEVT_to_GenEvent, "Convert HEPEVT to GenEvent\n\nC++: HepMC3::HEPEVT_Wrapper_Runtime::HEPEVT_to_GenEvent(class HepMC3::GenEvent *) const --> bool", pybind11::arg("evt"));
		cl.def("fix_daughters", (bool (HepMC3::HEPEVT_Wrapper_Runtime::*)()) &HepMC3::HEPEVT_Wrapper_Runtime::fix_daughters, "Tries to fix list of daughters \n\nC++: HepMC3::HEPEVT_Wrapper_Runtime::fix_daughters() --> bool");
		cl.def("allocate_internal_storage", (void (HepMC3::HEPEVT_Wrapper_Runtime::*)()) &HepMC3::HEPEVT_Wrapper_Runtime::allocate_internal_storage, "C++: HepMC3::HEPEVT_Wrapper_Runtime::allocate_internal_storage() --> void");
		cl.def("copy_to_internal_storage", (void (HepMC3::HEPEVT_Wrapper_Runtime::*)(char *, int)) &HepMC3::HEPEVT_Wrapper_Runtime::copy_to_internal_storage, "C++: HepMC3::HEPEVT_Wrapper_Runtime::copy_to_internal_storage(char *, int) --> void", pybind11::arg("c"), pybind11::arg("N"));
		cl.def("set_max_number_entries", (void (HepMC3::HEPEVT_Wrapper_Runtime::*)(unsigned int)) &HepMC3::HEPEVT_Wrapper_Runtime::set_max_number_entries, "C++: HepMC3::HEPEVT_Wrapper_Runtime::set_max_number_entries(unsigned int) --> void", pybind11::arg("size"));
		cl.def("set_hepevt_address", (void (HepMC3::HEPEVT_Wrapper_Runtime::*)(char *)) &HepMC3::HEPEVT_Wrapper_Runtime::set_hepevt_address, "C++: HepMC3::HEPEVT_Wrapper_Runtime::set_hepevt_address(char *) --> void", pybind11::arg("c"));
		cl.def("max_number_entries", (int (HepMC3::HEPEVT_Wrapper_Runtime::*)() const) &HepMC3::HEPEVT_Wrapper_Runtime::max_number_entries, "C++: HepMC3::HEPEVT_Wrapper_Runtime::max_number_entries() const --> int");
		cl.def("event_number", (int (HepMC3::HEPEVT_Wrapper_Runtime::*)() const) &HepMC3::HEPEVT_Wrapper_Runtime::event_number, "C++: HepMC3::HEPEVT_Wrapper_Runtime::event_number() const --> int");
		cl.def("number_entries", (int (HepMC3::HEPEVT_Wrapper_Runtime::*)() const) &HepMC3::HEPEVT_Wrapper_Runtime::number_entries, "C++: HepMC3::HEPEVT_Wrapper_Runtime::number_entries() const --> int");
		cl.def("status", (int (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::status, "C++: HepMC3::HEPEVT_Wrapper_Runtime::status(const int) const --> int", pybind11::arg("index"));
		cl.def("id", (int (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::id, "C++: HepMC3::HEPEVT_Wrapper_Runtime::id(const int) const --> int", pybind11::arg("index"));
		cl.def("first_parent", (int (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::first_parent, "C++: HepMC3::HEPEVT_Wrapper_Runtime::first_parent(const int) const --> int", pybind11::arg("index"));
		cl.def("last_parent", (int (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::last_parent, "C++: HepMC3::HEPEVT_Wrapper_Runtime::last_parent(const int) const --> int", pybind11::arg("index"));
		cl.def("first_child", (int (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::first_child, "C++: HepMC3::HEPEVT_Wrapper_Runtime::first_child(const int) const --> int", pybind11::arg("index"));
		cl.def("last_child", (int (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::last_child, "C++: HepMC3::HEPEVT_Wrapper_Runtime::last_child(const int) const --> int", pybind11::arg("index"));
		cl.def("px", (double (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::px, "C++: HepMC3::HEPEVT_Wrapper_Runtime::px(const int) const --> double", pybind11::arg("index"));
		cl.def("py", (double (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::py, "C++: HepMC3::HEPEVT_Wrapper_Runtime::py(const int) const --> double", pybind11::arg("index"));
		cl.def("pz", (double (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::pz, "C++: HepMC3::HEPEVT_Wrapper_Runtime::pz(const int) const --> double", pybind11::arg("index"));
		cl.def("e", (double (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::e, "C++: HepMC3::HEPEVT_Wrapper_Runtime::e(const int) const --> double", pybind11::arg("index"));
		cl.def("m", (double (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::m, "C++: HepMC3::HEPEVT_Wrapper_Runtime::m(const int) const --> double", pybind11::arg("index"));
		cl.def("x", (double (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::x, "C++: HepMC3::HEPEVT_Wrapper_Runtime::x(const int) const --> double", pybind11::arg("index"));
		cl.def("y", (double (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::y, "C++: HepMC3::HEPEVT_Wrapper_Runtime::y(const int) const --> double", pybind11::arg("index"));
		cl.def("z", (double (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::z, "C++: HepMC3::HEPEVT_Wrapper_Runtime::z(const int) const --> double", pybind11::arg("index"));
		cl.def("t", (double (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::t, "C++: HepMC3::HEPEVT_Wrapper_Runtime::t(const int) const --> double", pybind11::arg("index"));
		cl.def("number_parents", (int (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::number_parents, "C++: HepMC3::HEPEVT_Wrapper_Runtime::number_parents(const int) const --> int", pybind11::arg("index"));
		cl.def("number_children", (int (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::number_children, "C++: HepMC3::HEPEVT_Wrapper_Runtime::number_children(const int) const --> int", pybind11::arg("index"));
		cl.def("number_children_exact", (int (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int) const) &HepMC3::HEPEVT_Wrapper_Runtime::number_children_exact, "C++: HepMC3::HEPEVT_Wrapper_Runtime::number_children_exact(const int) const --> int", pybind11::arg("index"));
		cl.def("set_event_number", (void (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int)) &HepMC3::HEPEVT_Wrapper_Runtime::set_event_number, "C++: HepMC3::HEPEVT_Wrapper_Runtime::set_event_number(const int) --> void", pybind11::arg("evtno"));
		cl.def("set_number_entries", (void (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int)) &HepMC3::HEPEVT_Wrapper_Runtime::set_number_entries, "C++: HepMC3::HEPEVT_Wrapper_Runtime::set_number_entries(const int) --> void", pybind11::arg("noentries"));
		cl.def("set_status", (void (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int, const int)) &HepMC3::HEPEVT_Wrapper_Runtime::set_status, "C++: HepMC3::HEPEVT_Wrapper_Runtime::set_status(const int, const int) --> void", pybind11::arg("index"), pybind11::arg("status"));
		cl.def("set_id", (void (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int, const int)) &HepMC3::HEPEVT_Wrapper_Runtime::set_id, "C++: HepMC3::HEPEVT_Wrapper_Runtime::set_id(const int, const int) --> void", pybind11::arg("index"), pybind11::arg("id"));
		cl.def("set_parents", (void (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int, const int, const int)) &HepMC3::HEPEVT_Wrapper_Runtime::set_parents, "C++: HepMC3::HEPEVT_Wrapper_Runtime::set_parents(const int, const int, const int) --> void", pybind11::arg("index"), pybind11::arg("firstparent"), pybind11::arg("lastparent"));
		cl.def("set_children", (void (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int, const int, const int)) &HepMC3::HEPEVT_Wrapper_Runtime::set_children, "C++: HepMC3::HEPEVT_Wrapper_Runtime::set_children(const int, const int, const int) --> void", pybind11::arg("index"), pybind11::arg("firstchild"), pybind11::arg("lastchild"));
		cl.def("set_momentum", (void (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int, const double, const double, const double, const double)) &HepMC3::HEPEVT_Wrapper_Runtime::set_momentum, "C++: HepMC3::HEPEVT_Wrapper_Runtime::set_momentum(const int, const double, const double, const double, const double) --> void", pybind11::arg("index"), pybind11::arg("px"), pybind11::arg("py"), pybind11::arg("pz"), pybind11::arg("e"));
		cl.def("set_mass", (void (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int, double)) &HepMC3::HEPEVT_Wrapper_Runtime::set_mass, "C++: HepMC3::HEPEVT_Wrapper_Runtime::set_mass(const int, double) --> void", pybind11::arg("index"), pybind11::arg("mass"));
		cl.def("set_position", (void (HepMC3::HEPEVT_Wrapper_Runtime::*)(const int, const double, const double, const double, const double)) &HepMC3::HEPEVT_Wrapper_Runtime::set_position, "C++: HepMC3::HEPEVT_Wrapper_Runtime::set_position(const int, const double, const double, const double, const double) --> void", pybind11::arg("index"), pybind11::arg("x"), pybind11::arg("y"), pybind11::arg("z"), pybind11::arg("t"));
		cl.def("assign", (class HepMC3::HEPEVT_Wrapper_Runtime & (HepMC3::HEPEVT_Wrapper_Runtime::*)(const class HepMC3::HEPEVT_Wrapper_Runtime &)) &HepMC3::HEPEVT_Wrapper_Runtime::operator=, "C++: HepMC3::HEPEVT_Wrapper_Runtime::operator=(const class HepMC3::HEPEVT_Wrapper_Runtime &) --> class HepMC3::HEPEVT_Wrapper_Runtime &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_HEPEVT_Wrapper_Runtime_binder(cl);
	}
}
