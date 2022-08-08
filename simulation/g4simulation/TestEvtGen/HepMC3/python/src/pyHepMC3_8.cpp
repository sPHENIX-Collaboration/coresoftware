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

// HepMC3::GenPdfInfo file: line:32
struct PyCallBack_HepMC3_GenPdfInfo : public HepMC3::GenPdfInfo {
	using HepMC3::GenPdfInfo::GenPdfInfo;

	bool from_string(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::GenPdfInfo *>(this), "from_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return GenPdfInfo::from_string(a0);
	}
	bool to_string(std::string & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::GenPdfInfo *>(this), "to_string");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return GenPdfInfo::to_string(a0);
	}
	bool init() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::GenPdfInfo *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::GenPdfInfo *>(this), "init");
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

void bind_pyHepMC3_8(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // HepMC3::GenRunInfo file:HepMC3/GenRunInfo.h line:33
		pybind11::class_<HepMC3::GenRunInfo, std::shared_ptr<HepMC3::GenRunInfo>> cl(M("HepMC3"), "GenRunInfo", "Stores run-related information\n\n Manages run-related information.\n Contains run-wide attributes");
		cl.def( pybind11::init( [](){ return new HepMC3::GenRunInfo(); } ) );
		cl.def( pybind11::init( [](HepMC3::GenRunInfo const &o){ return new HepMC3::GenRunInfo(o); } ) );
		cl.def("attribute", (class std::shared_ptr<class HepMC3::GenHeavyIon> (HepMC3::GenRunInfo::*)(const std::string &) const) &HepMC3::GenRunInfo::attribute<HepMC3::GenHeavyIon>, "C++: HepMC3::GenRunInfo::attribute(const std::string &) const --> class std::shared_ptr<class HepMC3::GenHeavyIon>", pybind11::arg("name"));
		cl.def("attribute", (class std::shared_ptr<class HepMC3::GenPdfInfo> (HepMC3::GenRunInfo::*)(const std::string &) const) &HepMC3::GenRunInfo::attribute<HepMC3::GenPdfInfo>, "C++: HepMC3::GenRunInfo::attribute(const std::string &) const --> class std::shared_ptr<class HepMC3::GenPdfInfo>", pybind11::arg("name"));
		cl.def("attribute", (class std::shared_ptr<class HepMC3::GenCrossSection> (HepMC3::GenRunInfo::*)(const std::string &) const) &HepMC3::GenRunInfo::attribute<HepMC3::GenCrossSection>, "C++: HepMC3::GenRunInfo::attribute(const std::string &) const --> class std::shared_ptr<class HepMC3::GenCrossSection>", pybind11::arg("name"));
		cl.def("assign", (class HepMC3::GenRunInfo & (HepMC3::GenRunInfo::*)(const class HepMC3::GenRunInfo &)) &HepMC3::GenRunInfo::operator=, "Assignmet\n\nC++: HepMC3::GenRunInfo::operator=(const class HepMC3::GenRunInfo &) --> class HepMC3::GenRunInfo &", pybind11::return_value_policy::automatic, pybind11::arg("r"));
		cl.def("has_weight", (bool (HepMC3::GenRunInfo::*)(const std::string &) const) &HepMC3::GenRunInfo::has_weight, "Check if a weight name is present.\n\nC++: HepMC3::GenRunInfo::has_weight(const std::string &) const --> bool", pybind11::arg("name"));
		cl.def("weight_indices", (class std::map<std::string, int> (HepMC3::GenRunInfo::*)() const) &HepMC3::GenRunInfo::weight_indices, "Returns a copy of indices map.\n\nC++: HepMC3::GenRunInfo::weight_indices() const --> class std::map<std::string, int>");
		cl.def("weight_index", (int (HepMC3::GenRunInfo::*)(const std::string &) const) &HepMC3::GenRunInfo::weight_index, "Return the index corresponding to a weight name.\n \n\n -1 if name was not found\n\nC++: HepMC3::GenRunInfo::weight_index(const std::string &) const --> int", pybind11::arg("name"));
		cl.def("weight_names", (const class std::vector<std::string > & (HepMC3::GenRunInfo::*)() const) &HepMC3::GenRunInfo::weight_names, "Get the vector of weight names.\n\nC++: HepMC3::GenRunInfo::weight_names() const --> const class std::vector<std::string > &", pybind11::return_value_policy::automatic);
		cl.def("set_weight_names", (void (HepMC3::GenRunInfo::*)(const class std::vector<std::string > &)) &HepMC3::GenRunInfo::set_weight_names, "Set the names of the weights in this run.\n\n For consistency, the length of the vector should be the same as\n the number of weights in the events in the run.\n\nC++: HepMC3::GenRunInfo::set_weight_names(const class std::vector<std::string > &) --> void", pybind11::arg("names"));
		cl.def("add_attribute", (void (HepMC3::GenRunInfo::*)(const std::string &, const class std::shared_ptr<class HepMC3::Attribute> &)) &HepMC3::GenRunInfo::add_attribute, "add an attribute\n This will overwrite existing attribute if an attribute\n with the same name is present\n\nC++: HepMC3::GenRunInfo::add_attribute(const std::string &, const class std::shared_ptr<class HepMC3::Attribute> &) --> void", pybind11::arg("name"), pybind11::arg("att"));
		cl.def("remove_attribute", (void (HepMC3::GenRunInfo::*)(const std::string &)) &HepMC3::GenRunInfo::remove_attribute, "Remove attribute\n\nC++: HepMC3::GenRunInfo::remove_attribute(const std::string &) --> void", pybind11::arg("name"));
		cl.def("attribute_as_string", (std::string (HepMC3::GenRunInfo::*)(const std::string &) const) &HepMC3::GenRunInfo::attribute_as_string, "Get attribute of any type as string\n\nC++: HepMC3::GenRunInfo::attribute_as_string(const std::string &) const --> std::string", pybind11::arg("name"));
		cl.def("attribute_names", (class std::vector<std::string > (HepMC3::GenRunInfo::*)() const) &HepMC3::GenRunInfo::attribute_names, "Get list of attribute names\n\nC++: HepMC3::GenRunInfo::attribute_names() const --> class std::vector<std::string >");
		cl.def("attributes", (class std::map<std::string, class std::shared_ptr<class HepMC3::Attribute> > (HepMC3::GenRunInfo::*)() const) &HepMC3::GenRunInfo::attributes, "Get a copy of the list of attributes\n \n\n To avoid thread issues, this is returns a copy. Better solution may be needed.\n\nC++: HepMC3::GenRunInfo::attributes() const --> class std::map<std::string, class std::shared_ptr<class HepMC3::Attribute> >");
		cl.def("write_data", (void (HepMC3::GenRunInfo::*)(struct HepMC3::GenRunInfoData &) const) &HepMC3::GenRunInfo::write_data, "Fill GenRunInfoData object\n\nC++: HepMC3::GenRunInfo::write_data(struct HepMC3::GenRunInfoData &) const --> void", pybind11::arg("data"));
		cl.def("read_data", (void (HepMC3::GenRunInfo::*)(const struct HepMC3::GenRunInfoData &)) &HepMC3::GenRunInfo::read_data, "Fill GenRunInfo based on GenRunInfoData\n\nC++: HepMC3::GenRunInfo::read_data(const struct HepMC3::GenRunInfoData &) --> void", pybind11::arg("data"));

		 binder::custom_GenRunInfo_binder(cl);

		{ // HepMC3::GenRunInfo::ToolInfo file:HepMC3/GenRunInfo.h line:38
			auto & enclosing_class = cl;
			pybind11::class_<HepMC3::GenRunInfo::ToolInfo, std::shared_ptr<HepMC3::GenRunInfo::ToolInfo>> cl(enclosing_class, "ToolInfo", "Interrnal struct for keeping track of tools.");
			cl.def( pybind11::init( [](HepMC3::GenRunInfo::ToolInfo const &o){ return new HepMC3::GenRunInfo::ToolInfo(o); } ) );
			cl.def( pybind11::init( [](){ return new HepMC3::GenRunInfo::ToolInfo(); } ) );
			cl.def_readwrite("name", &HepMC3::GenRunInfo::ToolInfo::name);
			cl.def_readwrite("version", &HepMC3::GenRunInfo::ToolInfo::version);
			cl.def_readwrite("description", &HepMC3::GenRunInfo::ToolInfo::description);
			cl.def("assign", (struct HepMC3::GenRunInfo::ToolInfo & (HepMC3::GenRunInfo::ToolInfo::*)(const struct HepMC3::GenRunInfo::ToolInfo &)) &HepMC3::GenRunInfo::ToolInfo::operator=, "C++: HepMC3::GenRunInfo::ToolInfo::operator=(const struct HepMC3::GenRunInfo::ToolInfo &) --> struct HepMC3::GenRunInfo::ToolInfo &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		}

	}
	{ // HepMC3::GenParticleData file:HepMC3/Data/GenParticleData.h line:24
		pybind11::class_<HepMC3::GenParticleData, std::shared_ptr<HepMC3::GenParticleData>> cl(M("HepMC3"), "GenParticleData", "");
		cl.def( pybind11::init( [](){ return new HepMC3::GenParticleData(); } ) );
		cl.def( pybind11::init( [](HepMC3::GenParticleData const &o){ return new HepMC3::GenParticleData(o); } ) );
		cl.def_readwrite("pid", &HepMC3::GenParticleData::pid);
		cl.def_readwrite("status", &HepMC3::GenParticleData::status);
		cl.def_readwrite("is_mass_set", &HepMC3::GenParticleData::is_mass_set);
		cl.def_readwrite("mass", &HepMC3::GenParticleData::mass);
		cl.def_readwrite("momentum", &HepMC3::GenParticleData::momentum);
		cl.def("assign", (struct HepMC3::GenParticleData & (HepMC3::GenParticleData::*)(const struct HepMC3::GenParticleData &)) &HepMC3::GenParticleData::operator=, "C++: HepMC3::GenParticleData::operator=(const struct HepMC3::GenParticleData &) --> struct HepMC3::GenParticleData &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::GenPdfInfo file: line:32
		pybind11::class_<HepMC3::GenPdfInfo, std::shared_ptr<HepMC3::GenPdfInfo>, PyCallBack_HepMC3_GenPdfInfo, HepMC3::Attribute> cl(M("HepMC3"), "GenPdfInfo", "");
		cl.def( pybind11::init( [](PyCallBack_HepMC3_GenPdfInfo const &o){ return new PyCallBack_HepMC3_GenPdfInfo(o); } ) );
		cl.def( pybind11::init( [](HepMC3::GenPdfInfo const &o){ return new HepMC3::GenPdfInfo(o); } ) );
		cl.def( pybind11::init( [](){ return new HepMC3::GenPdfInfo(); }, [](){ return new PyCallBack_HepMC3_GenPdfInfo(); } ) );
		cl.def_readwrite("scale", &HepMC3::GenPdfInfo::scale);
		cl.def("from_string", (bool (HepMC3::GenPdfInfo::*)(const std::string &)) &HepMC3::GenPdfInfo::from_string, "Implementation of Attribute::from_string \n\nC++: HepMC3::GenPdfInfo::from_string(const std::string &) --> bool", pybind11::arg("att"));
		cl.def("to_string", (bool (HepMC3::GenPdfInfo::*)(std::string &) const) &HepMC3::GenPdfInfo::to_string, "Implementation of Attribute::to_string \n\nC++: HepMC3::GenPdfInfo::to_string(std::string &) const --> bool", pybind11::arg("att"));
		cl.def("set", [](HepMC3::GenPdfInfo &o, const int & a0, const int & a1, const double & a2, const double & a3, const double & a4, const double & a5, const double & a6) -> void { return o.set(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("parton_id1"), pybind11::arg("parton_id2"), pybind11::arg("x1"), pybind11::arg("x2"), pybind11::arg("scale_in"), pybind11::arg("xf1"), pybind11::arg("xf2"));
		cl.def("set", [](HepMC3::GenPdfInfo &o, const int & a0, const int & a1, const double & a2, const double & a3, const double & a4, const double & a5, const double & a6, const int & a7) -> void { return o.set(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("parton_id1"), pybind11::arg("parton_id2"), pybind11::arg("x1"), pybind11::arg("x2"), pybind11::arg("scale_in"), pybind11::arg("xf1"), pybind11::arg("xf2"), pybind11::arg("pdf_id1"));
		cl.def("set", (void (HepMC3::GenPdfInfo::*)(const int &, const int &, const double &, const double &, const double &, const double &, const double &, const int &, const int &)) &HepMC3::GenPdfInfo::set, "Set all fields \n\nC++: HepMC3::GenPdfInfo::set(const int &, const int &, const double &, const double &, const double &, const double &, const double &, const int &, const int &) --> void", pybind11::arg("parton_id1"), pybind11::arg("parton_id2"), pybind11::arg("x1"), pybind11::arg("x2"), pybind11::arg("scale_in"), pybind11::arg("xf1"), pybind11::arg("xf2"), pybind11::arg("pdf_id1"), pybind11::arg("pdf_id2"));
		cl.def("__eq__", (bool (HepMC3::GenPdfInfo::*)(const class HepMC3::GenPdfInfo &) const) &HepMC3::GenPdfInfo::operator==, "C++: HepMC3::GenPdfInfo::operator==(const class HepMC3::GenPdfInfo &) const --> bool", pybind11::arg(""));
		cl.def("__ne__", (bool (HepMC3::GenPdfInfo::*)(const class HepMC3::GenPdfInfo &) const) &HepMC3::GenPdfInfo::operator!=, "C++: HepMC3::GenPdfInfo::operator!=(const class HepMC3::GenPdfInfo &) const --> bool", pybind11::arg(""));
		cl.def("is_valid", (bool (HepMC3::GenPdfInfo::*)() const) &HepMC3::GenPdfInfo::is_valid, "C++: HepMC3::GenPdfInfo::is_valid() const --> bool");
		cl.def("assign", (class HepMC3::GenPdfInfo & (HepMC3::GenPdfInfo::*)(const class HepMC3::GenPdfInfo &)) &HepMC3::GenPdfInfo::operator=, "C++: HepMC3::GenPdfInfo::operator=(const class HepMC3::GenPdfInfo &) --> class HepMC3::GenPdfInfo &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::GenEvent file:HepMC3/GenEvent.h line:41
		pybind11::class_<HepMC3::GenEvent, std::shared_ptr<HepMC3::GenEvent>> cl(M("HepMC3"), "GenEvent", "Stores event-related information\n\n Manages event-related information.\n Contains lists of GenParticle and GenVertex objects");
		cl.def( pybind11::init( [](){ return new HepMC3::GenEvent(); } ), "doc" );
		cl.def( pybind11::init( [](enum HepMC3::Units::MomentumUnit const & a0){ return new HepMC3::GenEvent(a0); } ), "doc" , pybind11::arg("momentum_unit"));
		cl.def( pybind11::init<enum HepMC3::Units::MomentumUnit, enum HepMC3::Units::LengthUnit>(), pybind11::arg("momentum_unit"), pybind11::arg("length_unit") );

		cl.def( pybind11::init( [](class std::shared_ptr<class HepMC3::GenRunInfo> const & a0){ return new HepMC3::GenEvent(a0); } ), "doc" , pybind11::arg("run"));
		cl.def( pybind11::init( [](class std::shared_ptr<class HepMC3::GenRunInfo> const & a0, enum HepMC3::Units::MomentumUnit const & a1){ return new HepMC3::GenEvent(a0, a1); } ), "doc" , pybind11::arg("run"), pybind11::arg("momentum_unit"));
		cl.def( pybind11::init<class std::shared_ptr<class HepMC3::GenRunInfo>, enum HepMC3::Units::MomentumUnit, enum HepMC3::Units::LengthUnit>(), pybind11::arg("run"), pybind11::arg("momentum_unit"), pybind11::arg("length_unit") );

		cl.def( pybind11::init( [](HepMC3::GenEvent const &o){ return new HepMC3::GenEvent(o); } ) );
		cl.def("assign", (class HepMC3::GenEvent & (HepMC3::GenEvent::*)(const class HepMC3::GenEvent &)) &HepMC3::GenEvent::operator=, "Assignment operator\n\nC++: HepMC3::GenEvent::operator=(const class HepMC3::GenEvent &) --> class HepMC3::GenEvent &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("particles", (const class std::vector<class std::shared_ptr<class HepMC3::GenParticle> > & (HepMC3::GenEvent::*)()) &HepMC3::GenEvent::particles, "Get/set list of particles (non-const)\n\nC++: HepMC3::GenEvent::particles() --> const class std::vector<class std::shared_ptr<class HepMC3::GenParticle> > &", pybind11::return_value_policy::automatic);
		cl.def("vertices", (const class std::vector<class std::shared_ptr<class HepMC3::GenVertex> > & (HepMC3::GenEvent::*)()) &HepMC3::GenEvent::vertices, "Get/set list of vertices (non-const)\n\nC++: HepMC3::GenEvent::vertices() --> const class std::vector<class std::shared_ptr<class HepMC3::GenVertex> > &", pybind11::return_value_policy::automatic);
		cl.def("particles_size", (int (HepMC3::GenEvent::*)() const) &HepMC3::GenEvent::particles_size, "Particles size, HepMC2 compatiility\n\nC++: HepMC3::GenEvent::particles_size() const --> int");
		cl.def("particles_empty", (bool (HepMC3::GenEvent::*)() const) &HepMC3::GenEvent::particles_empty, "Particles empty, HepMC2 compatiility\n\nC++: HepMC3::GenEvent::particles_empty() const --> bool");
		cl.def("vertices_size", (int (HepMC3::GenEvent::*)() const) &HepMC3::GenEvent::vertices_size, "Vertices size, HepMC2 compatiility\n\nC++: HepMC3::GenEvent::vertices_size() const --> int");
		cl.def("vertices_empty", (bool (HepMC3::GenEvent::*)() const) &HepMC3::GenEvent::vertices_empty, "Vertices empty, HepMC2 compatiility\n\nC++: HepMC3::GenEvent::vertices_empty() const --> bool");
		cl.def("weights", (class std::vector<double> & (HepMC3::GenEvent::*)()) &HepMC3::GenEvent::weights, "Get event weights as a vector (non-const)\n\nC++: HepMC3::GenEvent::weights() --> class std::vector<double> &", pybind11::return_value_policy::automatic);
		cl.def("weight", [](HepMC3::GenEvent const &o) -> double { return o.weight(); }, "");
		cl.def("weight", (double (HepMC3::GenEvent::*)(const unsigned long &) const) &HepMC3::GenEvent::weight, "Get event weight accessed by index (or the canonical/first one if there is no argument)\n \n\n It's the user's responsibility to ensure that the given index exists!\n\nC++: HepMC3::GenEvent::weight(const unsigned long &) const --> double", pybind11::arg("index"));
		cl.def("weight", (double & (HepMC3::GenEvent::*)(const std::string &)) &HepMC3::GenEvent::weight, "Get event weight accessed by weight name\n \n\n Requires there to be an attached GenRunInfo, otherwise will throw an exception\n \n\n It's the user's responsibility to ensure that the given name exists!\n\nC++: HepMC3::GenEvent::weight(const std::string &) --> double &", pybind11::return_value_policy::automatic, pybind11::arg("name"));
		cl.def("weight_names", (const class std::vector<std::string > & (HepMC3::GenEvent::*)() const) &HepMC3::GenEvent::weight_names, "Get event weight names, if there are some\n \n\n Requires there to be an attached GenRunInfo with registered weight names, otherwise will throw an exception\n\nC++: HepMC3::GenEvent::weight_names() const --> const class std::vector<std::string > &", pybind11::return_value_policy::automatic);
		cl.def("run_info", (class std::shared_ptr<class HepMC3::GenRunInfo> (HepMC3::GenEvent::*)() const) &HepMC3::GenEvent::run_info, "Get a pointer to the the GenRunInfo object.\n\nC++: HepMC3::GenEvent::run_info() const --> class std::shared_ptr<class HepMC3::GenRunInfo>");
		cl.def("set_run_info", (void (HepMC3::GenEvent::*)(class std::shared_ptr<class HepMC3::GenRunInfo>)) &HepMC3::GenEvent::set_run_info, "Set the GenRunInfo object by smart pointer.\n\nC++: HepMC3::GenEvent::set_run_info(class std::shared_ptr<class HepMC3::GenRunInfo>) --> void", pybind11::arg("run"));
		cl.def("event_number", (int (HepMC3::GenEvent::*)() const) &HepMC3::GenEvent::event_number, "Get event number\n\nC++: HepMC3::GenEvent::event_number() const --> int");
		cl.def("set_event_number", (void (HepMC3::GenEvent::*)(const int &)) &HepMC3::GenEvent::set_event_number, "Set event number\n\nC++: HepMC3::GenEvent::set_event_number(const int &) --> void", pybind11::arg("num"));
		cl.def("momentum_unit", (const enum HepMC3::Units::MomentumUnit & (HepMC3::GenEvent::*)() const) &HepMC3::GenEvent::momentum_unit, "Get momentum unit\n\nC++: HepMC3::GenEvent::momentum_unit() const --> const enum HepMC3::Units::MomentumUnit &", pybind11::return_value_policy::automatic);
		cl.def("length_unit", (const enum HepMC3::Units::LengthUnit & (HepMC3::GenEvent::*)() const) &HepMC3::GenEvent::length_unit, "Get length unit\n\nC++: HepMC3::GenEvent::length_unit() const --> const enum HepMC3::Units::LengthUnit &", pybind11::return_value_policy::automatic);
		cl.def("set_units", (void (HepMC3::GenEvent::*)(enum HepMC3::Units::MomentumUnit, enum HepMC3::Units::LengthUnit)) &HepMC3::GenEvent::set_units, "Change event units\n Converts event from current units to new ones\n\nC++: HepMC3::GenEvent::set_units(enum HepMC3::Units::MomentumUnit, enum HepMC3::Units::LengthUnit) --> void", pybind11::arg("new_momentum_unit"), pybind11::arg("new_length_unit"));
		cl.def("heavy_ion", (class std::shared_ptr<class HepMC3::GenHeavyIon> (HepMC3::GenEvent::*)()) &HepMC3::GenEvent::heavy_ion, "Get heavy ion generator additional information\n\nC++: HepMC3::GenEvent::heavy_ion() --> class std::shared_ptr<class HepMC3::GenHeavyIon>");
		cl.def("set_heavy_ion", (void (HepMC3::GenEvent::*)(class std::shared_ptr<class HepMC3::GenHeavyIon>)) &HepMC3::GenEvent::set_heavy_ion, "Set heavy ion generator additional information\n\nC++: HepMC3::GenEvent::set_heavy_ion(class std::shared_ptr<class HepMC3::GenHeavyIon>) --> void", pybind11::arg("hi"));
		cl.def("pdf_info", (class std::shared_ptr<class HepMC3::GenPdfInfo> (HepMC3::GenEvent::*)()) &HepMC3::GenEvent::pdf_info, "Get PDF information\n\nC++: HepMC3::GenEvent::pdf_info() --> class std::shared_ptr<class HepMC3::GenPdfInfo>");
		cl.def("set_pdf_info", (void (HepMC3::GenEvent::*)(class std::shared_ptr<class HepMC3::GenPdfInfo>)) &HepMC3::GenEvent::set_pdf_info, "Set PDF information\n\nC++: HepMC3::GenEvent::set_pdf_info(class std::shared_ptr<class HepMC3::GenPdfInfo>) --> void", pybind11::arg("pi"));
		cl.def("cross_section", (class std::shared_ptr<class HepMC3::GenCrossSection> (HepMC3::GenEvent::*)()) &HepMC3::GenEvent::cross_section, "Get cross-section information\n\nC++: HepMC3::GenEvent::cross_section() --> class std::shared_ptr<class HepMC3::GenCrossSection>");
		cl.def("set_cross_section", (void (HepMC3::GenEvent::*)(class std::shared_ptr<class HepMC3::GenCrossSection>)) &HepMC3::GenEvent::set_cross_section, "Set cross-section information\n\nC++: HepMC3::GenEvent::set_cross_section(class std::shared_ptr<class HepMC3::GenCrossSection>) --> void", pybind11::arg("cs"));
		cl.def("event_pos", (const class HepMC3::FourVector & (HepMC3::GenEvent::*)() const) &HepMC3::GenEvent::event_pos, "Vertex representing the overall event position\n\nC++: HepMC3::GenEvent::event_pos() const --> const class HepMC3::FourVector &", pybind11::return_value_policy::automatic);
		cl.def("beams", (const class std::vector<class std::shared_ptr<class HepMC3::GenParticle> > & (HepMC3::GenEvent::*)()) &HepMC3::GenEvent::beams, "Vector of beam particles\n\nC++: HepMC3::GenEvent::beams() --> const class std::vector<class std::shared_ptr<class HepMC3::GenParticle> > &", pybind11::return_value_policy::automatic);
		cl.def("shift_position_by", (void (HepMC3::GenEvent::*)(const class HepMC3::FourVector &)) &HepMC3::GenEvent::shift_position_by, "Shift position of all vertices in the event by \n\nC++: HepMC3::GenEvent::shift_position_by(const class HepMC3::FourVector &) --> void", pybind11::arg("delta"));
		cl.def("shift_position_to", (void (HepMC3::GenEvent::*)(const class HepMC3::FourVector &)) &HepMC3::GenEvent::shift_position_to, "Shift position of all vertices in the event to \n\nC++: HepMC3::GenEvent::shift_position_to(const class HepMC3::FourVector &) --> void", pybind11::arg("newpos"));
		cl.def("boost", (bool (HepMC3::GenEvent::*)(const class HepMC3::FourVector &)) &HepMC3::GenEvent::boost, "Boost event using x,y,z components of  as velocities\n\nC++: HepMC3::GenEvent::boost(const class HepMC3::FourVector &) --> bool", pybind11::arg("v"));
		cl.def("rotate", (bool (HepMC3::GenEvent::*)(const class HepMC3::FourVector &)) &HepMC3::GenEvent::rotate, "Rotate event using x,y,z components of  as rotation angles\n\nC++: HepMC3::GenEvent::rotate(const class HepMC3::FourVector &) --> bool", pybind11::arg("v"));
		cl.def("reflect", (bool (HepMC3::GenEvent::*)(const int)) &HepMC3::GenEvent::reflect, "Change sign of \n\nC++: HepMC3::GenEvent::reflect(const int) --> bool", pybind11::arg("axis"));
		cl.def("add_attribute", [](HepMC3::GenEvent &o, const std::string & a0, const class std::shared_ptr<class HepMC3::Attribute> & a1) -> void { return o.add_attribute(a0, a1); }, "", pybind11::arg("name"), pybind11::arg("att"));
		cl.def("add_attribute", (void (HepMC3::GenEvent::*)(const std::string &, const class std::shared_ptr<class HepMC3::Attribute> &, const int &)) &HepMC3::GenEvent::add_attribute, "@{\n \n\n Add event attribute to event\n\n This will overwrite existing attribute if an attribute\n with the same name is present\n\nC++: HepMC3::GenEvent::add_attribute(const std::string &, const class std::shared_ptr<class HepMC3::Attribute> &, const int &) --> void", pybind11::arg("name"), pybind11::arg("att"), pybind11::arg("id"));
		cl.def("remove_attribute", [](HepMC3::GenEvent &o, const std::string & a0) -> void { return o.remove_attribute(a0); }, "", pybind11::arg("name"));
		cl.def("remove_attribute", (void (HepMC3::GenEvent::*)(const std::string &, const int &)) &HepMC3::GenEvent::remove_attribute, "Remove attribute\n\nC++: HepMC3::GenEvent::remove_attribute(const std::string &, const int &) --> void", pybind11::arg("name"), pybind11::arg("id"));
		cl.def("attribute_as_string", [](HepMC3::GenEvent const &o, const std::string & a0) -> std::string { return o.attribute_as_string(a0); }, "", pybind11::arg("name"));
		cl.def("attribute_as_string", (std::string (HepMC3::GenEvent::*)(const std::string &, const int &) const) &HepMC3::GenEvent::attribute_as_string, "Get attribute of any type as string\n\nC++: HepMC3::GenEvent::attribute_as_string(const std::string &, const int &) const --> std::string", pybind11::arg("name"), pybind11::arg("id"));
		cl.def("attribute_names", [](HepMC3::GenEvent const &o) -> std::vector<std::string > { return o.attribute_names(); }, "");
		cl.def("attribute_names", (class std::vector<std::string > (HepMC3::GenEvent::*)(const int &) const) &HepMC3::GenEvent::attribute_names, "Get list of attribute names\n\nC++: HepMC3::GenEvent::attribute_names(const int &) const --> class std::vector<std::string >", pybind11::arg("id"));
		cl.def("attributes", (class std::map<std::string, class std::map<int, class std::shared_ptr<class HepMC3::Attribute> > > (HepMC3::GenEvent::*)() const) &HepMC3::GenEvent::attributes, "Get a copy of the list of attributes\n \n\n To avoid thread issues, this is returns a copy. Better solution may be needed.\n\nC++: HepMC3::GenEvent::attributes() const --> class std::map<std::string, class std::map<int, class std::shared_ptr<class HepMC3::Attribute> > >");
		cl.def("add_particle", (void (HepMC3::GenEvent::*)(class std::shared_ptr<class HepMC3::GenParticle>)) &HepMC3::GenEvent::add_particle, "Add particle\n\nC++: HepMC3::GenEvent::add_particle(class std::shared_ptr<class HepMC3::GenParticle>) --> void", pybind11::arg("p"));
		cl.def("add_vertex", (void (HepMC3::GenEvent::*)(class std::shared_ptr<class HepMC3::GenVertex>)) &HepMC3::GenEvent::add_vertex, "Add vertex\n\nC++: HepMC3::GenEvent::add_vertex(class std::shared_ptr<class HepMC3::GenVertex>) --> void", pybind11::arg("v"));
		cl.def("remove_particle", (void (HepMC3::GenEvent::*)(class std::shared_ptr<class HepMC3::GenParticle>)) &HepMC3::GenEvent::remove_particle, "Remove particle from the event\n\n This function  will remove whole sub-tree starting from this particle\n if it is the only incoming particle of this vertex.\n It will also production vertex of this particle if this vertex\n has no more outgoing particles\n\nC++: HepMC3::GenEvent::remove_particle(class std::shared_ptr<class HepMC3::GenParticle>) --> void", pybind11::arg("v"));
		cl.def("remove_particles", (void (HepMC3::GenEvent::*)(class std::vector<class std::shared_ptr<class HepMC3::GenParticle> >)) &HepMC3::GenEvent::remove_particles, "Remove a set of particles\n\n This function follows rules of GenEvent::remove_particle to remove\n a list of particles from the event.\n\nC++: HepMC3::GenEvent::remove_particles(class std::vector<class std::shared_ptr<class HepMC3::GenParticle> >) --> void", pybind11::arg("v"));
		cl.def("remove_vertex", (void (HepMC3::GenEvent::*)(class std::shared_ptr<class HepMC3::GenVertex>)) &HepMC3::GenEvent::remove_vertex, "Remove vertex from the event\n\n This will remove all sub-trees of all outgoing particles of this vertex\n\nC++: HepMC3::GenEvent::remove_vertex(class std::shared_ptr<class HepMC3::GenVertex>) --> void", pybind11::arg("v"));
		cl.def("add_tree", (void (HepMC3::GenEvent::*)(const class std::vector<class std::shared_ptr<class HepMC3::GenParticle> > &)) &HepMC3::GenEvent::add_tree, "Add whole tree in topological order\n\n This function will find the beam particles (particles\n that have no production vertices or their production vertices\n have no particles) and will add the whole decay tree starting from\n these particles.\n\n \n Any particles on this list that do not belong to the tree\n       will be ignored.\n\nC++: HepMC3::GenEvent::add_tree(const class std::vector<class std::shared_ptr<class HepMC3::GenParticle> > &) --> void", pybind11::arg("particles"));
		cl.def("reserve", [](HepMC3::GenEvent &o, const unsigned long & a0) -> void { return o.reserve(a0); }, "", pybind11::arg("particles"));
		cl.def("reserve", (void (HepMC3::GenEvent::*)(const unsigned long &, const unsigned long &)) &HepMC3::GenEvent::reserve, "Reserve memory for particles and vertices\n\n Helps optimize event creation when size of the event is known beforehand\n\nC++: HepMC3::GenEvent::reserve(const unsigned long &, const unsigned long &) --> void", pybind11::arg("particles"), pybind11::arg("vertices"));
		cl.def("clear", (void (HepMC3::GenEvent::*)()) &HepMC3::GenEvent::clear, "Remove contents of this event\n\nC++: HepMC3::GenEvent::clear() --> void");
		cl.def("add_particle", (void (HepMC3::GenEvent::*)(class HepMC3::GenParticle *)) &HepMC3::GenEvent::add_particle, "Add particle by raw pointer\n \n\n Use GenEvent::add_particle( const GenParticlePtr& ) instead\n\nC++: HepMC3::GenEvent::add_particle(class HepMC3::GenParticle *) --> void", pybind11::arg("p"));
		cl.def("add_vertex", (void (HepMC3::GenEvent::*)(class HepMC3::GenVertex *)) &HepMC3::GenEvent::add_vertex, "Add vertex by raw pointer\n \n\n Use GenEvent::add_vertex( const GenVertexPtr& ) instead\n\nC++: HepMC3::GenEvent::add_vertex(class HepMC3::GenVertex *) --> void", pybind11::arg("v"));
		cl.def("set_beam_particles", (void (HepMC3::GenEvent::*)(class std::shared_ptr<class HepMC3::GenParticle>, class std::shared_ptr<class HepMC3::GenParticle>)) &HepMC3::GenEvent::set_beam_particles, "Set incoming beam particles\n \n\n Backward compatibility\n\nC++: HepMC3::GenEvent::set_beam_particles(class std::shared_ptr<class HepMC3::GenParticle>, class std::shared_ptr<class HepMC3::GenParticle>) --> void", pybind11::arg("p1"), pybind11::arg("p2"));
		cl.def("add_beam_particle", (void (HepMC3::GenEvent::*)(class std::shared_ptr<class HepMC3::GenParticle>)) &HepMC3::GenEvent::add_beam_particle, "Add  particle to root vertex\n\nC++: HepMC3::GenEvent::add_beam_particle(class std::shared_ptr<class HepMC3::GenParticle>) --> void", pybind11::arg("p1"));
		cl.def("write_data", (void (HepMC3::GenEvent::*)(struct HepMC3::GenEventData &) const) &HepMC3::GenEvent::write_data, "Fill GenEventData object\n\nC++: HepMC3::GenEvent::write_data(struct HepMC3::GenEventData &) const --> void", pybind11::arg("data"));
		cl.def("read_data", (void (HepMC3::GenEvent::*)(const struct HepMC3::GenEventData &)) &HepMC3::GenEvent::read_data, "Fill GenEvent based on GenEventData\n\nC++: HepMC3::GenEvent::read_data(const struct HepMC3::GenEventData &) --> void", pybind11::arg("data"));

		 binder::custom_GenEvent_binder(cl);
	}
}
