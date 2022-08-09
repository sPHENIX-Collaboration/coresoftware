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
#include <HepMC3/Writer.h>
#include <functional>
#include <ios>
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

// HepMC3::Writer file:HepMC3/Writer.h line:25
struct PyCallBack_HepMC3_Writer : public HepMC3::Writer {
	using HepMC3::Writer::Writer;

	void write_event(const class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::Writer *>(this), "write_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Writer::write_event\"");
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::Writer *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Writer::failed\"");
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::Writer *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Writer::close\"");
	}
	void set_run_info(class std::shared_ptr<class HepMC3::GenRunInfo> a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::Writer *>(this), "set_run_info");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::Writer *>(this), "run_info");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::Writer *>(this), "set_options");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::Writer *>(this), "get_options");
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

void bind_pyHepMC3_9(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // HepMC3::GenParticle file:HepMC3/GenParticle.h line:31
		pybind11::class_<HepMC3::GenParticle, std::shared_ptr<HepMC3::GenParticle>> cl(M("HepMC3"), "GenParticle", "");
		cl.def( pybind11::init( [](){ return new HepMC3::GenParticle(); } ), "doc" );
		cl.def( pybind11::init( [](const class HepMC3::FourVector & a0){ return new HepMC3::GenParticle(a0); } ), "doc" , pybind11::arg("momentum"));
		cl.def( pybind11::init( [](const class HepMC3::FourVector & a0, int const & a1){ return new HepMC3::GenParticle(a0, a1); } ), "doc" , pybind11::arg("momentum"), pybind11::arg("pid"));
		cl.def( pybind11::init<const class HepMC3::FourVector &, int, int>(), pybind11::arg("momentum"), pybind11::arg("pid"), pybind11::arg("status") );

		cl.def( pybind11::init<const struct HepMC3::GenParticleData &>(), pybind11::arg("data") );

		cl.def( pybind11::init( [](HepMC3::GenParticle const &o){ return new HepMC3::GenParticle(o); } ) );
		cl.def("in_event", (bool (HepMC3::GenParticle::*)() const) &HepMC3::GenParticle::in_event, "Check if this particle belongs to an event \n\nC++: HepMC3::GenParticle::in_event() const --> bool");
		cl.def("parent_event", (class HepMC3::GenEvent * (HepMC3::GenParticle::*)()) &HepMC3::GenParticle::parent_event, "C++: HepMC3::GenParticle::parent_event() --> class HepMC3::GenEvent *", pybind11::return_value_policy::automatic);
		cl.def("id", (int (HepMC3::GenParticle::*)() const) &HepMC3::GenParticle::id, "C++: HepMC3::GenParticle::id() const --> int");
		cl.def("data", (const struct HepMC3::GenParticleData & (HepMC3::GenParticle::*)() const) &HepMC3::GenParticle::data, "C++: HepMC3::GenParticle::data() const --> const struct HepMC3::GenParticleData &", pybind11::return_value_policy::automatic);
		cl.def("production_vertex", (class std::shared_ptr<class HepMC3::GenVertex> (HepMC3::GenParticle::*)()) &HepMC3::GenParticle::production_vertex, "C++: HepMC3::GenParticle::production_vertex() --> class std::shared_ptr<class HepMC3::GenVertex>");
		cl.def("end_vertex", (class std::shared_ptr<class HepMC3::GenVertex> (HepMC3::GenParticle::*)()) &HepMC3::GenParticle::end_vertex, "C++: HepMC3::GenParticle::end_vertex() --> class std::shared_ptr<class HepMC3::GenVertex>");
		cl.def("parents", (class std::vector<class std::shared_ptr<class HepMC3::GenParticle> > (HepMC3::GenParticle::*)()) &HepMC3::GenParticle::parents, "Convenience access to immediate incoming particles via production vertex\n \n\n Less efficient than via the vertex since return must be by value (in case there is no vertex)\n\nC++: HepMC3::GenParticle::parents() --> class std::vector<class std::shared_ptr<class HepMC3::GenParticle> >");
		cl.def("children", (class std::vector<class std::shared_ptr<class HepMC3::GenParticle> > (HepMC3::GenParticle::*)()) &HepMC3::GenParticle::children, "Convenience access to immediate outgoing particles via end vertex\n \n\n Less efficient than via the vertex since return must be by value (in case there is no vertex)\n\nC++: HepMC3::GenParticle::children() --> class std::vector<class std::shared_ptr<class HepMC3::GenParticle> >");
		cl.def("pid", (int (HepMC3::GenParticle::*)() const) &HepMC3::GenParticle::pid, "C++: HepMC3::GenParticle::pid() const --> int");
		cl.def("abs_pid", (int (HepMC3::GenParticle::*)() const) &HepMC3::GenParticle::abs_pid, "C++: HepMC3::GenParticle::abs_pid() const --> int");
		cl.def("status", (int (HepMC3::GenParticle::*)() const) &HepMC3::GenParticle::status, "C++: HepMC3::GenParticle::status() const --> int");
		cl.def("momentum", (const class HepMC3::FourVector & (HepMC3::GenParticle::*)() const) &HepMC3::GenParticle::momentum, "C++: HepMC3::GenParticle::momentum() const --> const class HepMC3::FourVector &", pybind11::return_value_policy::automatic);
		cl.def("is_generated_mass_set", (bool (HepMC3::GenParticle::*)() const) &HepMC3::GenParticle::is_generated_mass_set, "C++: HepMC3::GenParticle::is_generated_mass_set() const --> bool");
		cl.def("generated_mass", (double (HepMC3::GenParticle::*)() const) &HepMC3::GenParticle::generated_mass, "Get generated mass\n\n This function will return mass as set by a generator/tool.\n If not set, it will return momentum().m()\n\nC++: HepMC3::GenParticle::generated_mass() const --> double");
		cl.def("set_pid", (void (HepMC3::GenParticle::*)(int)) &HepMC3::GenParticle::set_pid, "C++: HepMC3::GenParticle::set_pid(int) --> void", pybind11::arg("pid"));
		cl.def("set_status", (void (HepMC3::GenParticle::*)(int)) &HepMC3::GenParticle::set_status, "C++: HepMC3::GenParticle::set_status(int) --> void", pybind11::arg("status"));
		cl.def("set_momentum", (void (HepMC3::GenParticle::*)(const class HepMC3::FourVector &)) &HepMC3::GenParticle::set_momentum, "C++: HepMC3::GenParticle::set_momentum(const class HepMC3::FourVector &) --> void", pybind11::arg("momentum"));
		cl.def("set_generated_mass", (void (HepMC3::GenParticle::*)(double)) &HepMC3::GenParticle::set_generated_mass, "C++: HepMC3::GenParticle::set_generated_mass(double) --> void", pybind11::arg("m"));
		cl.def("unset_generated_mass", (void (HepMC3::GenParticle::*)()) &HepMC3::GenParticle::unset_generated_mass, "C++: HepMC3::GenParticle::unset_generated_mass() --> void");
		cl.def("add_attribute", (bool (HepMC3::GenParticle::*)(const std::string &, class std::shared_ptr<class HepMC3::Attribute>)) &HepMC3::GenParticle::add_attribute, "Add an attribute to this particle\n\n  This will overwrite existing attribute if an attribute with\n  the same name is present. The attribute will be stored in the\n  parent_event(). \n\n false if there is no parent_event();\n\nC++: HepMC3::GenParticle::add_attribute(const std::string &, class std::shared_ptr<class HepMC3::Attribute>) --> bool", pybind11::arg("name"), pybind11::arg("att"));
		cl.def("attribute_names", (class std::vector<std::string > (HepMC3::GenParticle::*)() const) &HepMC3::GenParticle::attribute_names, "Get list of names of attributes assigned to this particle\n\nC++: HepMC3::GenParticle::attribute_names() const --> class std::vector<std::string >");
		cl.def("remove_attribute", (void (HepMC3::GenParticle::*)(const std::string &)) &HepMC3::GenParticle::remove_attribute, "Remove attribute\n\nC++: HepMC3::GenParticle::remove_attribute(const std::string &) --> void", pybind11::arg("name"));
		cl.def("attribute_as_string", (std::string (HepMC3::GenParticle::*)(const std::string &) const) &HepMC3::GenParticle::attribute_as_string, "Get attribute of any type as string\n\nC++: HepMC3::GenParticle::attribute_as_string(const std::string &) const --> std::string", pybind11::arg("name"));
		cl.def("pdg_id", (int (HepMC3::GenParticle::*)() const) &HepMC3::GenParticle::pdg_id, "Get PDG ID\n \n\n Use pid() instead\n\nC++: HepMC3::GenParticle::pdg_id() const --> int");
		cl.def("set_pdg_id", (void (HepMC3::GenParticle::*)(const int &)) &HepMC3::GenParticle::set_pdg_id, "Set PDG ID\n \n\n Use set_pid() instead\n\nC++: HepMC3::GenParticle::set_pdg_id(const int &) --> void", pybind11::arg("pidin"));
		cl.def("assign", (class HepMC3::GenParticle & (HepMC3::GenParticle::*)(const class HepMC3::GenParticle &)) &HepMC3::GenParticle::operator=, "C++: HepMC3::GenParticle::operator=(const class HepMC3::GenParticle &) --> class HepMC3::GenParticle &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_GenParticle_binder(cl);
	}
	{ // HepMC3::GenVertexData file: line:22
		pybind11::class_<HepMC3::GenVertexData, std::shared_ptr<HepMC3::GenVertexData>> cl(M("HepMC3"), "GenVertexData", "");
		cl.def( pybind11::init( [](){ return new HepMC3::GenVertexData(); } ) );
		cl.def( pybind11::init( [](HepMC3::GenVertexData const &o){ return new HepMC3::GenVertexData(o); } ) );
		cl.def_readwrite("status", &HepMC3::GenVertexData::status);
		cl.def_readwrite("position", &HepMC3::GenVertexData::position);
		cl.def("is_zero", (bool (HepMC3::GenVertexData::*)() const) &HepMC3::GenVertexData::is_zero, "Check if this struct fields are zero\n\nC++: HepMC3::GenVertexData::is_zero() const --> bool");
		cl.def("assign", (struct HepMC3::GenVertexData & (HepMC3::GenVertexData::*)(const struct HepMC3::GenVertexData &)) &HepMC3::GenVertexData::operator=, "C++: HepMC3::GenVertexData::operator=(const struct HepMC3::GenVertexData &) --> struct HepMC3::GenVertexData &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::GenVertex file:HepMC3/GenVertex.h line:26
		pybind11::class_<HepMC3::GenVertex, std::shared_ptr<HepMC3::GenVertex>> cl(M("HepMC3"), "GenVertex", "Stores vertex-related information");
		cl.def( pybind11::init( [](){ return new HepMC3::GenVertex(); } ), "doc" );
		cl.def( pybind11::init<const class HepMC3::FourVector &>(), pybind11::arg("position") );

		cl.def( pybind11::init<const struct HepMC3::GenVertexData &>(), pybind11::arg("data") );

		cl.def( pybind11::init( [](HepMC3::GenVertex const &o){ return new HepMC3::GenVertex(o); } ) );
		cl.def("parent_event", (class HepMC3::GenEvent * (HepMC3::GenVertex::*)()) &HepMC3::GenVertex::parent_event, "Get parent event\n\nC++: HepMC3::GenVertex::parent_event() --> class HepMC3::GenEvent *", pybind11::return_value_policy::automatic);
		cl.def("in_event", (bool (HepMC3::GenVertex::*)() const) &HepMC3::GenVertex::in_event, "Check if this vertex belongs to an event\n\nC++: HepMC3::GenVertex::in_event() const --> bool");
		cl.def("id", (int (HepMC3::GenVertex::*)() const) &HepMC3::GenVertex::id, "Get the vertex unique identifier\n\n \n This is not the same as id() in HepMC v2, which is now \n\nC++: HepMC3::GenVertex::id() const --> int");
		cl.def("set_id", (void (HepMC3::GenVertex::*)(int)) &HepMC3::GenVertex::set_id, "set the vertex identifier\n\nC++: HepMC3::GenVertex::set_id(int) --> void", pybind11::arg("id"));
		cl.def("status", (int (HepMC3::GenVertex::*)() const) &HepMC3::GenVertex::status, "Get vertex status code\n\nC++: HepMC3::GenVertex::status() const --> int");
		cl.def("set_status", (void (HepMC3::GenVertex::*)(int)) &HepMC3::GenVertex::set_status, "Set vertex status code\n\nC++: HepMC3::GenVertex::set_status(int) --> void", pybind11::arg("stat"));
		cl.def("data", (const struct HepMC3::GenVertexData & (HepMC3::GenVertex::*)() const) &HepMC3::GenVertex::data, "Get vertex data\n\nC++: HepMC3::GenVertex::data() const --> const struct HepMC3::GenVertexData &", pybind11::return_value_policy::automatic);
		cl.def("add_particle_in", (void (HepMC3::GenVertex::*)(class std::shared_ptr<class HepMC3::GenParticle>)) &HepMC3::GenVertex::add_particle_in, "Add incoming particle\n\nC++: HepMC3::GenVertex::add_particle_in(class std::shared_ptr<class HepMC3::GenParticle>) --> void", pybind11::arg("p"));
		cl.def("add_particle_out", (void (HepMC3::GenVertex::*)(class std::shared_ptr<class HepMC3::GenParticle>)) &HepMC3::GenVertex::add_particle_out, "Add outgoing particle\n\nC++: HepMC3::GenVertex::add_particle_out(class std::shared_ptr<class HepMC3::GenParticle>) --> void", pybind11::arg("p"));
		cl.def("remove_particle_in", (void (HepMC3::GenVertex::*)(class std::shared_ptr<class HepMC3::GenParticle>)) &HepMC3::GenVertex::remove_particle_in, "Remove incoming particle\n\nC++: HepMC3::GenVertex::remove_particle_in(class std::shared_ptr<class HepMC3::GenParticle>) --> void", pybind11::arg("p"));
		cl.def("remove_particle_out", (void (HepMC3::GenVertex::*)(class std::shared_ptr<class HepMC3::GenParticle>)) &HepMC3::GenVertex::remove_particle_out, "Remove outgoing particle\n\nC++: HepMC3::GenVertex::remove_particle_out(class std::shared_ptr<class HepMC3::GenParticle>) --> void", pybind11::arg("p"));
		cl.def("particles_in_size", (int (HepMC3::GenVertex::*)() const) &HepMC3::GenVertex::particles_in_size, "Number of incoming particles, HepMC2 compatiility\n\nC++: HepMC3::GenVertex::particles_in_size() const --> int");
		cl.def("particles_out_size", (int (HepMC3::GenVertex::*)() const) &HepMC3::GenVertex::particles_out_size, "Number of outgoing particles, HepMC2 compatiility\n\nC++: HepMC3::GenVertex::particles_out_size() const --> int");
		cl.def("particles_in", (const class std::vector<class std::shared_ptr<class HepMC3::GenParticle> > & (HepMC3::GenVertex::*)()) &HepMC3::GenVertex::particles_in, "Get list of incoming particles\n\nC++: HepMC3::GenVertex::particles_in() --> const class std::vector<class std::shared_ptr<class HepMC3::GenParticle> > &", pybind11::return_value_policy::automatic);
		cl.def("particles_out", (const class std::vector<class std::shared_ptr<class HepMC3::GenParticle> > & (HepMC3::GenVertex::*)()) &HepMC3::GenVertex::particles_out, "Get list of outgoing particles\n\nC++: HepMC3::GenVertex::particles_out() --> const class std::vector<class std::shared_ptr<class HepMC3::GenParticle> > &", pybind11::return_value_policy::automatic);
		cl.def("position", (const class HepMC3::FourVector & (HepMC3::GenVertex::*)() const) &HepMC3::GenVertex::position, "Get vertex position\n\n Returns the position of this vertex. If a position is not set on _this_ vertex,\n the production vertices of ancestors are searched to find the inherited position.\n FourVector(0,0,0,0) is returned if no position information is found.\n\nC++: HepMC3::GenVertex::position() const --> const class HepMC3::FourVector &", pybind11::return_value_policy::automatic);
		cl.def("has_set_position", (bool (HepMC3::GenVertex::*)() const) &HepMC3::GenVertex::has_set_position, "Check if position of this vertex is set\n\nC++: HepMC3::GenVertex::has_set_position() const --> bool");
		cl.def("set_position", (void (HepMC3::GenVertex::*)(const class HepMC3::FourVector &)) &HepMC3::GenVertex::set_position, "Set vertex position\n\nC++: HepMC3::GenVertex::set_position(const class HepMC3::FourVector &) --> void", pybind11::arg("new_pos"));
		cl.def("add_attribute", (bool (HepMC3::GenVertex::*)(const std::string &, class std::shared_ptr<class HepMC3::Attribute>)) &HepMC3::GenVertex::add_attribute, "Add event attribute to this vertex\n\n This will overwrite existing attribute if an attribute with\n the same name is present. The attribute will be stored in the\n parent_event(). \n\n false if there is no parent_event();\n\nC++: HepMC3::GenVertex::add_attribute(const std::string &, class std::shared_ptr<class HepMC3::Attribute>) --> bool", pybind11::arg("name"), pybind11::arg("att"));
		cl.def("attribute_names", (class std::vector<std::string > (HepMC3::GenVertex::*)() const) &HepMC3::GenVertex::attribute_names, "Get list of names of attributes assigned to this particle\n\nC++: HepMC3::GenVertex::attribute_names() const --> class std::vector<std::string >");
		cl.def("remove_attribute", (void (HepMC3::GenVertex::*)(const std::string &)) &HepMC3::GenVertex::remove_attribute, "Remove attribute\n\nC++: HepMC3::GenVertex::remove_attribute(const std::string &) --> void", pybind11::arg("name"));
		cl.def("attribute_as_string", (std::string (HepMC3::GenVertex::*)(const std::string &) const) &HepMC3::GenVertex::attribute_as_string, "Get attribute of any type as string\n\nC++: HepMC3::GenVertex::attribute_as_string(const std::string &) const --> std::string", pybind11::arg("name"));
		cl.def("add_particle_in", (void (HepMC3::GenVertex::*)(class HepMC3::GenParticle *)) &HepMC3::GenVertex::add_particle_in, "Add incoming particle by raw pointer\n \n\n Use GenVertex::add_particle_in( const GenParticlePtr &p ) instead\n\nC++: HepMC3::GenVertex::add_particle_in(class HepMC3::GenParticle *) --> void", pybind11::arg("p"));
		cl.def("add_particle_out", (void (HepMC3::GenVertex::*)(class HepMC3::GenParticle *)) &HepMC3::GenVertex::add_particle_out, "Add outgoing particle by raw pointer\n \n\n Use GenVertex::add_particle_out( const GenParticlePtr &p ) instead\n\nC++: HepMC3::GenVertex::add_particle_out(class HepMC3::GenParticle *) --> void", pybind11::arg("p"));
		cl.def("assign", (class HepMC3::GenVertex & (HepMC3::GenVertex::*)(const class HepMC3::GenVertex &)) &HepMC3::GenVertex::operator=, "C++: HepMC3::GenVertex::operator=(const class HepMC3::GenVertex &) --> class HepMC3::GenVertex &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_GenVertex_binder(cl);
	}
	{ // HepMC3::GenRunInfoData file:HepMC3/Data/GenRunInfoData.h line:23
		pybind11::class_<HepMC3::GenRunInfoData, std::shared_ptr<HepMC3::GenRunInfoData>> cl(M("HepMC3"), "GenRunInfoData", "");
		cl.def( pybind11::init( [](){ return new HepMC3::GenRunInfoData(); } ) );
		cl.def( pybind11::init( [](HepMC3::GenRunInfoData const &o){ return new HepMC3::GenRunInfoData(o); } ) );
		cl.def_readwrite("weight_names", &HepMC3::GenRunInfoData::weight_names);
		cl.def_readwrite("tool_name", &HepMC3::GenRunInfoData::tool_name);
		cl.def_readwrite("tool_version", &HepMC3::GenRunInfoData::tool_version);
		cl.def_readwrite("tool_description", &HepMC3::GenRunInfoData::tool_description);
		cl.def_readwrite("attribute_name", &HepMC3::GenRunInfoData::attribute_name);
		cl.def_readwrite("attribute_string", &HepMC3::GenRunInfoData::attribute_string);
		cl.def("assign", (struct HepMC3::GenRunInfoData & (HepMC3::GenRunInfoData::*)(const struct HepMC3::GenRunInfoData &)) &HepMC3::GenRunInfoData::operator=, "C++: HepMC3::GenRunInfoData::operator=(const struct HepMC3::GenRunInfoData &) --> struct HepMC3::GenRunInfoData &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // HepMC3::GenEventData file:HepMC3/Data/GenEventData.h line:26
		pybind11::class_<HepMC3::GenEventData, std::shared_ptr<HepMC3::GenEventData>> cl(M("HepMC3"), "GenEventData", "");
		cl.def( pybind11::init( [](){ return new HepMC3::GenEventData(); } ) );
		cl.def( pybind11::init( [](HepMC3::GenEventData const &o){ return new HepMC3::GenEventData(o); } ) );
		cl.def_readwrite("event_number", &HepMC3::GenEventData::event_number);
		cl.def_readwrite("momentum_unit", &HepMC3::GenEventData::momentum_unit);
		cl.def_readwrite("length_unit", &HepMC3::GenEventData::length_unit);
		cl.def_readwrite("particles", &HepMC3::GenEventData::particles);
		cl.def_readwrite("vertices", &HepMC3::GenEventData::vertices);
		cl.def_readwrite("weights", &HepMC3::GenEventData::weights);
		cl.def_readwrite("event_pos", &HepMC3::GenEventData::event_pos);
		cl.def_readwrite("links1", &HepMC3::GenEventData::links1);
		cl.def_readwrite("links2", &HepMC3::GenEventData::links2);
		cl.def_readwrite("attribute_id", &HepMC3::GenEventData::attribute_id);
		cl.def_readwrite("attribute_name", &HepMC3::GenEventData::attribute_name);
		cl.def_readwrite("attribute_string", &HepMC3::GenEventData::attribute_string);
		cl.def("assign", (struct HepMC3::GenEventData & (HepMC3::GenEventData::*)(const struct HepMC3::GenEventData &)) &HepMC3::GenEventData::operator=, "C++: HepMC3::GenEventData::operator=(const struct HepMC3::GenEventData &) --> struct HepMC3::GenEventData &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	// HepMC3::version() file: line:20
	M("HepMC3").def("version", (std::string (*)()) &HepMC3::version, "Get the HepMC library version string\n\nC++: HepMC3::version() --> std::string");

	// HepMC3::Print file: line:27
	binder::print_binder(M("HepMC3"));

	{ // HepMC3::Writer file:HepMC3/Writer.h line:25
		pybind11::class_<HepMC3::Writer, std::shared_ptr<HepMC3::Writer>, PyCallBack_HepMC3_Writer> cl(M("HepMC3"), "Writer", "");
		cl.def( pybind11::init( [](){ return new PyCallBack_HepMC3_Writer(); } ) );
		cl.def("write_event", (void (HepMC3::Writer::*)(const class HepMC3::GenEvent &)) &HepMC3::Writer::write_event, "Write event  to output target\n\nC++: HepMC3::Writer::write_event(const class HepMC3::GenEvent &) --> void", pybind11::arg("evt"));
		cl.def("failed", (bool (HepMC3::Writer::*)()) &HepMC3::Writer::failed, "Get file and/or stream error state \n\nC++: HepMC3::Writer::failed() --> bool");
		cl.def("close", (void (HepMC3::Writer::*)()) &HepMC3::Writer::close, "Close file and/or stream \n\nC++: HepMC3::Writer::close() --> void");
		cl.def("set_run_info", (void (HepMC3::Writer::*)(class std::shared_ptr<class HepMC3::GenRunInfo>)) &HepMC3::Writer::set_run_info, "Set the global GenRunInfo object.\n\nC++: HepMC3::Writer::set_run_info(class std::shared_ptr<class HepMC3::GenRunInfo>) --> void", pybind11::arg("run"));
		cl.def("run_info", (class std::shared_ptr<class HepMC3::GenRunInfo> (HepMC3::Writer::*)() const) &HepMC3::Writer::run_info, "Get the global GenRunInfo object.\n\nC++: HepMC3::Writer::run_info() const --> class std::shared_ptr<class HepMC3::GenRunInfo>");
		cl.def("set_options", (void (HepMC3::Writer::*)(const class std::map<std::string, std::string > &)) &HepMC3::Writer::set_options, "Set options\n\nC++: HepMC3::Writer::set_options(const class std::map<std::string, std::string > &) --> void", pybind11::arg("options"));
		cl.def("get_options", (class std::map<std::string, std::string > (HepMC3::Writer::*)() const) &HepMC3::Writer::get_options, "Set options\n\nC++: HepMC3::Writer::get_options() const --> class std::map<std::string, std::string >");
	}
}
