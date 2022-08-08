#include <map>
#include <memory>
#include <stdexcept>
#include <functional>
#include <string>

#include <pybind11/pybind11.h>

typedef std::function< pybind11::module & (std::string const &) > ModuleGetter;

void bind_pyHepMC3_0(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_1(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_2(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_3(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_4(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_5(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_6(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_7(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_8(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_9(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_10(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_11(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_12(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_13(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_14(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_15(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_16(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_17(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_pyHepMC3_18(std::function< pybind11::module &(std::string const &namespace_) > &M);


PYBIND11_MODULE(pyHepMC3, root_module) {
	root_module.doc() = "pyHepMC3 module";

	std::map <std::string, pybind11::module> modules;
	ModuleGetter M = [&](std::string const &namespace_) -> pybind11::module & {
		auto it = modules.find(namespace_);
		if( it == modules.end() ) throw std::runtime_error("Attempt to access pybind11::module for namespace " + namespace_ + " before it was created!!!");
		return it->second;
	};

	modules[""] = root_module;

	std::vector< std::pair<std::string, std::string> > sub_modules {
		{"", "HepMC3"},
		{"", "LHEF"},
		{"", "std"},
	};
	for(auto &p : sub_modules ) modules[p.first.size() ? p.first+"::"+p.second : p.second] = modules[p.first].def_submodule(p.second.c_str(), ("Bindings for " + p.first + "::" + p.second + " namespace").c_str() );

	//pybind11::class_<std::shared_ptr<void>>(M(""), "_encapsulated_data_");

	bind_pyHepMC3_0(M);
	bind_pyHepMC3_1(M);
	bind_pyHepMC3_2(M);
	bind_pyHepMC3_3(M);
	bind_pyHepMC3_4(M);
	bind_pyHepMC3_5(M);
	bind_pyHepMC3_6(M);
	bind_pyHepMC3_7(M);
	bind_pyHepMC3_8(M);
	bind_pyHepMC3_9(M);
	bind_pyHepMC3_10(M);
	bind_pyHepMC3_11(M);
	bind_pyHepMC3_12(M);
	bind_pyHepMC3_13(M);
	bind_pyHepMC3_14(M);
	bind_pyHepMC3_15(M);
	bind_pyHepMC3_16(M);
	bind_pyHepMC3_17(M);
	bind_pyHepMC3_18(M);

}
