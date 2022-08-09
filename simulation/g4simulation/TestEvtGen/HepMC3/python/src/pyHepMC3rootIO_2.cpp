#include <HepMC3/ReaderRoot.h>
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
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/WriterRoot.h>
#include <HepMC3/WriterRootTree.h>
#include <HepMC3/ReaderRoot.h>
#include <HepMC3/ReaderRootTree.h>
#include <src/root_binders.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

void bind_pyHepMC3rootIO_2(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// HepMC3::ReaderRoot file:HepMC3/ReaderRoot.h line:33
	 binder::ReaderRoot_binder(M("HepMC3"));

}
