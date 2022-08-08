#ifndef ROOT_BINDERS_H
#define ROOT_BINDERS_H
#include <HepMC3/WriterRoot.h>
#include <HepMC3/WriterRootTree.h>
#include <HepMC3/ReaderRoot.h>
#include <HepMC3/ReaderRootTree.h>
#include <pybind11/pybind11.h>
namespace binder {
	void	WriterRootTree_binder(pybind11::module &M);
	void	WriterRoot_binder(pybind11::module &M);
	void	ReaderRootTree_binder(pybind11::module &M);
	void	ReaderRoot_binder(pybind11::module &M);
} // namespace binder
#endif // _INCLUDED_root_binders_hpp_
