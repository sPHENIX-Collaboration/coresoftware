#ifndef PROTOBUF_BINDERS_H
#define PROTOBUF_BINDERS_H
#include <HepMC3/Writerprotobuf.h>
#include <HepMC3/Readerprotobuf.h>
#include <pybind11/pybind11.h>
namespace binder {
	void	Writerprotobuf_binder(pybind11::module &M);
	void	Readerprotobuf_binder(pybind11::module &M);
} // namespace binder
#endif // _INCLUDED_protobuf_binders_hpp_
