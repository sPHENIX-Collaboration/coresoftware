#include "protobuf_binders.h"

namespace binder {



	void	Writerprotobuf_binder(pybind11::module &M)
	{
			pybind11::class_<HepMC3::Writerprotobuf, std::shared_ptr<HepMC3::Writerprotobuf>,  HepMC3::Writer> cl(M, "Writerprotobuf", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](const class std::basic_string<char> & a0){ return new HepMC3::Writerprotobuf(a0); }/*, [](const class std::basic_string<char> & a0){ return new PyCallBack_HepMC3_Writerprotobuf(a0); } */), "doc");
		cl.def( pybind11::init<const std::string &, class std::shared_ptr<class HepMC3::GenRunInfo>>(), pybind11::arg("filename"), pybind11::arg("run") );

		cl.def("write_event", (void (HepMC3::Writerprotobuf::*)(const class HepMC3::GenEvent &)) &HepMC3::Writerprotobuf::write_event, "Write event to file\n\n  \n Event to be serialized\n\nC++: HepMC3::Writerprotobuf::write_event(const class HepMC3::GenEvent &) --> void", pybind11::arg("evt"));
		cl.def("write_run_info", (void (HepMC3::Writerprotobuf::*)()) &HepMC3::Writerprotobuf::write_run_info, "Write the GenRunInfo object to file. \n\nC++: HepMC3::Writerprotobuf::write_run_info() --> void");
		cl.def("close", (void (HepMC3::Writerprotobuf::*)()) &HepMC3::Writerprotobuf::close, "Close file stream \n\nC++: HepMC3::Writerprotobuf::close() --> void");
		cl.def("failed", (bool (HepMC3::Writerprotobuf::*)()) &HepMC3::Writerprotobuf::failed, "Get stream error state flag \n\nC++: HepMC3::Writerprotobuf::failed() --> bool");
	}


	void	Readerprotobuf_binder(pybind11::module &M)
	{
			pybind11::class_<HepMC3::Readerprotobuf, std::shared_ptr<HepMC3::Readerprotobuf>,  HepMC3::Reader> cl(M, "Readerprotobuf", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<const std::string &>(), pybind11::arg("filename") );

		cl.def("skip", (bool (HepMC3::Readerprotobuf::*)(const int)) &HepMC3::Readerprotobuf::skip, "skip events\n\nC++: HepMC3::Readerprotobuf::skip(const int) --> bool", pybind11::arg(""));

		cl.def("read_event", (bool (HepMC3::Readerprotobuf::*)(class HepMC3::GenEvent &)) &HepMC3::Readerprotobuf::read_event, "Read event from file\n\n  \n Contains parsed event\n\nC++: HepMC3::Readerprotobuf::read_event(class HepMC3::GenEvent &) --> bool", pybind11::arg("evt"));
		cl.def("close", (void (HepMC3::Readerprotobuf::*)()) &HepMC3::Readerprotobuf::close, "Close file stream \n\nC++: HepMC3::Readerprotobuf::close() --> void");
		cl.def("failed", (bool (HepMC3::Readerprotobuf::*)()) &HepMC3::Readerprotobuf::failed, "Get stream error state \n\nC++: HepMC3::Readerprotobuf::failed() --> bool");

}

} // namespace binder
