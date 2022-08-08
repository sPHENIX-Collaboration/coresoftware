#include "root_binders.h"

namespace binder {


	void	WriterRootTree_binder(pybind11::module &M)
	{
		pybind11::class_<HepMC3::WriterRootTree, std::shared_ptr<HepMC3::WriterRootTree>, HepMC3::Writer> cl(M, "WriterRootTree", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](const class std::basic_string<char> & a0){ return new HepMC3::WriterRootTree(a0); }/*, [](const class std::basic_string<char> & a0){ return new PyCallBack_HepMC3_WriterRootTree(a0); }*/ ), "doc");
		cl.def( pybind11::init<const std::string &, class std::shared_ptr<class HepMC3::GenRunInfo>>(), pybind11::arg("filename"), pybind11::arg("run") );

		cl.def( pybind11::init( [](const class std::basic_string<char> & a0, const class std::basic_string<char> & a1, const class std::basic_string<char> & a2){ return new HepMC3::WriterRootTree(a0, a1, a2); }/*, [](const class std::basic_string<char> & a0, const class std::basic_string<char> & a1, const class std::basic_string<char> & a2){ return new PyCallBack_HepMC3_WriterRootTree(a0, a1, a2); }*/ ), "doc");
		cl.def( pybind11::init<const std::string &, const std::string &, const std::string &, class std::shared_ptr<class HepMC3::GenRunInfo>>(), pybind11::arg("filename"), pybind11::arg("treename"), pybind11::arg("branchname"), pybind11::arg("run") );

		cl.def("write_event", (void (HepMC3::WriterRootTree::*)(const class HepMC3::GenEvent &)) &HepMC3::WriterRootTree::write_event, "Write event to file\n\n  \n Event to be serialized\n\nC++: HepMC3::WriterRootTree::write_event(const class HepMC3::GenEvent &) --> void", pybind11::arg("evt"));
		cl.def("write_run_info", (void (HepMC3::WriterRootTree::*)()) &HepMC3::WriterRootTree::write_run_info, "Write the GenRunInfo object to file. \n\nC++: HepMC3::WriterRootTree::write_run_info() --> void");
		cl.def("close", (void (HepMC3::WriterRootTree::*)()) &HepMC3::WriterRootTree::close, "Close file stream \n\nC++: HepMC3::WriterRootTree::close() --> void");
		cl.def("failed", (bool (HepMC3::WriterRootTree::*)()) &HepMC3::WriterRootTree::failed, "Get stream error state flag \n\nC++: HepMC3::WriterRootTree::failed() --> bool");
	}

	void	WriterRoot_binder(pybind11::module &M)
	{
			pybind11::class_<HepMC3::WriterRoot, std::shared_ptr<HepMC3::WriterRoot>,  HepMC3::Writer> cl(M, "WriterRoot", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](const class std::basic_string<char> & a0){ return new HepMC3::WriterRoot(a0); }/*, [](const class std::basic_string<char> & a0){ return new PyCallBack_HepMC3_WriterRoot(a0); } */), "doc");
		cl.def( pybind11::init<const std::string &, class std::shared_ptr<class HepMC3::GenRunInfo>>(), pybind11::arg("filename"), pybind11::arg("run") );

		cl.def("write_event", (void (HepMC3::WriterRoot::*)(const class HepMC3::GenEvent &)) &HepMC3::WriterRoot::write_event, "Write event to file\n\n  \n Event to be serialized\n\nC++: HepMC3::WriterRoot::write_event(const class HepMC3::GenEvent &) --> void", pybind11::arg("evt"));
		cl.def("write_run_info", (void (HepMC3::WriterRoot::*)()) &HepMC3::WriterRoot::write_run_info, "Write the GenRunInfo object to file. \n\nC++: HepMC3::WriterRoot::write_run_info() --> void");
		cl.def("close", (void (HepMC3::WriterRoot::*)()) &HepMC3::WriterRoot::close, "Close file stream \n\nC++: HepMC3::WriterRoot::close() --> void");
		cl.def("failed", (bool (HepMC3::WriterRoot::*)()) &HepMC3::WriterRoot::failed, "Get stream error state flag \n\nC++: HepMC3::WriterRoot::failed() --> bool");
	}
	void	ReaderRootTree_binder(pybind11::module &M)
	{
			pybind11::class_<HepMC3::ReaderRootTree, std::shared_ptr<HepMC3::ReaderRootTree>,  HepMC3::Reader> cl(M, "ReaderRootTree", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<const std::string &>(), pybind11::arg("filename") );

		cl.def( pybind11::init<const std::string &, const std::string &, const std::string &>(), pybind11::arg("filename"), pybind11::arg("treename"), pybind11::arg("branchname") );

		cl.def("skip", (bool (HepMC3::ReaderRootTree::*)(const int)) &HepMC3::ReaderRootTree::skip, "skip events\n\nC++: HepMC3::ReaderRootTree::skip(const int) --> bool", pybind11::arg(""));

		cl.def("read_event", (bool (HepMC3::ReaderRootTree::*)(class HepMC3::GenEvent &)) &HepMC3::ReaderRootTree::read_event, "Read event from file\n\n  \n Contains parsed event\n\nC++: HepMC3::ReaderRootTree::read_event(class HepMC3::GenEvent &) --> bool", pybind11::arg("evt"));
		cl.def("close", (void (HepMC3::ReaderRootTree::*)()) &HepMC3::ReaderRootTree::close, "Close file \n\nC++: HepMC3::ReaderRootTree::close() --> void");
		cl.def("failed", (bool (HepMC3::ReaderRootTree::*)()) &HepMC3::ReaderRootTree::failed, "Get file  error state \n\nC++: HepMC3::ReaderRootTree::failed() --> bool");

}

	void	ReaderRoot_binder(pybind11::module &M)
	{
			pybind11::class_<HepMC3::ReaderRoot, std::shared_ptr<HepMC3::ReaderRoot>,  HepMC3::Reader> cl(M, "ReaderRoot", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<const std::string &>(), pybind11::arg("filename") );

		cl.def("skip", (bool (HepMC3::ReaderRoot::*)(const int)) &HepMC3::ReaderRoot::skip, "skip events\n\nC++: HepMC3::ReaderRoot::skip(const int) --> bool", pybind11::arg(""));

		cl.def("read_event", (bool (HepMC3::ReaderRoot::*)(class HepMC3::GenEvent &)) &HepMC3::ReaderRoot::read_event, "Read event from file\n\n  \n Contains parsed event\n\nC++: HepMC3::ReaderRoot::read_event(class HepMC3::GenEvent &) --> bool", pybind11::arg("evt"));
		cl.def("close", (void (HepMC3::ReaderRoot::*)()) &HepMC3::ReaderRoot::close, "Close file stream \n\nC++: HepMC3::ReaderRoot::close() --> void");
		cl.def("failed", (bool (HepMC3::ReaderRoot::*)()) &HepMC3::ReaderRoot::failed, "Get stream error state \n\nC++: HepMC3::ReaderRoot::failed() --> bool");

}

} // namespace binder
