#include <HepMC3/LHEF.h>
#include <ios>
#include <istream>
#include <iterator>
#include <memory>
#include <ostream>
#include <sstream> // __str__
#include <streambuf>
#include <string>

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

void bind_pyHepMC3_16(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // LHEF::Reader file:HepMC3/LHEF.h line:2742
		pybind11::class_<LHEF::Reader, std::shared_ptr<LHEF::Reader>> cl(M("LHEF"), "Reader", "The Reader class is initialized with a stream from which to read a\n version 1/2 Les Houches Accord event file. In the constructor of\n the Reader object the optional header information is read and then\n the mandatory init is read. After this the whole header block\n including the enclosing lines with tags are available in the public\n headerBlock member variable. Also the information from the init\n block is available in the heprup member variable and any additional\n comment lines are available in initComments. After each successful\n call to the readEvent() function the standard Les Houches Accord\n information about the event is available in the hepeup member\n variable and any additional comments in the eventComments\n variable. A typical reading sequence would look as follows:\n\n ");
		cl.def( pybind11::init<std::string>(), pybind11::arg("filename") );

		cl.def_readwrite("version", &LHEF::Reader::version);
		cl.def_readwrite("outsideBlock", &LHEF::Reader::outsideBlock);
		cl.def_readwrite("headerBlock", &LHEF::Reader::headerBlock);
		cl.def_readwrite("heprup", &LHEF::Reader::heprup);
		cl.def_readwrite("initComments", &LHEF::Reader::initComments);
		cl.def_readwrite("hepeup", &LHEF::Reader::hepeup);
		cl.def_readwrite("eventComments", &LHEF::Reader::eventComments);
		cl.def_readwrite("currevent", &LHEF::Reader::currevent);
		cl.def_readwrite("curreventfile", &LHEF::Reader::curreventfile);
		cl.def_readwrite("currfileevent", &LHEF::Reader::currfileevent);
		cl.def_readwrite("dirpath", &LHEF::Reader::dirpath);
		cl.def("readEvent", (bool (LHEF::Reader::*)()) &LHEF::Reader::readEvent, "Read an event from the file and store it in the hepeup\n object. Optional comment lines are stored i the eventComments\n member variable.\n \n\n true if the read sas successful.\n\nC++: LHEF::Reader::readEvent() --> bool");
		cl.def("openeventfile", (void (LHEF::Reader::*)(int)) &LHEF::Reader::openeventfile, "Open the efentfile with index ifile. If another eventfile is\n being read, its remaining contents is discarded. This is a noop\n if current read session is not a multi-file run.\n\nC++: LHEF::Reader::openeventfile(int) --> void", pybind11::arg("ifile"));
	}
	{ // LHEF::Writer file:HepMC3/LHEF.h line:3090
		pybind11::class_<LHEF::Writer, std::shared_ptr<LHEF::Writer>> cl(M("LHEF"), "Writer", "The Writer class is initialized with a stream to which to write a\n version 1.0 Les Houches Accord event file. In the constructor of\n the Writer object the main XML tag is written out, with the\n corresponding end tag is written in the destructor. After a Writer\n object has been created, it is possible to assign standard init\n information in the heprup member variable. In addition any XML\n formatted information can be added to the headerBlock member\n variable (directly or via the addHeader() function). Further\n comment line (beginning with a # character) can be\n added to the initComments variable (directly or with the\n addInitComment() function). After this information is set, it\n should be written out to the file with the init() function.\n\n Before each event is written out with the writeEvent() function,\n the standard event information can then be assigned to the hepeup\n variable and optional comment lines (beginning with a\n # character) may be given to the eventComments\n variable (directly or with the addEventComment() function).\n\n ");
		cl.def( pybind11::init<std::string>(), pybind11::arg("filename") );

		cl.def_readwrite("heprup", &LHEF::Writer::heprup);
		cl.def_readwrite("hepeup", &LHEF::Writer::hepeup);
		cl.def("headerBlock", (void (LHEF::Writer::*)(const std::string &)) &LHEF::Writer::headerBlock, "Add header lines consisting of XML code with this stream.\n\nC++: LHEF::Writer::headerBlock(const std::string &) --> void", pybind11::arg("a"));
		cl.def("initComments", (void (LHEF::Writer::*)(const std::string &)) &LHEF::Writer::initComments, "Add comment lines to the init block with this stream.\n\nC++: LHEF::Writer::initComments(const std::string &) --> void", pybind11::arg("a"));
		cl.def("eventComments", (void (LHEF::Writer::*)(const std::string &)) &LHEF::Writer::eventComments, "Add comment lines to the next event to be written out with this stream.\n\nC++: LHEF::Writer::eventComments(const std::string &) --> void", pybind11::arg("a"));
		cl.def("init", (void (LHEF::Writer::*)()) &LHEF::Writer::init, "Initialize the writer.\n\nC++: LHEF::Writer::init() --> void");
		cl.def("openeventfile", (bool (LHEF::Writer::*)(int)) &LHEF::Writer::openeventfile, "Open a new event file, possibly closing a previous opened one.\n\nC++: LHEF::Writer::openeventfile(int) --> bool", pybind11::arg("ifile"));
		cl.def("writeinit", (void (LHEF::Writer::*)()) &LHEF::Writer::writeinit, "Write out an optional header block followed by the standard init\n block information together with any comment lines.\n\nC++: LHEF::Writer::writeinit() --> void");
		cl.def("writeEvent", (void (LHEF::Writer::*)()) &LHEF::Writer::writeEvent, "Write the current HEPEUP object to the stream;\n\nC++: LHEF::Writer::writeEvent() --> void");
	}
}
