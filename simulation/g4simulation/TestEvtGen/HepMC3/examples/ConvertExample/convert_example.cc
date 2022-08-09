// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
/// @example convert_example.cc
/// @brief Utility to convert between different types of event records
///
#include "HepMC3/Print.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/Reader.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/WriterHEPEVT.h"
#include "HepMC3/WriterPlugin.h"
#include "HepMC3/ReaderHEPEVT.h"
#include "HepMC3/ReaderLHEF.h"
#include "HepMC3/ReaderPlugin.h"
#include "HepMC3/ReaderFactory.h"

#ifdef HEPMC3_ROOTIO
#include "HepMC3/ReaderRoot.h"
#include "HepMC3/WriterRoot.h"
#include "HepMC3/ReaderRootTree.h"
#include "HepMC3/WriterRootTree.h"
#endif
#if HEPMC3_USE_COMPRESSION
#include "HepMC3/ReaderGZ.h"
#include "HepMC3/WriterGZ.h"
#endif
#ifdef HEPMC3_PROTOBUFIO
#include "HepMC3/Readerprotobuf.h"
#include "HepMC3/Writerprotobuf.h"
#endif

/* Extension example*/
#ifdef HEPMCCONVERT_EXTENSION_ROOTTREEOPAL
#ifndef HEPMC3_ROOTIO
#warning "HEPMCCONVERT_EXTENSION_ROOTTREEOPAL requires  compilation with of HepMC with ROOT, i.e. HEPMC3_ROOTIO.This extension will be disabled."
#undef HEPMCCONVERT_EXTENSION_ROOTTREEOPAL
#else
#include "WriterRootTreeOPAL.h"
#endif
#endif
#ifdef HEPMCCONVERT_EXTENSION_HEPEVTZEUS
#include "WriterHEPEVTZEUS.h"
#endif
#ifdef HEPMCCONVERT_EXTENSION_DOT
#include "WriterDOT.h"
#endif
#ifdef HEPMCCONVERT_EXTENSION_UPROOTTREEREADER
#include "ReaderuprootTree.h"
#endif


#include "cmdline.h"
using namespace HepMC3;
enum formats {autodetect, hepmc2, hepmc3, hpe,root, treeroot, treerootopal, hpezeus, lhef, dump, dot, uproot, plugin, none, proto};

template <class T>
std::shared_ptr<Reader> get_input_file(const char* name, const bool input_is_stdin, const bool use_compression) {
    std::string n(name);
#if  HEPMC3_USE_COMPRESSION
    if (use_compression) {
        return (input_is_stdin?std::make_shared<ReaderGZ<T> >(std::cin):std::make_shared<ReaderGZ<T> >(n));
    }
#endif
    return (input_is_stdin?std::make_shared<T>(std::cin):std::make_shared<T>(n));
}
template <class T>
std::shared_ptr<Writer> get_output_file(const char* name, const char* use_compression) {
    std::string n(name);
#if HEPMC3_USE_COMPRESSION
    if (std::string(use_compression) == "z" )  return std::make_shared< WriterGZ<T,Compression::z> >(n);
    if (std::string(use_compression) == "lzma" )  return std::make_shared< WriterGZ<T,Compression::lzma> >(n);
    if (std::string(use_compression) == "bz2" )  return std::make_shared< WriterGZ<T,Compression::bz2> >(n);
#endif
    return std::make_shared<T>(n);
}

int main(int argc, char** argv)
{
    gengetopt_args_info ai;
    if (cmdline_parser (argc, argv, &ai) != 0) {
        exit(1);
    }
    if ( !( ( ai.inputs_num == 2 && ( std::string(ai.output_format_arg) !=  "none" )) || ( ai.inputs_num == 1 && ( std::string(ai.output_format_arg) ==  "none") ) )   )
    {
        printf("Exactly two arguments are requred: the name of input and output files if the output format in not \"none\"\n");
        printf("In case the output format is \"none\" exactly one argument should be given: the name of input file.\n");
        exit(1);
    }
    std::map<std::string,formats> format_map;
    format_map.insert(std::pair<std::string,formats> ( "auto", autodetect ));
    format_map.insert(std::pair<std::string,formats> ( "hepmc2", hepmc2 ));
    format_map.insert(std::pair<std::string,formats> ( "hepmc3", hepmc3 ));
    format_map.insert(std::pair<std::string,formats> ( "hpe", hpe  ));
    format_map.insert(std::pair<std::string,formats> ( "root", root ));
    format_map.insert(std::pair<std::string,formats> ( "treeroot", treeroot ));
    format_map.insert(std::pair<std::string,formats> ( "treerootopal", treerootopal ));
    format_map.insert(std::pair<std::string,formats> ( "hpezeus", hpezeus ));
    format_map.insert(std::pair<std::string,formats> ( "lhef", lhef ));
    format_map.insert(std::pair<std::string,formats> ( "dump", dump ));
    format_map.insert(std::pair<std::string,formats> ( "dot", dot ));
    format_map.insert(std::pair<std::string,formats> ( "uproot", uproot ));
    format_map.insert(std::pair<std::string,formats> ( "plugin", plugin ));
    format_map.insert(std::pair<std::string,formats> ( "none", none ));
    format_map.insert(std::pair<std::string,formats> ( "proto", proto ));
    std::map<std::string, std::string> options;
    for (size_t i=0; i<ai.extensions_given; i++)
    {
        std::string optarg=std::string(ai.extensions_arg[i]);
        size_t pos = optarg.find_first_of('=');
        if ( pos < optarg.length() )
            options[std::string(optarg,0,pos)] = std::string(optarg, pos+1, optarg.length());
    }
    long int  events_parsed = 0;
    long int  events_limit = ai.events_limit_arg;
    long int  first_event_number = ai.first_event_number_arg;
    long int  last_event_number = ai.last_event_number_arg;
    long int  print_each_events_parsed = ai.print_every_events_parsed_arg;
    std::string InputPluginLibrary;
    std::string InputPluginName;

    std::string OutputPluginLibrary;
    std::string OutputPluginName;

    std::shared_ptr<Reader>      input_file;
    bool input_is_stdin = (std::string(ai.inputs[0]) == std::string("-"));
    if (input_is_stdin) std::ios_base::sync_with_stdio(false);
    bool ignore_writer = false;
    switch (format_map.at(std::string(ai.input_format_arg)))
    {
    case autodetect:
        input_file = (input_is_stdin?deduce_reader(std::cin):deduce_reader(ai.inputs[0]));
        if (!input_file)
        {
            input_is_stdin?printf("Input format  detection for std input has failed\n"):printf("Input format  detection for file %s has failed\n",ai.inputs[0]);
            exit(2);
        }
        break;
    case hepmc2:
        input_file = get_input_file<ReaderAsciiHepMC2>(ai.inputs[0], input_is_stdin, ai.compressed_input_flag);
        break;
    case hepmc3:
        input_file = get_input_file<ReaderAscii>(ai.inputs[0], input_is_stdin, ai.compressed_input_flag);
        break;
    case hpe:
        input_file = get_input_file<ReaderHEPEVT>(ai.inputs[0], input_is_stdin,ai.compressed_input_flag);
        break;
    case lhef:
        input_file = get_input_file<ReaderLHEF>(ai.inputs[0], input_is_stdin, ai.compressed_input_flag);
        break;
    case uproot:
#ifdef HEPMCCONVERT_EXTENSION_UPROOTTREEREADER
        input_file = std::make_shared<ReaderuprootTree>(ai.inputs[0]);
        break;
#else
        printf("Input format %s  is not supported\n", ai.input_format_arg);
        exit(2);
#endif
    case treeroot:
#ifdef HEPMC3_ROOTIO
        input_file = std::make_shared<ReaderRootTree>(ai.inputs[0]);
        break;
#else
        printf("Input format %s  is not supported\n", ai.input_format_arg);
        exit(2);
#endif
    case root:
#ifdef HEPMC3_ROOTIO
        input_file = std::make_shared<ReaderRoot>(ai.inputs[0]);
        break;
#else
        printf("Input format %s  is not supported\n", ai.input_format_arg);
        exit(2);
#endif
    case proto:
#ifdef HEPMC3_PROTOBUFIO
        input_file = std::make_shared<Readerprotobuf>(ai.inputs[0]);
        break;
#else
        printf("Input format %s  is not supported\n", ai.input_format_arg);
        exit(2);
#endif
    case plugin:
        if (options.find("InputPluginLibrary") == options.end())         {
            printf("InputPluginLibrary option required\n");
            exit(2);
        }
        else InputPluginLibrary = options.at("InputPluginLibrary");
        if (options.find("InputPluginName") == options.end())            {
            printf("InputPluginName option required\n");
            exit(2);
        }
        else InputPluginName = options.at("InputPluginName");
        input_file = std::make_shared<ReaderPlugin>(std::string(ai.inputs[0]), InputPluginLibrary, InputPluginName);
        if (input_file->failed()) {
            printf("Plugin initialization failed\n");
            exit(2);
        }
        break;
    default:
        printf("Input format %s  is not known\n", ai.input_format_arg);
        exit(2);
        break;
    }
    std::shared_ptr<Writer>      output_file;
    switch (format_map.at(std::string(ai.output_format_arg)))
    {
    case hepmc2:
        output_file = get_output_file<WriterAsciiHepMC2>(ai.inputs[1], ai.compressed_output_arg);
        break;
    case hepmc3:
        output_file = get_output_file<WriterAscii>(ai.inputs[1], ai.compressed_output_arg);
        break;
    case hpe:
        output_file = get_output_file<WriterHEPEVT>(ai.inputs[1], ai.compressed_output_arg);
        break;
    case root:
#ifdef HEPMC3_ROOTIO
        output_file = std::make_shared<WriterRoot>(ai.inputs[1]);
        break;
#else
        printf("Output format %s  is not supported\n", ai.output_format_arg);
        exit(2);
#endif
    case proto:
#ifdef HEPMC3_PROTOBUFIO
        output_file = std::make_shared<Writerprotobuf>(ai.inputs[1]);
        break;
#else
        printf("Output format %s  is not supported\n", ai.output_format_arg);
        exit(2);
#endif
    case treeroot:
#ifdef HEPMC3_ROOTIO
        output_file = std::make_shared<WriterRootTree>(ai.inputs[1]);
        break;
#else
        printf("Output format %s  is not supported\n",ai.output_format_arg);
        exit(2);
#endif
    /* Extension example*/
    case treerootopal:
#ifdef HEPMCCONVERT_EXTENSION_ROOTTREEOPAL
        output_file = std::make_shared<WriterRootTreeOPAL>(ai.inputs[1]);
        (std::dynamic_pointer_cast<WriterRootTreeOPAL>(output_file))->init_branches();
        if (options.find("Run") != options.end()) (std::dynamic_pointer_cast<WriterRootTreeOPAL>(output_file))->set_run_number(std::atoi(options.at("Run").c_str()));
        break;
#else
        printf("Output format %s  is not supported\n",ai.output_format_arg);
        exit(2);
        break;
#endif
    case hpezeus:
#ifdef HEPMCCONVERT_EXTENSION_HEPEVTZEUS
        output_file = std::make_shared<WriterHEPEVTZEUS>(ai.inputs[1]);
        break;
#else
        printf("Output format %s  is not supported\n",ai.output_format_arg);
        exit(2);
#endif
    case dot:
#ifdef HEPMCCONVERT_EXTENSION_DOT
        output_file = std::make_shared<WriterDOT>(ai.inputs[1]);
        if (options.find("Style") != options.end()) (std::dynamic_pointer_cast<WriterDOT>(output_file))->set_style(std::atoi(options.at("Style").c_str()));
        break;
#else
        printf("Output format %s  is not supported\n",ai.output_format_arg);
        exit(2);
        break;
#endif
    case plugin:
        if (options.find("OutputPluginLibrary") == options.end())         {
            printf("OutputPluginLibrary option required, e.g. OutputPluginLibrary=libAnalysis.so\n");
            exit(2);
        }
        else OutputPluginLibrary = options.at("OutputPluginLibrary");
        if (options.find("OutputPluginName") == options.end())            {
            printf("OutputPluginName option required, e.g. OutputPluginName=newAnalysisExamplefile\n");
            exit(2);
        }
        else OutputPluginName = options.at("OutputPluginName");
        output_file = std::make_shared<WriterPlugin>(std::string(ai.inputs[1]), OutputPluginLibrary, OutputPluginName);
        if (output_file->failed()) {
            printf("Plugin initialization failed\n");
            exit(2);
        }
        break;
    case dump:
        output_file = NULL;
        break;
    case none:
        output_file = NULL;
        ignore_writer = true;
        break;
    default:
        printf("Output format %s  is not known\n", ai.output_format_arg);
        exit(2);
        break;
    }
    while( !input_file->failed() )
    {
        GenEvent evt(Units::GEV, Units::MM);
        input_file->read_event(evt);
        if( input_file->failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        if (evt.event_number() < first_event_number) continue;
        if (evt.event_number() > last_event_number) continue;
        evt.set_run_info(input_file->run_info());
        //Note the difference between ROOT and Ascii readers. The former read GenRunInfo before first event and the later at the same time as first event.
        if (!ignore_writer)
        {
            if (output_file)
            {
                output_file->write_event(evt);
            }
            else
            {
                Print::content(evt);
            }
        }
        evt.clear();
        ++events_parsed;
        if( events_parsed%print_each_events_parsed == 0 ) printf("Events parsed: %li\n", events_parsed);
        if( events_parsed >= events_limit ) {
            printf("Event limit reached:->events_parsed(%li) >= events_limit(%li)<-. Exit.\n", events_parsed, events_limit);
            break;
        }
    }

    if (input_file)   input_file->close();
    if (output_file)  output_file->close();
    cmdline_parser_free(&ai);
    return EXIT_SUCCESS;
}
