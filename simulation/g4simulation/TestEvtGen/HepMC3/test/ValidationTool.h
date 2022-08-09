// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef VALIDATION_TOOL_H
#define VALIDATION_TOOL_H

#ifdef HEPMC2

// Ignore HepMC3 code, use HepMC2 code
#define HEPMC2CODE( x ) x
#define HEPMC3CODE( x )
using namespace HepMC;

#else

// Ignore HepMC2 code, use HepMC3 code
#define HEPMC2CODE( x )
#define HEPMC3CODE( x ) x
using namespace HepMC3;

#endif // ifdef HEPMC2
/// @class ValidationTool
/// @brief  Virtual Interface to validation tools
class ValidationTool {
//
// Constructors
//
public:
    /** Virtual destructor */
    virtual ~ValidationTool() {};

//
// Abstract functions
//
public:
    /** @brief Get information if this tool modifies the event
     *
     *  Tools that do not modify event will be ignored during event printing
     *  and momentum conservation checks
     */
    virtual bool  tool_modifies_event()    = 0;

    /** @brief Get name of the tool */
    virtual const std::string name()       = 0;

    virtual void  initialize()             = 0; //!< Initialize
    virtual int   process(GenEvent &hepmc) = 0; //!< Process event
    virtual void  finalize()               = 0; //!< Finalize

//
// Virtual functions
//
public:
    /** @brief Get long name of the tool */
    virtual const std::string long_name()  { return name(); }

    /** @brief Get timer for this tool (if this tool is being timed)
     *
     *  Note that normally the tool itself should not use the timer it provides
     *  However, if one want to exclude some part of initialization
     *  timer()->start() can be used to restart the timer per each event
     */
    virtual class Timer* timer() { return NULL; }
};

#endif
