// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_HEPEVT_WRAPPER_RUNTIME_STATIC_H
#define HEPMC3_HEPEVT_WRAPPER_RUNTIME_STATIC_H
#include <iostream>
#include <cstdio>
#include <set>
#include <map>
#include <cstring> // memset
#include <algorithm> //min max for VS2017
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/HEPEVT_Helpers.h"

/**
 *  @file HEPEVT_Wrapper_Runtime_Static.h
 *  @brief Definition of \b class HEPEVT_Wrapper_Runtime_Static
 *
 *  @class HepMC3::HEPEVT_Wrapper_Runtime_Static
 *  @brief A static interface to HEPEVT common block implemented to deal with varying block size in runtime
 */
namespace HepMC3
{


class HEPEVT_Wrapper_Runtime_Static
{

//
// Functions
//
public:
    /** @brief Print information from HEPEVT common block */
    static void print_hepevt( std::ostream& ostr = std::cout );
    /** @brief Print particle information */
    static void print_hepevt_particle( int index, std::ostream& ostr = std::cout );
    /** @brief Set all entries in HEPEVT to zero */
    static void zero_everything();
    /** @brief Convert GenEvent to HEPEVT*/
    static bool GenEvent_to_HEPEVT( const GenEvent* evt ) { return GenEvent_to_HEPEVT_static<HEPEVT_Wrapper_Runtime_Static>(evt);};
    /** @brief Convert HEPEVT to GenEvent*/
    static bool HEPEVT_to_GenEvent( GenEvent* evt ) { return HEPEVT_to_GenEvent_static<HEPEVT_Wrapper_Runtime_Static>(evt);};
    /** @brief Tries to fix list of daughters */
    static bool fix_daughters();
private:
    /** @brief  Fortran common block HEPEVT */
    HEPMC3_EXPORT_API static std::shared_ptr<struct HEPEVT_Pointers<double> >  m_hepevtptr;
    /** @brief Block size */
    HEPMC3_EXPORT_API static int m_max_particles;
//
// Accessors
//
public:
    static void   set_max_number_entries( unsigned int size ) { m_max_particles = size; }//!< Set block size
    static void   set_hepevt_address(char *c); //!< Set Fortran block address
    static int    max_number_entries()      { return m_max_particles;                              } //!< Block size
    static int    event_number()            { return *(m_hepevtptr->nevhep);             } //!< Get event number
    static int    number_entries()          { return *(m_hepevtptr->nhep);               } //!< Get number of entries
    static int    status(const int index )       { return m_hepevtptr->isthep[index-1];    } //!< Get status code
    static int    id(const int index )           { return m_hepevtptr->idhep[index-1];     } //!< Get PDG particle id
    static int    first_parent(const int index ) { return m_hepevtptr->jmohep[2*(index-1)]; } //!< Get index of 1st mother
    static int    last_parent(const int index )  { return m_hepevtptr->jmohep[2*(index-1)+1]; } //!< Get index of last mother
    static int    first_child(const int index )  { return m_hepevtptr->jdahep[2*(index-1)]; } //!< Get index of 1st daughter
    static int    last_child(const int index )   { return m_hepevtptr->jdahep[2*(index-1)+1]; } //!< Get index of last daughter
    static double px(const int index )           { return m_hepevtptr->phep[5*(index-1)];   } //!< Get X momentum
    static double py(const int index )           { return m_hepevtptr->phep[5*(index-1)+1];   } //!< Get Y momentum
    static double pz(const int index )           { return m_hepevtptr->phep[5*(index-1)+2];   } //!< Get Z momentum
    static double e(const int index )            { return m_hepevtptr->phep[5*(index-1)+3];   } //!< Get Energy
    static double m(const int index )            { return m_hepevtptr->phep[5*(index-1)+4];   } //!< Get generated mass
    static double x(const int index )            { return m_hepevtptr->vhep[4*(index-1)];   } //!< Get X Production vertex
    static double y(const int index )            { return m_hepevtptr->vhep[4*(index-1)+1];   } //!< Get Y Production vertex
    static double z(const int index )            { return m_hepevtptr->vhep[4*(index-1)+2];   } //!< Get Z Production vertex
    static double t(const int index )            { return m_hepevtptr->vhep[4*(index-1)+3];   } //!< Get production time
    static int    number_parents(const int index );                                   //!< Get number of parents
    static int    number_children(const int index );                                  //!< Get number of children from the range of daughters
    static int    number_children_exact(const int index );                            //!< Get number of children by counting
    static void set_event_number( const int evtno )       { *(m_hepevtptr->nevhep) = evtno;         } //!< Set event number
    static void set_number_entries( const int noentries ) { *(m_hepevtptr->nhep) = noentries;       } //!< Set number of entries
    static void set_status( const int index, const int status ) { m_hepevtptr->isthep[index-1] = status; } //!< Set status code
    static void set_id(const int index, const int id )         { m_hepevtptr->idhep[index-1] = id;      } //!< Set PDG particle id
    static void set_parents( const int index, const int firstparent, const int lastparent );              //!< Set parents
    static void set_children( const int index, const int firstchild, const int lastchild );               //!< Set children
    static void set_momentum( const int index, const double px, const double py, const double pz, const double e );   //!< Set 4-momentum
    static void set_mass( const int index, double mass );                                     //!< Set mass
    static void set_position( const int index, const double x, const double y, const double z, const double t );      //!< Set position in time-space
};

//
// inline definitions
//

inline void HEPEVT_Wrapper_Runtime_Static::set_hepevt_address(char *c) {
    m_hepevtptr = std::make_shared<struct HEPEVT_Pointers<double> >();
    char* x = c;
    m_hepevtptr->nevhep = (int*)x;
    x += sizeof(int);
    m_hepevtptr->nhep = (int*)(x);
    x += sizeof(int);
    m_hepevtptr->isthep = (int*)(x);
    x += sizeof(int)*m_max_particles;
    m_hepevtptr->idhep = (int*)(x);
    x += sizeof(int)*m_max_particles;
    m_hepevtptr->jmohep = (int*)(x);
    x += sizeof(int)*m_max_particles*2;
    m_hepevtptr->jdahep = (int*)(x);
    x += sizeof(int)*m_max_particles*2;
    m_hepevtptr->phep = (double*)(x);
    x += sizeof(double)*m_max_particles*5;
    m_hepevtptr->vhep = (double*)(x);
}

inline void HEPEVT_Wrapper_Runtime_Static::print_hepevt( std::ostream& ostr )
{
    ostr << " Event No.: " << *(m_hepevtptr->nevhep) << std::endl;
    ostr << "  Nr   Type   Parent(s)  Daughter(s)      Px       Py       Pz       E    Inv. M." << std::endl;
    for ( int i = 1; i <= *(m_hepevtptr->nhep); ++i )
    {
        print_hepevt_particle( i, ostr );
    }
}

inline void HEPEVT_Wrapper_Runtime_Static::print_hepevt_particle( int index, std::ostream& ostr )
{
    char buf[255];//Note: the format is fixed, so no reason for complicated treatment

    sprintf(buf, "%5i %6i", index, m_hepevtptr->idhep[index-1]);
    ostr << buf;
    sprintf(buf, "%4i - %4i  ", m_hepevtptr->jmohep[2*(index-1)], m_hepevtptr->jmohep[2*(index-1)+1]);
    ostr << buf;
    sprintf(buf, "%4i - %4i ", m_hepevtptr->jdahep[2*(index-1)], m_hepevtptr->jdahep[2*(index-1)+1]);
    ostr << buf;
    sprintf(buf, "%8.2f %8.2f %8.2f %8.2f %8.2f", m_hepevtptr->phep[5*(index-1)], m_hepevtptr->phep[5*(index-1)+1], m_hepevtptr->phep[5*(index-1)+2],
            m_hepevtptr->phep[5*(index-1)+3], m_hepevtptr->phep[5*(index-1)+4]);
    ostr << buf << std::endl;
}

inline void HEPEVT_Wrapper_Runtime_Static::zero_everything()
{
    *(m_hepevtptr->nevhep) = 0;
    *(m_hepevtptr->nhep) = 0;
    memset(m_hepevtptr->isthep, 0, sizeof(int)*m_max_particles);
    memset(m_hepevtptr->idhep, 0, sizeof(int)*m_max_particles);
    memset(m_hepevtptr->jmohep, 0, sizeof(int)*m_max_particles*2);
    memset(m_hepevtptr->jdahep, 0, sizeof(int)*m_max_particles*2);
    memset(m_hepevtptr->phep, 0, sizeof(double)*m_max_particles*5);
    memset(m_hepevtptr->vhep, 0, sizeof(double)*m_max_particles*4);
}

inline int HEPEVT_Wrapper_Runtime_Static::number_parents( const int index )
{
    return (m_hepevtptr->jmohep[2*(index-1)]) ? (m_hepevtptr->jmohep[2*(index-1)+1]) ? m_hepevtptr->jmohep[2*(index-1)+1]
           -m_hepevtptr->jmohep[2*(index-1)] : 1 : 0;
}

inline int HEPEVT_Wrapper_Runtime_Static::number_children( const int index )
{
    return (m_hepevtptr->jdahep[2*(index-1)]) ? (m_hepevtptr->jdahep[2*(index-1)+1]) ? m_hepevtptr->jdahep[2*(index-1)+1]-m_hepevtptr->jdahep[2*(index-1)] : 1 : 0;
}

inline int HEPEVT_Wrapper_Runtime_Static::number_children_exact( const int index )
{
    int nc = 0;
    for ( int i = 1; i <= *(m_hepevtptr->nhep); ++i )
        if (((m_hepevtptr->jmohep[2*(i-1)] <= index && m_hepevtptr->jmohep[2*(i-1)+1] >= index)) || (m_hepevtptr->jmohep[2*(i-1)] == index) ||
                (m_hepevtptr->jmohep[2*(index-1)+1]==index)) nc++;
    return nc;
}


inline void HEPEVT_Wrapper_Runtime_Static::set_parents( const int index,  const int firstparent, const int lastparent )
{
    m_hepevtptr->jmohep[2*(index-1)] = firstparent;
    m_hepevtptr->jmohep[2*(index-1)+1] = lastparent;
}

inline void HEPEVT_Wrapper_Runtime_Static::set_children(  const int index,  const int  firstchild,  const int lastchild )
{
    m_hepevtptr->jdahep[2*(index-1)] = firstchild;
    m_hepevtptr->jdahep[2*(index-1)+1] = lastchild;
}

inline void HEPEVT_Wrapper_Runtime_Static::set_momentum( const int index, const double px, const double py, const double pz, const double e )
{
    m_hepevtptr->phep[5*(index-1)] = px;
    m_hepevtptr->phep[5*(index-1)+1] = py;
    m_hepevtptr->phep[5*(index-1)+2] = pz;
    m_hepevtptr->phep[5*(index-1)+3] = e;
}

inline void HEPEVT_Wrapper_Runtime_Static::set_mass( const int index, double mass )
{
    m_hepevtptr->phep[5*(index-1)+4] = mass;
}

inline void HEPEVT_Wrapper_Runtime_Static::set_position( const int index, const double x, const double y, const double z, const double t )
{
    m_hepevtptr->vhep[4*(index-1)] = x;
    m_hepevtptr->vhep[4*(index-1)+1] = y;
    m_hepevtptr->vhep[4*(index-1)+2] = z;
    m_hepevtptr->vhep[4*(index-1)+3] = t;
}

inline bool HEPEVT_Wrapper_Runtime_Static::fix_daughters()
{
    /*AV The function should be called  for a record that has correct particle ordering and mother ids.
    As a result it produces a record with ranges where the daughters can be found.
    Not every particle in the range will be a daughter. It is true only for proper events.
    The return tells if the record was fixed succesfully.
    */
    for ( int i = 1; i <= number_entries(); i++ )
        for ( int k=1; k <= number_entries(); k++ ) if (i != k)
                if ((first_parent(k) <= i) && (i <= last_parent(k)))
                    set_children(i, (first_child(i) == 0 ? k : std::min(first_child(i), k)), (last_child(i) == 0 ? k : std::max(last_child(i), k)));
    bool is_fixed = true;
    for ( int i = 1; i <= number_entries(); i++ )
        is_fixed = (is_fixed && (number_children_exact(i) == number_children(i)));
    return is_fixed;
}

} // namespace HepMC3
#endif
