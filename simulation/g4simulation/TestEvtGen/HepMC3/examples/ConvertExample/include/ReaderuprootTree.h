#ifndef READERUPROOTTREE_H
#define READERUPROOTTREE_H
#include "HepMC3/GenEvent.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/Print.h"
#include "HepMC3/Reader.h"
#include "HepMC3/Data/GenEventData.h"
#include "HepMC3/Data/GenRunInfoData.h"
#include <iostream>
#include "HepMC3/Units.h"
#include "HepMC3/Version.h"
#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

namespace HepMC3
{
/** @brief  ReaderuprootTree */
class ReaderuprootTree : public Reader
{
public:
    /** @brief Constructor with tree and branch names*/
    ReaderuprootTree(const std::string &filename,const std::string &treename="hepmc3_tree",const std::string &branchname="hepmc3_event");

    /// @brief skip events
    bool skip(const int)  override;

    /** @brief Read event from file
     *
     *  @param[out] evt Contains parsed event
     */
    bool read_event(GenEvent &evt)   override;

    /** @brief Close file */
    void close()  override;

    /** @brief Get file  error state */
    bool failed()  override;

    ~ReaderuprootTree();
private:
    /** @brief init routine */
    bool init(const std::string &filename);

    int   m_events_count; //!< Events count. Needed to read the tree
    GenEventData* m_event_data; //!< Pointer to structure that holds event data
    GenRunInfoData* m_run_info_data; //!< Pointer to structure that holds run info data
    std::string m_tree_name; //!< Name of TTree
    std::string m_branch_name; //!< Name of TBranch in TTree

    //PyThreadState* m_thread_state;
    PyObject* m_file;                       //!< Python file handler

    PyObject* m_tree;                       //!< Python tree handler.

    PyObject* m_genruninfo;                 //!< Python runInfo handler.

    PyObject* m_access_function;            //!< Python access function for arrays

    PyObject* m_python_module;              //!< Python module

    long int m_tree_getEntries;             //!< number of processed events

    PyObject* get_function(PyObject*, const std::string& );   //!< Get python functions

    PyObject* init_python_module(const std::string&);         //!< Init python module

    template <class T> std::vector<T> get_vector(PyObject * file_name,const std::string& array_name,std::string desired_type=""); //!< Get arrays
};

}
#endif
