#include "ReaderuprootTree.h"
namespace HepMC3
{

HEPMC3_DECLARE_READER_FILE(ReaderuprootTree)

/// @brief obtain vector of objects using name and type
template <class T> std::vector<T>
ReaderuprootTree::get_vector(PyObject * file_name, const std::string& array_name, std::string desired_type) {
    int i = m_events_count;
    PyObject *pFunc = m_access_function;
    PyObject * pArgs = PyTuple_New(4);
    PyTuple_SetItem(pArgs, 0, file_name);
    PyTuple_SetItem(pArgs, 1, Py_BuildValue("s#", array_name.c_str(), array_name.length()));
    PyTuple_SetItem(pArgs, 2, Py_BuildValue("i", i));
    PyTuple_SetItem(pArgs, 3, Py_BuildValue("s#", desired_type.c_str(), desired_type.length()));
    PyObject *pReturn = PyObject_CallObject(pFunc, pArgs);
    PyArrayObject *np_ret = reinterpret_cast<PyArrayObject*>(pReturn);
    std::vector<T> out;
    int len0 = 0;
    if (np_ret) len0 = PyArray_SHAPE(np_ret)[0];
    if (len0 > 0) {
        int len = PyArray_SHAPE(np_ret)[1];
        T*  c_out = reinterpret_cast<T*>(PyArray_DATA(np_ret));
        for (int i = 0; i < len; i++) out.push_back(c_out[i]);
    }
    Py_DECREF(pArgs);
    if (np_ret) Py_DECREF(np_ret);
    return out;
}

/// @brief obtain vector of objects using name and type, specified for std::string
template <>
std::vector<std::string>
ReaderuprootTree::get_vector<std::string>(PyObject * file_name, const std::string& array_name, std::string desired_type) {
    if (desired_type.length() == 0) desired_type = "U500";
    int i = m_events_count;
    PyObject *pFunc = m_access_function;
    PyObject * pArgs = PyTuple_New(4);
    PyTuple_SetItem(pArgs, 0, file_name);
    PyTuple_SetItem(pArgs, 1, Py_BuildValue("s#", array_name.c_str(), array_name.length()));
    PyTuple_SetItem(pArgs, 2, Py_BuildValue("i", i));
    PyTuple_SetItem(pArgs, 3, Py_BuildValue("s#", desired_type.c_str(), desired_type.length()));

    PyObject *pReturn = PyObject_CallObject(pFunc, pArgs);
    PyArrayObject *np_ret = reinterpret_cast<PyArrayObject*>(pReturn);
    std::vector<std::string> out;
    int len0 = 0;
    if (np_ret) len0 = PyArray_SHAPE(np_ret)[0];
    if (len0>0) {
        int len = PyArray_SHAPE(np_ret)[1];
        typedef wchar_t wc500[500];
        wc500* c_out = reinterpret_cast<wc500*>(PyArray_DATA(np_ret));

        for (int i = 0; i < len; i++) {
            std::wstring wa((c_out[i]));
            std::string ret(wa.begin(), wa.end() );
            out.push_back(ret);
        }
    }
    Py_DECREF(pArgs);
    if (np_ret) Py_DECREF(np_ret);
    return out;
}

PyObject* ReaderuprootTree::get_function(PyObject* m_python_module, const std::string& name)
{
    if (!m_python_module) return nullptr;
    PyObject* pFuncInitFile = PyObject_GetAttrString(m_python_module, name.c_str());
    if (!pFuncInitFile || !PyCallable_Check(pFuncInitFile)) {
        Py_XDECREF(pFuncInitFile);
        std::cout << name<<"is null or not callable" << std::endl;
        return nullptr;
    }
    return pFuncInitFile;
}


PyObject* ReaderuprootTree::init_python_module(const std::string& code)
{
    const char *SomeModuleName = "uproot4forhepmc3";
    const char *SomeModuleCode = code.c_str();
    PyObject *m_python_module = PyModule_New(SomeModuleName);
    PyModule_AddStringConstant(m_python_module, "__file__", "");
    PyObject *localDict = PyModule_GetDict(m_python_module);
    PyObject *builtins = PyEval_GetBuiltins();
    PyDict_SetItemString(localDict, "__builtins__", builtins);

    PyObject *pyValue = PyRun_String(SomeModuleCode, Py_file_input, localDict, localDict);
    if (pyValue == nullptr) {
        return nullptr;
    }
    else
    {
        Py_DECREF(pyValue);
    }
    return m_python_module;
}

ReaderuprootTree::ReaderuprootTree(const std::string &filename,const std::string &treename,const std::string &branchname):
    m_events_count(0),m_tree_name(treename.c_str()), m_branch_name(branchname.c_str()),m_tree(nullptr)
{
    if (!init(filename)) return;
}

bool ReaderuprootTree::init(const std::string &filename)
{

    m_event_data = new GenEventData();

    m_run_info_data = new GenRunInfoData();

    set_run_info(std::make_shared<GenRunInfo>());

    Py_Initialize();
    import_array()

    m_python_module = init_python_module(

                          R"EOT(
import uproot
import numpy
def init_file(filename):
    rootfile=uproot.open(str(filename))
#    print(rootfile.keys())
    return rootfile

def close_file(filename):
    return filename.close()

def init_tree(rootfile,treename,branchname):
    tree=rootfile[str(treename)+"/"+str(branchname)]
#    print(tree.keys())
    return tree

def init_genruninfo(rootfile,treename,branchname):
    tree=rootfile[str(treename)+"/"+str(branchname)]
#    print(tree.keys())
    return tree

def get_number_of_entries_in_tree(tree):
    return tree.num_entries

def get_array_from_tree(tree,branch,i,destype):
    result=tree[str(branch)].array(library="np")[i]
    if len(destype.strip()) == 0:
     output=numpy.array([result])
    else:
     output=numpy.array([result], dtype=destype)
#    print("a.shape={}, a.dtype={}".format(output.shape, output.dtype))
#    print(branch,output)
    return output
)EOT"
);
     bool result =false;
     if (!m_python_module) {
     HEPMC3_ERROR( "ReaderuprootTree: cannot initialize python modulr. Please check your uproot and/or numpy instalation.")
     return result;
     }
     PyObject *pFuncInitFile = nullptr;
     PyObject *pFuncInitTree = nullptr;
     PyObject* pFuncEntries = nullptr;


     PyObject * pArgsFile = nullptr;
     PyObject * pArgsTree = nullptr;
     PyObject * pArgsEntries = nullptr;
     PyObject * pArgsGenRunInfo = nullptr;

     m_access_function = get_function(m_python_module, "get_array_from_tree");
     pFuncInitFile = get_function(m_python_module, "init_file");
     pFuncInitTree = get_function(m_python_module, "init_tree");
     pFuncEntries = get_function(m_python_module, "get_number_of_entries_in_tree");



if (m_access_function && pFuncInitFile && pFuncInitTree && pFuncEntries) {
    pArgsFile = PyTuple_New(1);
    PyTuple_SetItem(pArgsFile, 0, Py_BuildValue("s#", filename.c_str(), filename.length()));
    m_file=PyObject_CallObject(pFuncInitFile,pArgsFile);


    if (m_file){
    pArgsTree = PyTuple_New(3);
    PyTuple_SetItem(pArgsTree, 0, m_file);
    PyTuple_SetItem(pArgsTree, 1, Py_BuildValue("s#", m_tree_name.c_str(), m_tree_name.length()));
    PyTuple_SetItem(pArgsTree, 2, Py_BuildValue("s#", m_branch_name.c_str(), m_branch_name.length()));
    m_tree = PyObject_CallObject(pFuncInitTree, pArgsTree);

    pArgsGenRunInfo = PyTuple_New(3);
    PyTuple_SetItem(pArgsGenRunInfo, 0, m_file);
    PyTuple_SetItem(pArgsGenRunInfo, 1, Py_BuildValue("s#", m_tree_name.c_str(), m_tree_name.length()));
    PyTuple_SetItem(pArgsGenRunInfo, 2, Py_BuildValue("s#", "GenRunInfo", 10));
    m_genruninfo = PyObject_CallObject(pFuncInitTree, pArgsGenRunInfo);
    }

    if (m_tree){
    m_tree_getEntries = 0;
    pArgsEntries = PyTuple_New(1);
    PyTuple_SetItem(pArgsEntries, 0, m_tree);
    PyObject * ret = PyObject_CallObject(pFuncEntries,pArgsEntries);
    m_tree_getEntries = PyLong_AsLong(ret);
    Py_DECREF(ret);
    }
   if (m_tree && m_file && m_genruninfo && m_tree_getEntries) result=true;
}

    Py_DECREF(pFuncInitFile);
    Py_DECREF(pFuncInitTree);
    Py_DECREF(pArgsEntries);

    Py_DECREF(pFuncEntries);
    Py_DECREF(pArgsFile);
    return result;
}

bool ReaderuprootTree::skip(const int n)
{
    m_events_count+=n;
    if (m_events_count>m_tree_getEntries) return false;
    return true;
}



bool ReaderuprootTree::read_event(GenEvent& evt)
{
    if (!m_python_module) return false;
    if (m_events_count >= m_tree_getEntries) { m_events_count++; return false;}
    m_event_data->particles.clear();
    m_event_data->vertices.clear();
    m_event_data->links1.clear();
    m_event_data->links2.clear();
    m_event_data->attribute_id.clear();
    m_event_data->attribute_name.clear();
    m_event_data->attribute_string.clear();

    auto event_number_v  = get_vector<int>(m_tree, "event_number");
    if (event_number_v.size() == 0) { m_events_count++; return false;}
    auto weights = get_vector<double>(m_tree, "weights");
    auto event_pos_1_v = get_vector<double>(m_tree, "event_pos/event_pos.m_v1");
    if (event_pos_1_v.size() == 0) { m_events_count++; return false;}
    auto event_pos_2_v = get_vector<double>(m_tree,"event_pos/event_pos.m_v2");
    if (event_pos_2_v.size() == 0) { m_events_count++; return false;}
    auto event_pos_3_v = get_vector<double>(m_tree, "event_pos/event_pos.m_v3");
    if (event_pos_3_v.size() == 0) { m_events_count++; return false;}
    auto event_pos_4_v = get_vector<double>(m_tree, "event_pos/event_pos.m_v4");
    if (event_pos_4_v.size() == 0) { m_events_count++; return false;}
    auto momentum_unit_v = get_vector<int>(m_tree, "momentum_unit");
    if (momentum_unit_v.size() == 0) { m_events_count++; return false;}
    auto length_unit_v = get_vector<int>(m_tree, "length_unit");
    if (length_unit_v.size() == 0) { m_events_count++; return false;}

    auto event_number    = event_number_v.at(0);
    auto event_pos_1     = event_pos_1_v.at(0);
    auto event_pos_2     = event_pos_2_v.at(0);
    auto event_pos_3     = event_pos_3_v.at(0);
    auto event_pos_4     = event_pos_4_v.at(0);
    auto momentum_unit   = momentum_unit_v.at(0);
    auto length_unit     = length_unit_v.at(0);

    auto links1                   = get_vector<int>(m_tree, "links1");
    auto links2                   = get_vector<int>(m_tree, "links2");
    auto attribute_id             = get_vector<int>(m_tree, "attribute_id");
    auto attribute_name           = get_vector<std::string>(m_tree, "attribute_name");
    auto attribute_string         = get_vector<std::string>(m_tree, "attribute_string");
    auto particlesmomentumm_v1    = get_vector<double>(m_tree, "particles/particles.momentum.m_v1");
    auto particlesmomentumm_v2    = get_vector<double>(m_tree, "particles/particles.momentum.m_v2");
    auto particlesmomentumm_v3    = get_vector<double>(m_tree, "particles/particles.momentum.m_v3");
    auto particlesmomentumm_v4    = get_vector<double>(m_tree, "particles/particles.momentum.m_v4");
    auto particlesmass            = get_vector<double>(m_tree, "particles/particles.mass");
    auto particlesis_mass_set     = get_vector<bool>(m_tree, "particles/particles.is_mass_set");
    auto particlesparticlespid    = get_vector<int>(m_tree, "particles/particles.pid");
    auto particlesparticlesstatus = get_vector<int>(m_tree, "particles/particles.status");

    auto verticespositionm_v1    = get_vector<double>(m_tree, "vertices/vertices.position.m_v1");
    auto verticespositionm_v2    = get_vector<double>(m_tree, "vertices/vertices.position.m_v2");
    auto verticespositionm_v3    = get_vector<double>(m_tree, "vertices/vertices.position.m_v3");
    auto verticespositionm_v4    = get_vector<double>(m_tree, "vertices/vertices.position.m_v4");
    auto verticesverticesstatus  = get_vector<int>(m_tree, "vertices/vertices.status");

    m_event_data->event_number = event_number;
    m_event_data->momentum_unit = momentum_unit == 0?HepMC3::Units::MEV:HepMC3::Units::GEV;
    m_event_data->length_unit = length_unit == 0?HepMC3::Units::MM:HepMC3::Units::CM;
    m_event_data->event_pos = HepMC3::FourVector(event_pos_1, event_pos_2, event_pos_3, event_pos_4) ;
    m_event_data->links1 = links1;
    m_event_data->links2 = links2;

    for (size_t k=0; k < particlesparticlespid.size(); k++)
    {
        HepMC3::GenParticleData p = { particlesparticlespid[k], particlesparticlesstatus[k], particlesis_mass_set[k], particlesmass[k],
                                     HepMC3::FourVector(particlesmomentumm_v1[k], particlesmomentumm_v1[k], particlesmomentumm_v1[k], particlesmomentumm_v1[k])
                                    };
        m_event_data->particles.push_back(p);
    }

    for (size_t k=0; k < verticesverticesstatus.size(); k++)
    {
        HepMC3::GenVertexData v = { verticesverticesstatus[k], HepMC3::FourVector(verticespositionm_v1[k], verticespositionm_v2[k], verticespositionm_v3[k], verticespositionm_v4[k])};
        m_event_data->vertices.push_back(v);
    }
    m_event_data->weights=weights;
    m_event_data->attribute_id=attribute_id;
    m_event_data->attribute_name=attribute_name;
    m_event_data->attribute_string=attribute_string;
    evt.read_data(*m_event_data);

    m_run_info_data->weight_names.clear();
    m_run_info_data->tool_name.clear();
    m_run_info_data->tool_version.clear();
    m_run_info_data->tool_description.clear();
    m_run_info_data->attribute_name.clear();
    m_run_info_data->attribute_string.clear();

    m_run_info_data->weight_names       =  get_vector<std::string>(m_genruninfo, "weight_names");
    m_run_info_data->tool_name          =  get_vector<std::string>(m_genruninfo, "tool_name");
    m_run_info_data->tool_version       =  get_vector<std::string>(m_genruninfo, "tool_version");
    m_run_info_data->tool_description   =  get_vector<std::string>(m_genruninfo, "tool_description");
    m_run_info_data->attribute_name     =  get_vector<std::string>(m_genruninfo, "attribute_name");
    m_run_info_data->attribute_string   =  get_vector<std::string>(m_genruninfo, "attribute_string");

    run_info()->read_data(*m_run_info_data);
    evt.set_run_info(run_info());
    m_events_count++;
    return true;
}

void ReaderuprootTree::close()
{
    Py_DECREF(m_genruninfo);
    Py_DECREF(m_tree);
    Py_DECREF(m_file); //This should close the file, right?
    Py_DECREF(m_access_function);
    Py_DECREF(m_python_module);

    m_file = nullptr;
    m_tree = nullptr;
    m_genruninfo = nullptr;
    m_access_function = nullptr;
    m_python_module = nullptr;
    Py_DECREF( PyImport_ImportModule("threading")); //If someone, at some point would document it in CPython...
    Py_Finalize();
}

bool ReaderuprootTree::failed()
{
    if (!m_python_module) return true;
    if (m_events_count >= m_tree_getEntries) return true;
    return false;
}
ReaderuprootTree::~ReaderuprootTree()
{
    if (m_event_data) {delete m_event_data; m_event_data=nullptr;}
    if (m_run_info_data) {delete m_run_info_data; m_run_info_data=nullptr;}
}

} // namespace HepMC3


