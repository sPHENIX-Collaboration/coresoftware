#ifndef INTT_LOADABLE_H
#define INTT_LOADABLE_H

/// author: Joseph Bertaux
/// email:  jbertau@purdue.edu

#include <string>

class CDBTTree;

/// Base class for INTT calibrations loaded from CDBTTrees
class InttLoadable
{
public:
	InttLoadable() = default;
	virtual ~InttLoadable() = default;

	/// Returns true if the internal calibration is stored
	bool IsLoaded() {return m_is_loaded;}

	/// Loads a calibration from either the CDB or from a file
	/// If the argument contains the extension ".root", it treats it as the full path to a ROOT file containing a CDBTTree
	/// If the argument does not contain the extension, it treats it as tag and loads the CDBTTree from the CDB
	int Load(std::string const&);

protected:
	/// Member function which should be overwritten by derived classes assuming specific structure of the CDBTTree
	/// This function should return 0 on success and 1 on failure
	virtual int LoadFromCdbTTree(CDBTTree&);

private:
	bool m_is_loaded{false};
};

#endif//INTT_LOADABLE_H
