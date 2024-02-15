#ifndef INTT_DEAD_MAP_H
#define INTT_DEAD_MAP_H

#include "InttMap.h"

#include <filesystem>
#include <iostream>
#include <set>
#include <string>

#include <cdbobjects/CDBTTree.h>
#include <g4detectors/PHG4CellDefs.h>
#include <ffamodules/CDBInterface.h>
#include <phool/PHObject.h>

class InttDeadMap : public PHObject {
public:
	InttDeadMap();
	~InttDeadMap() override;

	virtual void identify(std::ostream& = std::cout) const override;
	virtual std::size_t size() const;

	int LoadFromFile(std::string const& = "InttDeadMap.root");
	int LoadFromCDB(std::string const& = "InttDeadMap");

	bool IsDeadChannel(int const&, int const&, int const&, int const&, int const&) const;
	virtual bool IsDeadChannel(InttMap::Offline_s const&) const;

protected:
	virtual int v_LoadFromCDBTTree(CDBTTree&);

private:
	ClassDefOverride(InttDeadMap, 1)
};

#endif//INTT_DEAD_MAP_H
