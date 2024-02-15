#ifndef INTT_DEAD_MAPv1_H
#define INTT_DEAD_MAPv1_H

#include "InttDeadMap.h"
#include "InttMap.h"

#include <iostream>
#include <set>
#include <string>

#include <phool/PHObject.h>

class InttDeadMapv1 : public InttDeadMap {
public:
	InttDeadMapv1();
	~InttDeadMapv1() override;

	void identify(std::ostream& = std::cout) const override;
	std::size_t size() const override;

	bool IsDeadChannel(InttMap::Offline_s const&) const override;

protected:
	int v_LoadFromCDBTTree(CDBTTree&) override;

private:
	typedef std::set<InttMap::Offline_s, InttMap::OfflineWildcardComparator> Set_t;
	Set_t* m_HotChannelSet = nullptr;

	ClassDefOverride(InttDeadMapv1, 1)
};

#endif//INTT_DEAD_MAPv1_H
