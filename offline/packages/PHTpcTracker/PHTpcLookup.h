/*!
 *  \file       PHTpcLookup.h
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

#ifndef PHTPCLOOKUP_H_
#define PHTPCLOOKUP_H_

#include <trackbase/TrkrClusterContainer.h>
#include "externals/kdfinder.hpp"
#include <vector>

/// \class PHTpcLookup
///
/// \brief 
///
class PHTpcLookup
{
	public:
		PHTpcLookup();
		~PHTpcLookup();

		void init( TrkrClusterContainer* cluster_map );
		void clear();

		std::vector<std::vector<double>*> find( double x, double y, double z, double radius, size_t& nMatches );

	protected:
		TrkrClusterContainer* mClusterMap;
		std::vector<std::vector<double> > mKDhits;
		kdfinder::KDPointCloud<double> mCloud;
		nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, kdfinder::KDPointCloud<double> >,
			kdfinder::KDPointCloud<double>,3>* mKDindex;
	
	private:

};

#endif /* PHTPCLOOKUP_H_ */