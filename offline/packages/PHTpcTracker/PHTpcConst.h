/*!
 *  \file       PHTpcConst.h
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

#ifndef PHTPCCONST_H_
#define PHTPCCONST_H_

namespace PHTpcConst {
    // ----- TPC geometry constants ----- 
    const int TPC_LAYERS_MAX = 48;
    const double TPC_LAYERS_RADIUS[ TPC_LAYERS_MAX ] = {
		// 16 layers 0.625 apart
		30.3125, 30.9375, 31.5625, 32.1875,
		32.8125, 33.4375, 34.0625, 34.6875,
		35.3125, 35.9375, 36.5625, 37.1875,
		37.8125, 38.4375, 39.0625, 39.6875,
		// 32 layers 1.25 apart
		40.625, 41.875, 43.125, 44.375,
		45.625, 46.875, 48.125, 49.375,
		50.625, 51.875, 53.125, 54.375,
		55.625, 56.875, 58.125, 59.375,
		60.5313, 61.5938, 62.6562, 63.7187,
		64.7812, 65.8438, 66.9063, 67.9688,
		69.0312, 70.0937, 71.1563, 72.2187,
		73.2812, 74.3437, 75.4063, 76.4687
	};

    const double TPC_LAYERS_DELTAR[ TPC_LAYERS_MAX ] = {
		0.625, 0.625, 0.625, 0.625,
		0.625, 0.625, 0.625, 0.625,
		0.625, 0.625, 0.625, 0.625,
		0.625, 0.625, 0.625, 0.625,
		1.25, 1.25, 1.25, 1.25,
		1.25, 1.25, 1.25, 1.25,
		1.25, 1.25, 1.25, 1.25,
		1.25, 1.25, 1.25, 1.25,
		1.25, 1.25, 1.25, 1.25,
		1.25, 1.25, 1.25, 1.25,
		1.25, 1.25, 1.25, 1.25,
		1.25, 1.25, 1.25, 1.25
	};
    const double TPC_RADIUS_MIN = TPC_LAYERS_RADIUS[ 0 ] - TPC_LAYERS_DELTAR[ 0 ];
    const double TPC_RADIUS_MAX = TPC_LAYERS_RADIUS[ TPC_LAYERS_MAX - 1 ] + TPC_LAYERS_DELTAR[ TPC_LAYERS_MAX - 1 ];
}

#endif // PHTpcConst