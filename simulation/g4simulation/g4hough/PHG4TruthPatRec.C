/*!
 *  \file		PHG4TruthPatRec.C
 *  \brief		Truth Pattern Recognition
 *  \details	Using generate, GEANT level information to mimic ideal pattern recognition
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderCell_MAPS.h>
#include <g4detectors/PHG4CylinderGeom_MAPS.h>
#include <g4detectors/PHG4CylinderGeom_Siladders.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4VtxPointv1.h>
#include <GenFit/FieldManager.h>
#include <GenFit/GFRaveVertex.h>
#include <GenFit/GFRaveVertexFactory.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/StateOnPlane.h>
#include <GenFit/Track.h>
#include <phgenfit/Fitter.h>
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/SpacepointMeasurement.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phgeom/PHGeomUtility.h>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "TClonesArray.h"
#include "TMatrixDSym.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TRotation.h"

#include "SvtxCluster.h"
#include "SvtxClusterMap.h"
#include "SvtxTrackState_v1.h"
#include "SvtxHit_v1.h"
#include "SvtxHitMap.h"
#include "SvtxTrack.h"
#include "SvtxTrack_v1.h"
#include "SvtxVertex_v1.h"
#include "SvtxTrackMap.h"
#include "SvtxTrackMap_v1.h"
#include "SvtxVertexMap_v1.h"
#include "PHG4TruthPatRec.h"

#define LogDebug(exp)		std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"
#define LogError(exp)		std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"
#define LogWarning(exp)	std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"

#define WILD_FLOAT -9999.

using namespace std;

