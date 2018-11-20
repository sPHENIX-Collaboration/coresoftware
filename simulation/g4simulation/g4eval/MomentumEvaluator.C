#include "MomentumEvaluator.h"


#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrack.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <TNtuple.h>
#include <TFile.h>

#include <vector>
#include <iostream>

using namespace std;

class TrivialTrack
{
	public:
		float px,py,pz;
		float dcax, dcay, dcaz;
		float quality;
		TrivialTrack(float x,float y,float z,float dx,float dy,float dz, float qual=0.) : px(x), py(y), pz(z), dcax(dx), dcay(dy), dcaz(dz), quality(qual) {}
		~TrivialTrack(){}
};


class RecursiveMomentumContainer
{
	protected:
		float px_lo, px_hi;
		float py_lo, py_hi;
		float pz_lo, pz_hi;

		int level;
		int maxlevel;

		unsigned int x_pos,y_pos,z_pos;

		RecursiveMomentumContainer* containers[2][2][2];

	public:
		RecursiveMomentumContainer( float PX_LO, float PX_HI, float PY_LO, float PY_HI, float PZ_LO, float PZ_HI, int MLEV, int LEV=0 ) : px_lo(PX_LO), px_hi(PX_HI), py_lo(PY_LO), py_hi(PY_HI), pz_lo(PZ_LO), pz_hi(PZ_HI), level(LEV), maxlevel(MLEV), x_pos(0), y_pos(0), z_pos(0)
		{
			for(unsigned int i=0;i<2;++i){
				for(unsigned int j=0;j<2;++j){
					for(unsigned int k=0;k<2;++k){containers[i][j][k]=nullptr;}}}
		}
		virtual ~RecursiveMomentumContainer()
		{
			for(unsigned int i=0;i<2;++i){
				for(unsigned int j=0;j<2;++j){
					for(unsigned int k=0;k<2;++k){
						if(containers[i][j][k]!=nullptr)
						{
							delete containers[i][j][k];
						}
					}}}
		}
		
		virtual bool insert( TrivialTrack& track );

		virtual TrivialTrack* begin()
		{
			x_pos = 0;y_pos = 0;z_pos = 0;
			while(true)
			{
				if( containers[x_pos][y_pos][z_pos] == nullptr )
				{
					if(z_pos==0){z_pos=1;continue;}
					else
					{
						if(y_pos==0){z_pos=0;y_pos=1;continue;}
						else
						{
							if(x_pos==0){z_pos=0;y_pos=0;x_pos=1;continue;}
							else{return nullptr;}
						}
					}
				}
				else
				{
					return containers[x_pos][y_pos][z_pos]->begin();
				}
			}
		}

		virtual TrivialTrack* next()
		{
			bool block_changed = false;
			while(true)
			{
				if( containers[x_pos][y_pos][z_pos] == nullptr )
				{
					block_changed = true;
					if(z_pos==0){z_pos=1;continue;}
					else
					{
						if(y_pos==0){z_pos=0;y_pos=1;continue;}
						else
						{
							if(x_pos==0){z_pos=0;y_pos=0;x_pos=1;continue;}
							else{return nullptr;}
						}
					}
				}
				TrivialTrack* val = nullptr;
				if(block_changed == true)
				{
					val = containers[x_pos][y_pos][z_pos]->begin();	
				}
				else
				{
					val = containers[x_pos][y_pos][z_pos]->next();
				}
				
				if(val == nullptr)
				{
					block_changed = true;
					if(z_pos==0){z_pos=1;continue;}
					else
					{
						if(y_pos==0){z_pos=0;y_pos=1;continue;}
						else
						{
							if(x_pos==0){z_pos=0;y_pos=0;x_pos=1;continue;}
							else{return nullptr;}
						}
					}
				}
				else{ return val; }
			}
		}

		virtual void append_list( vector<TrivialTrack*>& track_list, float PX_LO, float PX_HI, float PY_LO, float PY_HI, float PZ_LO, float PZ_HI )
		{
			for(unsigned int i=0;i<2;++i){
				for(unsigned int j=0;j<2;++j){
					for(unsigned int k=0;k<2;++k){
						if(containers[i][j][k]==nullptr){continue;}
						
						if( (containers[i][j][k]->px_hi<PX_LO) || (containers[i][j][k]->px_lo>PX_HI) || (containers[i][j][k]->py_hi<PY_LO) || (containers[i][j][k]->py_lo>PY_HI) || (containers[i][j][k]->pz_hi<PZ_LO) || (containers[i][j][k]->pz_lo>PZ_HI) )
						{
							continue;
						}

						containers[i][j][k]->append_list( track_list, PX_LO, PX_HI, PY_LO, PY_HI, PZ_LO, PZ_HI );
					}}}
			}
};

class RecursiveMomentumContainerEnd : public RecursiveMomentumContainer
{
	public:
		RecursiveMomentumContainerEnd( float PX_LO, float PX_HI, float PY_LO, float PY_HI, float PZ_LO, float PZ_HI, int MLEV, int LEV=0 ) : RecursiveMomentumContainer( PX_LO,PX_HI,PY_LO,PY_HI,PZ_LO,PZ_HI,MLEV,LEV )
		{

		}

		virtual ~RecursiveMomentumContainerEnd()
		{

		}

		virtual bool insert( TrivialTrack& track )
		{
			tracks.push_back(track);
			return true;
		}

		virtual TrivialTrack* begin()
		{
			x_pos = 0;
			return ( &(tracks.at(0)) );
		}

		virtual TrivialTrack* next()
		{
			if( x_pos >= (tracks.size()-1) ){return nullptr;}
			else
			{
				x_pos += 1;
				return (&(tracks[x_pos]));
			}
		}

		virtual void append_list( vector<TrivialTrack*>& track_list, float PX_LO, float PX_HI, float PY_LO, float PY_HI, float PZ_LO, float PZ_HI )
		{
			for(unsigned int i=0;i<tracks.size();++i)
			{
				if( (tracks[i].px<PX_LO) || (tracks[i].px>PX_HI) || (tracks[i].py<PY_LO) || (tracks[i].py>PY_HI) || (tracks[i].pz<PZ_LO) || (tracks[i].pz>PZ_HI) ){continue;}
				track_list.push_back( &(tracks[i]) );
			}
		}

	protected:
		vector<TrivialTrack> tracks;
};


bool RecursiveMomentumContainer::insert( TrivialTrack& track )
{
	if( (track.px < px_lo) || (track.py < py_lo) || (track.pz < pz_lo) || (track.px > px_hi) || (track.py > py_hi) || (track.pz > pz_hi)  )
	{
		return false;
	}

	int x_ind = 0;
	if(track.px > (px_lo + 0.5*(px_hi-px_lo))){x_ind=1;}
	int y_ind = 0;
	if(track.py > (py_lo + 0.5*(py_hi-py_lo))){y_ind=1;}
	int z_ind = 0;
	if(track.pz > (pz_lo + 0.5*(pz_hi-pz_lo))){z_ind=1;}

	if( containers[x_ind][y_ind][z_ind] == nullptr )
	{
		float px_lo_new = px_lo + (float(x_ind))*0.5*(px_hi-px_lo);
		float px_hi_new = px_lo_new + 0.5*(px_hi-px_lo);

		float py_lo_new = py_lo + (float(y_ind))*0.5*(py_hi-py_lo);
		float py_hi_new = py_lo_new + 0.5*(py_hi-py_lo);

		float pz_lo_new = pz_lo + (float(z_ind))*0.5*(pz_hi-pz_lo);
		float pz_hi_new = pz_lo_new + 0.5*(pz_hi-pz_lo);

		if(level < maxlevel)
		{
			containers[x_ind][y_ind][z_ind] = new RecursiveMomentumContainer( px_lo_new,px_hi_new, py_lo_new,py_hi_new, pz_lo_new,pz_hi_new, maxlevel, level+1 );
		}
		else
		{
			containers[x_ind][y_ind][z_ind] = new RecursiveMomentumContainerEnd( px_lo_new,px_hi_new, py_lo_new,py_hi_new, pz_lo_new,pz_hi_new, maxlevel, level+1 );
		}
	}
	return containers[x_ind][y_ind][z_ind]->insert( track );
}


MomentumEvaluator::MomentumEvaluator( std::string fname, float pt_s, float pz_s, unsigned int n_l, unsigned int n_i, unsigned int n_r, float i_z, float o_z ) : ntp_true(nullptr), ntp_reco(nullptr), pt_search_scale(pt_s), pz_search_scale(pz_s), event_counter(0), file_name(fname), n_layers(n_l), n_inner_layers(n_i), n_required_layers(n_r), inner_z_length(i_z), outer_z_length(o_z) {}
MomentumEvaluator::~MomentumEvaluator()
{
	if(ntp_true != nullptr){delete ntp_true;}
	if(ntp_reco != nullptr){delete ntp_reco;}
}

int MomentumEvaluator::Init( PHCompositeNode *topNode )
{
	if(ntp_true != nullptr){delete ntp_true;}
	if(ntp_reco != nullptr){delete ntp_reco;}

	ntp_true = new TNtuple( "ntp_true", "true simulated tracks", "event:px:py:pz:dcax:dcay:dcaz:r_px:r_py:r_pz:r_dcax:r_dcay:r_dcaz:quality" );
	ntp_reco = new TNtuple( "ntp_reco", "reconstructed tracks", "event:px:py:pz:dcax:dcay:dcaz:t_px:t_py:t_pz:t_dcax:t_dcay:t_dcaz:quality" );
	event_counter = 0;

	return Fun4AllReturnCodes::EVENT_OK;
}

int MomentumEvaluator::process_event( PHCompositeNode *topNode )
{
	PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");

	PHG4HitContainer* g4hits = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_SVTX");
	if(g4hits == nullptr){cout<<"can't find PHG4HitContainer"<<endl;exit(1);}
	PHG4HitContainer::ConstRange g4range = g4hits->getHits();

	// set<int> trkids;
	map<int, pair<unsigned int,unsigned int> > trkids;

	for( PHG4HitContainer::ConstIterator iter = g4range.first; iter != g4range.second; ++iter )
	{
		PHG4Hit* hit = iter->second;

		int layer = hit->get_layer();
		float length = outer_z_length;
		if(((unsigned int)layer)<n_inner_layers){length=inner_z_length;}
		if(fabs(hit->get_z(0))>length){continue;}

		int trk_id = hit->get_trkid();
		if(trkids.find(trk_id) == trkids.end())
		{
			trkids[trk_id].first = 0;
			trkids[trk_id].second = 0;
		}
		if( hit->get_layer() < 32 )
		{
		  trkids[trk_id].first = (trkids[trk_id].first | (1<<(hit->get_layer())));
		}
		else
		{
		  trkids[trk_id].second = (trkids[trk_id].second | (1<<(hit->get_layer()-32)));
		}
		
		// cout<<"trk_id = "<<trk_id<<endl;
		// cout<<"layer = "<<hit->get_layer()<<endl;
		// cout<<"nlayer = "<<__builtin_popcount(trkids[trk_id].first)+__builtin_popcount(trkids[trk_id].second)<<endl<<endl;
		// trkids.insert(trk_id);
	}


	SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");

	PHG4VtxPoint *gvertex = truthinfo->GetPrimaryVtx( truthinfo->GetPrimaryVertexIndex() );
	float gvx = gvertex->get_x();
	float gvy = gvertex->get_y();
	float gvz = gvertex->get_z();

	RecursiveMomentumContainer true_sorted( -20., 20., -20., 20., -20., 20., 10 );

	// PHG4TruthInfoContainer::Map primarymap = truthinfo->GetPrimaryMap();
	PHG4TruthInfoContainer::Map primarymap = truthinfo->GetMap();
   for(PHG4TruthInfoContainer::Iterator iter = primarymap.begin();iter != primarymap.end();++iter)
   {
   	PHG4Particle *particle = iter->second;

   	float vx = truthinfo->GetVtx(particle->get_vtx_id())->get_x();
   	float vy = truthinfo->GetVtx(particle->get_vtx_id())->get_y();
   	float vz = truthinfo->GetVtx(particle->get_vtx_id())->get_z();
   	
   	TrivialTrack track( particle->get_px(), particle->get_py(), particle->get_pz(), vx-gvx, vy-gvy, vz-gvz );

   	if( ( (track.px * track.px) + (track.py * track.py) ) < (0.1*0.1) ){continue;}

   	if( trkids.find(particle->get_track_id()) == trkids.end() )
   	{
   		continue;
   	}

   	// cout<<"trk, nhits = "<<particle->get_track_id()<<" "<<__builtin_popcount(trkids[particle->get_track_id()].first)+__builtin_popcount(trkids[particle->get_track_id()].second)<<endl;

   	if( __builtin_popcount(trkids[particle->get_track_id()].first)+__builtin_popcount(trkids[particle->get_track_id()].second) < (int)n_required_layers )
   	{
   		continue;
   	}

   	true_sorted.insert( track );
   }


   RecursiveMomentumContainer reco_sorted( -20., 20., -20., 20., -20., 20., 10 );
   for(SvtxTrackMap::Iter iter = trackmap->begin();iter != trackmap->end();++iter)
   {
   	SvtxTrack* track = iter->second;

   	TrivialTrack ttrack( track->get_px(), track->get_py(), track->get_pz(), track->get_x()-gvx, track->get_y()-gvy, track->get_z()-gvz, track->get_quality() );
   	reco_sorted.insert(ttrack);
   }


   TrivialTrack* t_track = true_sorted.begin();
   vector<TrivialTrack*> pointer_list;
   while(t_track != nullptr)
   {
   	pointer_list.clear();

   	float pt = sqrt((t_track->px * t_track->px) + (t_track->py * t_track->py));
   	float pt_diff = pt*pt_search_scale;
   	float px_lo = t_track->px - pt_diff;
   	float px_hi = t_track->px + pt_diff;
   	float py_lo = t_track->py - pt_diff;
   	float py_hi = t_track->py + pt_diff;
   	float pz_diff = fabs( t_track->pz )*pz_search_scale;
   	float pz_lo = t_track->pz - pz_diff;
   	float pz_hi = t_track->pz + pz_diff;

   	reco_sorted.append_list( pointer_list, px_lo,px_hi, py_lo,py_hi, pz_lo,pz_hi );

   	if(pointer_list.size() > 0)
   	{
   		float mom_true = sqrt(pt*pt + (t_track->pz)*(t_track->pz));
   		float best_ind = 0;
   		float mom_reco = sqrt( (pointer_list[0]->px)*(pointer_list[0]->px) + (pointer_list[0]->py)*(pointer_list[0]->py) + (pointer_list[0]->pz)*(pointer_list[0]->pz) );
   		float best_mom = mom_reco;
   		for(unsigned int i=1;i<pointer_list.size();++i)
   		{
   			mom_reco = sqrt( (pointer_list[i]->px)*(pointer_list[i]->px) + (pointer_list[i]->py)*(pointer_list[i]->py) + (pointer_list[i]->pz)*(pointer_list[i]->pz) );
   			if( fabs( mom_true - mom_reco ) < fabs( mom_true - best_mom )  )
   			{
   				best_mom = mom_reco;
   				best_ind = i;
   			}
   		}
   		
   		float ntp_data[14] = { (float) event_counter, t_track->px, t_track->py, t_track->pz, t_track->dcax, t_track->dcay, t_track->dcaz, pointer_list[best_ind]->px, pointer_list[best_ind]->py, pointer_list[best_ind]->pz, pointer_list[best_ind]->dcax, pointer_list[best_ind]->dcay, pointer_list[best_ind]->dcaz, pointer_list[best_ind]->quality };
   		ntp_true->Fill(ntp_data);
   	}
   	else
   	{
   		float ntp_data[14] = { (float) event_counter, t_track->px, t_track->py, t_track->pz, t_track->dcax, t_track->dcay, t_track->dcaz, -9999.,-9999.,-9999.,-9999.,-9999.,-9999., -9999. };
   		ntp_true->Fill(ntp_data);
   	}

   	t_track = true_sorted.next();
   }

   TrivialTrack* r_track = reco_sorted.begin();
   while(r_track != nullptr)
   {
   	pointer_list.clear();

   	float pt = sqrt((r_track->px * r_track->px) + (r_track->py * r_track->py));
   	float pt_diff = pt*pt_search_scale;
   	float px_lo = r_track->px - pt_diff;
   	float px_hi = r_track->px + pt_diff;
   	float py_lo = r_track->py - pt_diff;
   	float py_hi = r_track->py + pt_diff;
   	float pz_diff = fabs( r_track->pz )*pz_search_scale;
   	float pz_lo = r_track->pz - pz_diff;
   	float pz_hi = r_track->pz + pz_diff;

   	true_sorted.append_list( pointer_list, px_lo,px_hi, py_lo,py_hi, pz_lo,pz_hi );

   	if(pointer_list.size() > 0)
   	{
   		float mom_reco = sqrt(pt*pt + (r_track->pz)*(r_track->pz));
   		float best_ind = 0;
   		float mom_true = sqrt( (pointer_list[0]->px)*(pointer_list[0]->px) + (pointer_list[0]->py)*(pointer_list[0]->py) + (pointer_list[0]->pz)*(pointer_list[0]->pz) );
   		float best_mom = mom_true;
   		for(unsigned int i=1;i<pointer_list.size();++i)
   		{
   			mom_true = sqrt( (pointer_list[i]->px)*(pointer_list[i]->px) + (pointer_list[i]->py)*(pointer_list[i]->py) + (pointer_list[i]->pz)*(pointer_list[i]->pz) );
   			if( fabs( mom_reco - mom_true ) < fabs( mom_reco - best_mom )  )
   			{
   				best_mom = mom_true;
   				best_ind = i;
   			}
   		}
   		
   		float ntp_data[14] = { (float) event_counter, r_track->px, r_track->py, r_track->pz, r_track->dcax, r_track->dcay, r_track->dcaz, pointer_list[best_ind]->px, pointer_list[best_ind]->py, pointer_list[best_ind]->pz, pointer_list[best_ind]->dcax, pointer_list[best_ind]->dcay, pointer_list[best_ind]->dcaz, r_track->quality };
   		ntp_reco->Fill(ntp_data);
   	}
   	else
   	{
   		float ntp_data[14] = { (float) event_counter, r_track->px, r_track->py, r_track->pz, r_track->dcax, r_track->dcay, r_track->dcaz, -9999.,-9999.,-9999.,-9999.,-9999.,-9999., r_track->quality };
   		ntp_reco->Fill(ntp_data);
   	}

   	r_track = reco_sorted.next();
   }


   event_counter += 1;
   return Fun4AllReturnCodes::EVENT_OK;
}

int MomentumEvaluator::End( PHCompositeNode *topNode )
{
	TFile outfile(file_name.c_str(), "recreate");
	outfile.cd();
	ntp_true->Write();
	ntp_reco->Write();
	outfile.Close();

	return Fun4AllReturnCodes::EVENT_OK;
}







