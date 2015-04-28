#include "PHG4Hit.h"

using namespace std;

ClassImp(PHG4Hit)

void
PHG4Hit::Copy(PHG4Hit const &g4hit)
{
  for (int i =0; i<2; i++)
    {
      set_x(i,g4hit.get_x(i));
      set_y(i,g4hit.get_y(i));
      set_z(i,g4hit.get_z(i));
      set_px(i,g4hit.get_px(i));
      set_py(i,g4hit.get_py(i));
      set_pz(i,g4hit.get_pz(i));
      set_t(i,g4hit.get_t(i));
    }
  set_edep(g4hit.get_edep());
  set_eion(g4hit.get_eion());
  set_path_length(g4hit.get_path_length());
  set_light_yield(g4hit.get_light_yield());
  set_layer(g4hit.get_layer());
  set_hit_id(g4hit.get_hit_id());
  set_scint_id(g4hit.get_scint_id());
  set_trkid(g4hit.get_trkid());
  set_strip_z_index(g4hit.get_strip_z_index());
  set_strip_y_index(g4hit.get_strip_y_index());
  set_ladder_z_index(g4hit.get_ladder_z_index());
  set_ladder_phi_index(g4hit.get_ladder_phi_index());
}


void
PHG4Hit::identify(ostream& os) const
{
  cout << "Class " << this->ClassName() << endl;
  cout << "x0: " << get_x(0)
       << ", y0: " << get_y(0)
       << ", z0: " << get_z(0)
       << ", t0: " << get_t(0) << endl;
  cout << "x1: " << get_x(1)
       << ", y1: " << get_y(1)
       << ", z1: " << get_z(1)
       << ", t1: " << get_t(1) << endl;
  cout     << "trackid: " << get_trkid() << ", edep: " << get_edep() << endl;
  cout     << "strip_z_index: " << get_strip_z_index() << ", strip_y_index: " << get_strip_y_index() << endl;
  cout     << "ladder_z_index: " << get_ladder_z_index() << ", ladder_phi_index: " << get_ladder_phi_index() << endl;
  cout << "layer id: " << get_layer() << ", scint_id: " << get_scint_id() << endl;
  return;
}

ostream& operator<<(ostream& stream, const PHG4Hit * hit){
  stream <<endl<< "(x,y,z) = "<< "("<<hit->get_avg_x()<<", "<<hit->get_avg_y()<<", "<<hit->get_avg_z()<<")"<<endl;
	stream << "trackid: " << hit->get_trkid()<<" hitid: "<<hit->get_hit_id()<<" layer: "<<hit->get_layer()<<endl;
  return stream;
}

