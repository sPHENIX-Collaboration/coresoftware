// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 * \file PHG4SectorConstructor.h
 * \brief Generalized sector detectors
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.3 $
 * \date $Date: 2014/05/01 19:03:26 $
 */

#ifndef G4DETECTORS_PHG4SECTORCONSTRUCTOR_H
#define G4DETECTORS_PHG4SECTORCONSTRUCTOR_H

#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>  // for G4int

class G4LogicalVolume;
class G4PVPlacement;
class G4VSolid;
class PHG4SectorDisplayAction;
class PHG4Subsystem;

#include <map>
#include <utility>

#include <cassert>
#include <cmath>
#include <string>
#include <vector>

namespace PHG4Sector
{
  //! layer description data
  //! use GEANT units!
  class Layer
  {
   public:
    Layer(  //! name base for this layer
        const std::string &_name,

        //! material name in G4
        const std::string &_material,

        //! depth in G4 units
        double _depth,

        //! percentage filled
        double _percentage_filled,

        //! whether it is active
        bool _active)
      : name(_name)
      , material(_material)
      , depth(_depth)
      , percentage_filled(
            _percentage_filled)
      , active(_active)
    {
    }

   public:
    //! name base for this layer
    std::string name;

    //! material name in G4
    std::string material;

    //! depth in G4 units
    double depth;

    //! percentage filled
    double percentage_filled;

    //! whether it is active
    bool active;
  };

  //! geometry data
  //! use GEANT units! Use
  class Sector_Geometry
  {
   public:
    Sector_Geometry()
    {
      SetDefault();
    }

    void
    SetDefault();

    int get_N_Sector() const
    {
      return N_Sector;
    }

    std::vector<Layer> &
    get_layer_list()
    {
      return layer_list;
    }

    std::string
    get_material() const
    {
      return material;
    }

    double
    get_max_polar_angle() const
    {
      return max_polar_angle;
    }

    double
    get_min_polar_angle() const
    {
      return min_polar_angle;
    }

    double
    get_normal_polar_angle() const
    {
      return normal_polar_angle;
    }

    double
    get_normal_start() const
    {
      return normal_start;
    }

    void
    set_N_Sector(int sector)
    {
      assert(sector >= 1);

      N_Sector = sector;
    }

    void
    set_layer_list(const std::vector<Layer> &layerList)
    {
      layer_list = layerList;
    }

    void
    set_material(const std::string &_material)
    {
      material = _material;
    }

    void
    set_max_polar_angle(double maxPolarAngle)
    {
      assert(maxPolarAngle >= 0);
      assert(maxPolarAngle <= pi);

      max_polar_angle = maxPolarAngle;
    }

    void
    set_min_polar_angle(double minPolarAngle)
    {
      assert(minPolarAngle >= 0);
      assert(minPolarAngle <= pi);

      min_polar_angle = minPolarAngle;
    }

    void
    set_normal_polar_angle(double normalPolarAngle)
    {
      assert(normalPolarAngle >= 0);
      assert(normalPolarAngle <= pi);

      normal_polar_angle = normalPolarAngle;
    }

    void
    set_normal_start(double normalZStart)
    {
      normal_start = normalZStart;
    }

    // derivative constants

    //   ! max radius from IP
    double
    get_max_R() const;

    double
    get_total_thickness() const;

    //!Unit
    static double
    Unit_cm()
    {
      return cm;
    }

    //!Intercept certain z point at certain polar angle
    void
    set_normal_start(const double z_intercept, const double angle_intercept)
    {
      set_normal_start(
          z_intercept / cos(angle_intercept) * cos(normal_polar_angle - angle_intercept));
    }

    //! Pseudorapidity
    static double
    eta_to_polar_angle(const double eta)
    {
      return 2. * atan(exp(-eta));
    }

    // layer descriptions
   public:
    typedef std::vector<Layer> t_layer_list;
    t_layer_list layer_list;

    int GetNumActiveLayers() const
    {
      int n = 0;
      for (t_layer_list::const_iterator it = layer_list.begin();
           it != layer_list.end(); ++it)
        if ((*it).active)
          n++;
      return n;
    }

    void
    AddLayer(                            //
        const std::string &_name,        //! name base for this layer
        const std::string &_material,    //! material name in G4
        double _depth,                   //! depth in G4 units
        bool _active = false,            //! active detector element for sensitive detector?
        double _percentage_filled = 100  //! percentage filled//
    )
    {
      layer_list.push_back(
          Layer(_name, _material, _depth, _percentage_filled, _active));
    }

    //! add Entrace window and drift volume
    //! Ref: P. Abbon et al. The COMPASS experiment at CERN. Nucl. Instrum. Meth., A577:455
    //! 518, 2007. arXiv:hep-ex/0703049, doi:10.1016/j.nima.2007.03.026. 3
    void
    AddLayers_DriftVol_COMPASS(const double drift_vol_thickness = 3 * mm);

    //! add HBD GEM to the layer list,
    //! Ref: W. Anderson et al. Design, Construction, Operation and Performance of a Hadron
    //! Blind Detector for the PHENIX Experiment. Nucl. Instrum. Meth., A646:35 58, 2011.
    //! arXiv:1103.4277, doi:10.1016/j.nima.2011.04.015. 3.5.1
    void
    AddLayers_HBD_GEM(const int n_GEM_layers = 3);

    //! add HBD readout,
    //! Ref: W. Anderson et al. Design, Construction, Operation and Performance of a Hadron
    //! Blind Detector for the PHENIX Experiment. Nucl. Instrum. Meth., A646:35 58, 2011.
    //! arXiv:1103.4277, doi:10.1016/j.nima.2011.04.015. 3.5.1
    void
    AddLayers_HBD_Readout();

    //! Rough AeroGel detector
    //! Ref: T. Iijima et al. A Novel type of proximity focusing RICH counter with multiple
    //! refractive index aerogel radiator. Nucl. Instrum. Meth., A548:383 390, 2005. arXiv:
    //! physics/0504220, doi:10.1016/j.nima.2005.05.030
    void
    AddLayers_AeroGel_ePHENIX(const double radiator_length = 2 * cm,
                              const double expansion_length = 18 * cm, std::string radiator = "Default");

   public:
    typedef enum
    {

      //! cone cut for the polar edge
      kConeEdge = 0,

      //! flat line edge in the azimuthal direction and along the normal_polar_angle
      kFlatEdge = 1

    } e_edge_typ;

    static e_edge_typ
    ConeEdge()
    {
      return kConeEdge;
    }

    static e_edge_typ
    FlatEdge()
    {
      return kFlatEdge;
    }

    e_edge_typ
    get_max_polar_edge() const
    {
      return max_polar_edge;
    }

    e_edge_typ
    get_min_polar_edge() const
    {
      return min_polar_edge;
    }

    void
    set_max_polar_edge(e_edge_typ maxPolarEdge)
    {
      max_polar_edge = maxPolarEdge;
    }

    void
    set_min_polar_edge(e_edge_typ minPolarEdge)
    {
      min_polar_edge = minPolarEdge;
    }

   private:
    //! number of sectors
    int N_Sector;

    //! polar angle for the normal vector
    double normal_polar_angle;

    //! polar angle for edges
    double min_polar_angle;

    //! edge type
    e_edge_typ min_polar_edge;

    //! polar angle for edges
    double max_polar_angle;

    //! edge type
    e_edge_typ max_polar_edge;

    //! distance that detector starts from the normal direction
    double normal_start;

    //! base material, usually the gas. Will fill between layers
    std::string material;
  };

  //! \brief Generalized detector which use sectors of flat panels to cover full azimuthal acceptance
  class PHG4SectorConstructor
  {
   public:
    PHG4SectorConstructor(const std::string &name, PHG4Subsystem *subsys);
    virtual ~PHG4SectorConstructor() {}

    void
    Construct_Sectors(G4LogicalVolume *WorldLog);

    void
    OverlapCheck(bool check)
    {
      overlapcheck_sector = check;
    }

    void Verbosity(int v) { m_Verbosity = v; }
    int Verbosity() const { return m_Verbosity; }

   protected:
    bool overlapcheck_sector;

    G4VSolid *
    Construct_Sectors_Plane(           //
        const std::string &name,       //
        const double start_z,          //
        const double thickness,        //
        G4VSolid *SecConeBoundary_Det  //
    );

   public:
    // properties

    std::string name_base;

    Sector_Geometry geom;

   private:
    PHG4SectorDisplayAction *m_DisplayAction;
    int m_Verbosity;

   protected:
    G4LogicalVolume *
    RegisterLogicalVolume(G4LogicalVolume *);

    typedef std::map<G4String, G4LogicalVolume *> map_log_vol_t;
    map_log_vol_t map_log_vol;

    G4PVPlacement *
    RegisterPhysicalVolume(G4PVPlacement *v, const bool active = false);

    typedef std::pair<G4String, G4int> phy_vol_idx_t;
    typedef std::map<phy_vol_idx_t, G4PVPlacement *> map_phy_vol_t;
    map_phy_vol_t map_phy_vol;         //! all physics volume
    map_phy_vol_t map_active_phy_vol;  //! active physics volume
  };

}  // namespace PHG4Sector
#endif /* PHG4SectorConstructor_H_ */
