// $$Id: ePHENIXRICHConstruction.h,v 1.7 2014/05/01 19:02:45 phnxbld Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov> and Nils Feege <nils.feege@stonybrook.edu>
 * \version $$Revision: 1.7 $$
 * \date $$Date: 2014/05/01 19:02:45 $$
 */

#ifndef EPHENIXRICHCONSTRUCTION_H_
#define EPHENIXRICHCONSTRUCTION_H_

#include <Geant4/G4String.hh>

#ifndef __CINT__

#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4VUserDetectorConstruction.hh>

#include <map>
#include <set>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4PVPlacement;
#endif

class G4OpticalSurface;

namespace ePHENIXRICH
{
/**
   * \brief This class provides geometry and material parameters to construct the ePHENIX RICH
   * detector in ePHENIXRICHConstruction.
   *
   * \see ePHENIXRICH::RICH_Geometry
   * \see ePHENIXRICH::ePHENIXRICHConstruction
   * \see PHG4RICHDetector
   * \see PHG4RICHSteppingAction
   * \see PHG4RICHSubsystem
   *
   * \b Note: use GEANT units!
   */
class RICH_Geometry
{
 public:
  /**
     * Default constructor.
     */
  RICH_Geometry()
  {
    SetDefault();
    CreateOpticalSurfaces();
  }

  /**
     * Set default parameters for the RICH geometry.
     */
  void
  SetDefault();

  /**
     * Create optical surfaces (mirror, photocathode) to properly propagate optical photons in the RICH.
     */
  void
  CreateOpticalSurfaces();

  /**
     * Unit.
     */
  static double Unit_cm() { return cm; }
  /** @name Get Parameters
     *  Group of functions to get parameters. 
     */
  ///@{
  double
  get_R_max() const
  {
    return R_mirror_ref + dR_mirror + dR_mirror_spt + dR_backwindow;
  }

  double
  get_cone_size_z() const
  {
    return (z_shift + get_R_max()) * 2;
  }

  double
  get_R_frontwindow() const;

  double
  get_half_angle_HBD() const;

  double
  get_RZ_Seg1_HBD() const;

  double
  get_RZ_Seg2_HBD() const;

  double
  get_R_Tip_HBD() const;

  double
  get_Z_Tip_HBD() const;

  double
  get_Rotation_HBD() const;

  double
  get_R_beam_pipe() const
  {
    return R_beam_pipe;
  }

  double
  get_frontwindow_DisplaceRatio() const
  {
    return frontwindow_DisplaceRatio;
  }

  double
  get_min_eta() const
  {
    return min_eta;
  }

  double
  get_R_mirror_ref() const
  {
    return R_mirror_ref;
  }

  double
  get_dR_backwindow() const
  {
    return dR_backwindow;
  }

  double
  get_dR_frontwindow() const
  {
    return dR_frontwindow;
  }

  double
  get_dR_frontwindow_shrink() const
  {
    return dR_frontwindow_shrink;
  }

  double
  get_dR_mirror() const
  {
    return dR_mirror;
  }

  double
  get_dR_mirror_spt() const
  {
    return dR_mirror_spt;
  }

  G4String
  get_RICH_gas_mat() const
  {
    return RICH_gas_mat;
  }

  G4String
  get_RICH_Gas_Window_mat() const
  {
    return RICH_Gas_Window_mat;
  }

  G4String
  get_RICH_Mirror_mat() const
  {
    return RICH_Mirror_mat;
  }

  G4OpticalSurface*
  get_RICH_Mirror_OpticalSurface() const
  {
    return RICH_Mirror_OpticalSurface;
  }

  int get_N_RICH_Sector() const
  {
    return N_RICH_Sector;
  }

  double
  get_z_shift() const
  {
    return z_shift;
  }

  double
  get_R_shift() const
  {
    return R_shift;
  }

  int get_n_GEM_layers() const
  {
    return n_GEM_layers;
  }

  double
  get_HBD_thickness() const
  {
    return HBD_thickness;
  }

  G4OpticalSurface*
  get_RICH_Photocathode_OpticalSurface() const
  {
    return RICH_Photocathode_OpticalSurface;
  }
  ///@}

  /** @name Set Parameters
     *  Group of functions to set parameters. 
     */
  ///@{
  void
  set_R_beam_pipe(double beamPipe)
  {
    R_beam_pipe = beamPipe;
  }

  void
  set_frontwindow_DisplaceRatio(double frontwindowDisplaceRatio)
  {
    frontwindow_DisplaceRatio = frontwindowDisplaceRatio;
  }

  void
  set_min_eta(double minEta)
  {
    min_eta = minEta;
  }

  void
  set_R_mirror_ref(double mirrorRef)
  {
    R_mirror_ref = mirrorRef;
  }

  void
  set_dR_backwindow(double rBackwindow)
  {
    dR_backwindow = rBackwindow;
  }

  void
  set_dR_frontwindow(double rFrontwindow)
  {
    dR_frontwindow = rFrontwindow;
  }

  void
  set_dR_frontwindow_shrink(double rFrontwindowShrink)
  {
    dR_frontwindow_shrink = rFrontwindowShrink;
  }

  void
  set_dR_mirror(double rMirror)
  {
    dR_mirror = rMirror;
  }

  void
  set_dR_mirror_spt(double rMirrorSpt)
  {
    dR_mirror_spt = rMirrorSpt;
  }

  void
  set_RICH_gas_mat(G4String richGasMat)
  {
    RICH_gas_mat = richGasMat;
  }

  void
  set_RICH_Gas_Window_mat(G4String richGasWindowMat)
  {
    RICH_Gas_Window_mat = richGasWindowMat;
  }

  void
  set_RICH_Mirror_mat(G4String richMirrorMat)
  {
    RICH_Mirror_mat = richMirrorMat;
  }

  void
  set_N_RICH_Sector(int richSector)
  {
    N_RICH_Sector = richSector;
  }

  void
  set_z_shift(double shift)
  {
    z_shift = shift;
  }

  void
  set_R_shift(double shift)
  {
    R_shift = shift;
  }

  void
  set_n_GEM_layers(int gemLayers)
  {
    n_GEM_layers = gemLayers;
  }

  void
  set_HBD_thickness(double hbdThickness)
  {
    HBD_thickness = hbdThickness;
  }
  ///@}

 private:
  int N_RICH_Sector;

  double min_eta;
  double R_beam_pipe;

  double z_shift;
  double R_shift;
  double frontwindow_DisplaceRatio;  // Displace R,Z and radius simultainously
  double dR_frontwindow_shrink;
  double R_mirror_ref;

  double dR_mirror;
  double dR_mirror_spt;
  double dR_backwindow;
  double dR_frontwindow;

  int n_GEM_layers;
  double HBD_thickness;

  G4String RICH_gas_mat;
  G4String RICH_Mirror_mat;
  G4String RICH_Gas_Window_mat;

  G4OpticalSurface* RICH_Mirror_OpticalSurface;
  G4OpticalSurface* RICH_Photocathode_OpticalSurface;
};

#ifndef __CINT__

/**
   * \brief This class creates the ePHENIX RICH volumes for Geant4 based on the geometry
   * information in ePHENIXRICH::RICH_Geometry.
   *
   * \see ePHENIXRICH::RICH_Geometry
   * \see ePHENIXRICH::ePHENIXRICHConstruction
   * \see PHG4RICHDetector
   * \see PHG4RICHSteppingAction
   * \see PHG4RICHSubsystem
   *
   */
class ePHENIXRICHConstruction
{
 public:
  virtual ~ePHENIXRICHConstruction() {}
  ePHENIXRICHConstruction();
  ePHENIXRICHConstruction(const RICH_Geometry& g);

  virtual void
  OverlapCheck(bool check)
  {
    overlapcheck_rich = check;
  }

  G4LogicalVolume*
  Construct_RICH(G4LogicalVolume*);

  G4LogicalVolume*
  Construct_HBD(G4LogicalVolume* RICHSecLog);

  G4LogicalVolume*
  Construct_HBD_Layers(G4LogicalVolume* RICHHBDLog, const G4String name,
                       const G4String material, const double start_z, const double thickness);

  RICH_Geometry geom;

  /**
     * Checks if volume is sector volume
     */
  int is_in_sector(G4VPhysicalVolume*) const;

 protected:
  G4LogicalVolume*
  RegisterLogicalVolume(G4LogicalVolume*);

  typedef std::map<G4String, G4LogicalVolume*> map_log_vol_t;
  map_log_vol_t map_log_vol;

  G4PVPlacement*
  RegisterPhysicalVolume(G4PVPlacement*);

  typedef std::pair<G4String, G4int> phy_vol_idx_t;
  typedef std::map<phy_vol_idx_t, G4PVPlacement*> map_phy_vol_t;
  map_phy_vol_t map_phy_vol;

  bool overlapcheck_rich;

  /**
     * Set of volumes
     */
  std::set<G4VPhysicalVolume*> sector_vec;
};

#endif

}  //namespace ePHENIXRICH

#endif /* EPHENIXRICHCONSTRUCTION_H_ */
