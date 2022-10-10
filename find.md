/direct/sphenix+u/dstewart/coresoftware

Locations:

    simulation/g4simulation/g4tpc/PHG4TpcElectronDrift.h
    simulation/g4simulation/g4tpc/PHG4TpcElectronDrift.cc
    offline/packages/trackbase/TrkrHitSetContainer.h

    /direct/sphenix+u/dstewart/coresoftware/offline/packages/trackbase/TrkrClusterv4.cc
    /direct/sphenix+u/dstewart/coresoftware/offline/packages/trackbase/TrkrClusterv4.h
    /direct/sphenix+u/dstewart/coresoftware/simulation/g4simulation/g4tpc/PHG4TpcPadPlaneReadout.cc

Implementation of the TrkrClusterv4:

Good to know:
    gf - Edit existing file under cursor in same window
    C-Wf - Edit existing file under cursor in split window
    C-WC-F - Edit existing file under cursor in split window
    C-Wgf - Edit existing file under cursor in new tabpage


# Notes:
    (a) The single_hitsetcontainer has a key (contians information about the layers), and the energy 
    (neffelectrons) from ./simulation/g4simulation/g4tpc/PHG4TpcPadPlaneReadout.cc :: 310
    But -- I also need the Zbin (got it from the key), the pad-bin (phi-bin?!?), and the range...

    It also called populate_zigzag_phibins, which has the delta-phi information, I believe

    So... : ~/coresoftware/simulation/g4simulation/g4tpc/PHG4TpcPadPlaneReadout.cc :: 370-376

    In: ~/coresoftware/simulation/g4simulation/g4tpc/PHG4TpcPadPlaneReadout.cc
    void PHG4TpcPadPlaneReadout { 
        :: 129
            |
            :: 200 populate_zigzag_phibins( ... ) --> Get high and low phi bin
            | 
            :: 220 populate_tbins ( ... )  -->
            |
        :: 351
```
  const double philim_low = phi - (_nsigmas * cloud_sig_rp / radius) - phistepsize;
  const double philim_high = phi + (_nsigmas * cloud_sig_rp / radius) + phistepsize;

  // Find the pad range that covers this phi range
  int phibin_low  = LayerGeom->get_phibin(philim_low);
  int phibin_high = LayerGeom->get_phibin(philim_high);
  int npads = phibin_high - phibin_low; // note: pads = phibins!
```
kjk
