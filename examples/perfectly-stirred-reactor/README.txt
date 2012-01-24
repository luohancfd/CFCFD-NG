Brendan's notes
---------------
The h2-chemkin, drm19, drm22, grimech30, jazbec_etal, li_and_williams 
tested and confirmed working with the following caveats:

- h2-chemkin (and other mechanisms) *must* use GRI-Mech3.0 thermodynamic data to match with
  Chemkin-II output. This data can be generated using kinetic_mechanisms/chemkin2libgas.py,
  but it's included as h2-chemkin-thermo-data.tar.gz for this case.
- jazbec_etal does not match exactly with the Chemkin-II output, even with GRI-Mech thermodynamic
  data. This is not an error in lib/gas, the finite-rate chemistry module nor chemkin2libgas
  as these have all been verified. Jazbec_etal also does nothave PD reactions. Thus, the
  difference must be in the input. Right now, don't have time to redo the Chemkin simulation.
- the li_and_williams chemkin output has some bad data, probably due to my poor fortran
  programming skills. Also the good data differs slightly. This mechanism does have PD reactions.
