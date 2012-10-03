FACETDSECT_lit.m is a "sectorized" lattice. All inputs to the deck are passed via the global PARM and LINAC structs. The ENLD and LEFF values in des_amp_and_phase.m come from LinacEnergy_20120803.txt. This file list all klystrons in the Linac. It was provided by MDW and was last updated August 2012.

To run, copy FACETDSECT_lit.m to the Beamlines directory and then execute run_designLattice.

Selecting the "uniform" lattice yields bunches that have sigma < 20 um and peak current ~ 20 kA for the parameters in run_designLattice.m. Selecting the "decker" lattice yields bunch lengths that have sigma ~ 78 um and peak current ~ 5 kA for the parameters in  run_designLattice.m.